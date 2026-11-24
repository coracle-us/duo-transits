"""
Step 2 of the blend-ceiling test: the labelled long-period TOI sample (rules in RULES.md, section 1).

Input:  data/toi_2026-11-24.csv (ExoFOP TOI table, unedited)
        ../../singles/data/tess_orbit_times_latest.csv (TESS orbit start/end, UTC)
Cache:  data/lc_search.csv (MAST light-curve search results per TIC; metadata only, no light curves)
Output: data/sample.csv (every TOI passing the label and period cuts, kept or with a drop_reason)
        data/windows.csv (each usable predicted transit time: TOI, sector, author, T, sigma_T)
Usage:  .venv/bin/python -I build_sample.py [--workers 4]
"""
import argparse
import os
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
TOI = os.path.join(HERE, "data", "toi_2026-11-24.csv")
ORBITS = os.path.join(ROOT, "singles", "data", "tess_orbit_times_latest.csv")
SEARCH = os.path.join(HERE, "data", "lc_search.csv")
UTC_TO_BTJD = 69.184 / 86400 - 2457000.0  # JD_UTC + (TT-UTC) - 2457000; barycentric offset ignored (RULES 1)


def load_toi():
    d = pd.read_csv(TOI)
    d["disp"] = d["TFOPWG Disposition"].fillna("").str.strip()
    P = d["Period (days)"]
    d["single"] = P.isna() | (P == 0)
    keep = d.disp.isin(["CP", "KP", "FP", "FA"]) & ((P > 10) | d.single)
    d = d[keep].copy()
    d["label"] = d.disp.isin(["CP", "KP"]).astype(int)
    d = d.rename(columns={"TIC ID": "tic", "Period (days)": "P", "Period (days) err": "P_err",
                          "Epoch (BJD) err": "E_err", "Depth (ppm)": "depth", "Duration (hours)": "dur",
                          "TESS Mag": "tmag"})
    d["epoch_btjd"] = d["Epoch (BJD)"] - 2457000.0
    return d


def search_one(tic):
    import lightkurve as lk
    for attempt in range(4):
        try:
            sr = lk.search_lightcurve(f"TIC {tic}", author=["SPOC", "QLP"])
            t = sr.table
            if len(t) == 0:
                return [dict(tic=tic, sector=-1, author="none", exptime=np.nan, target_name="")]
            return [dict(tic=tic, sector=int(r["mission"].split()[-1]), author=str(r["author"]),
                         exptime=float(r["exptime"]), target_name=str(r["target_name"])) for r in t]
        except Exception as e:  # MAST hiccups: retry, then record the failure
            err = str(e)
            time.sleep(5 * (attempt + 1))
    return [dict(tic=tic, sector=-2, author="error", exptime=np.nan, target_name=err[:200])]


def lc_search(tics, workers):
    done = pd.read_csv(SEARCH) if os.path.exists(SEARCH) else pd.DataFrame(columns=["tic"])
    done = done[done.sector != -2] if len(done) else done  # retry earlier failures
    todo = sorted(set(tics) - set(done.tic))
    print(f"MAST search: {len(todo)} TICs to query, {len(set(done.tic))} cached")
    rows = []
    with ThreadPoolExecutor(workers) as ex:
        for i, r in enumerate(ex.map(search_one, todo), 1):
            rows += r
            if i % 50 == 0:
                print(f"  {i}/{len(todo)}")
    out = pd.concat([done, pd.DataFrame(rows)], ignore_index=True) if rows else done
    out.to_csv(SEARCH, index=False)
    return out


def orbits():
    o = pd.read_csv(ORBITS, comment="#")
    jd = lambda s: pd.to_datetime(s).map(lambda x: x.to_julian_date()) + UTC_TO_BTJD
    return pd.DataFrame(dict(sector=o["Sector"], start=jd(o["Start of Orbit"]), end=jd(o["End of Orbit"])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    d = load_toi()
    srch = lc_search(d["tic"].unique(), a.workers)
    orb = orbits()
    # per TIC and sector: SPOC 120 s if present, else QLP (RULES 1)
    s = srch[srch.sector > 0].copy()
    s["rank"] = np.where((s.author == "SPOC") & (s.exptime == 120), 0, np.where(s.author == "QLP", 1, 9))
    s = s[s["rank"] < 9].sort_values(["tic", "sector", "rank"]).drop_duplicates(["tic", "sector"])
    by_tic = {t: g for t, g in s.groupby("tic")}
    failed = set(srch.loc[srch.sector == -2, "tic"])

    rows, wins = [], []
    for r in d.itertuples(index=False):
        tic = r.tic
        rec = dict(tic=tic, toi=r.TOI, disposition=r.disp, label=r.label,
                   period_d=r.P if not r.single else np.nan,
                   single=r.single, epoch_btjd=r.epoch_btjd, depth_ppm=r.depth, duration_h=r.dur,
                   tmag=r.tmag, ra=r.RA, dec=r.Dec)
        if tic in failed:
            rows.append(dict(rec, n_windows=0, sectors_used="", drop_reason="MAST search failed"))
            continue
        g = by_tic.get(tic)
        if g is None:
            rows.append(dict(rec, n_windows=0, sectors_used="", drop_reason="no SPOC 2-min or QLP light curve"))
            continue
        sE = r.E_err if np.isfinite(r.E_err) else 0.0
        sP = r.P_err if np.isfinite(r.P_err) else 0.0
        n_any, n_sig = 0, 0
        for sec in g.itertuples():
            ob = orb[orb.sector == sec.sector]
            if ob.empty:
                continue
            lo, hi = ob.start.min(), ob.end.max()
            if r.single:
                ks = [0] if lo <= r.epoch_btjd <= hi else []
            else:
                P = r.P
                ks = range(int(np.ceil((lo - r.epoch_btjd) / P)), int(np.floor((hi - r.epoch_btjd) / P)) + 1)
            for k in ks:
                T = r.epoch_btjd + k * (0 if r.single else r.P)
                if not ((ob.start <= T) & (T <= ob.end)).any():
                    continue  # in the gap between orbits
                n_any += 1
                sig_h = 24 * np.hypot(sE, k * sP)
                if 3 * sig_h >= 25:
                    continue
                n_sig += 1
                wins.append(dict(tic=tic, toi=r.TOI, label=r.label, sector=sec.sector, author=sec.author,
                                 exptime=sec.exptime, k=k, t_pred_btjd=T, sigma_t_h=sig_h))
        used = sorted({w["sector"] for w in wins if w["toi"] == r.TOI})
        reason = "" if n_sig else ("no predicted transit inside a sector with a light curve" if n_any == 0
                                   else "all predicted transits have 3 sigma_T >= 25 h")
        rows.append(dict(rec, n_windows=n_sig, sectors_used=" ".join(map(str, used)), drop_reason=reason))

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(HERE, "data", "sample.csv"), index=False)
    pd.DataFrame(wins).to_csv(os.path.join(HERE, "data", "windows.csv"), index=False)
    k = out[out.drop_reason == ""]
    print("\nbefore coverage:", out.disposition.value_counts().to_dict())
    print("kept:", k.disposition.value_counts().to_dict(), f"-> {int(k.label.sum())} positives, "
          f"{int((k.label == 0).sum())} negatives, {k.tic.nunique()} stars, {int(k.n_windows.sum())} windows")
    print("drop reasons:\n" + out[out.drop_reason != ""].groupby(["drop_reason", "disposition"]).size().to_string())


if __name__ == "__main__":
    main()
