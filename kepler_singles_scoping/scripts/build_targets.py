"""
Merge the per-source rows (data/events_by_source.csv) into unique planet candidates, cross-match KIC -> TIC,
record Kepler quarters observed (lightkurve/MAST) and TESS sectors (tess-point).

Output: data/kepler_singles_targets.csv
Run:    .venv/Scripts/python -I kepler_singles_scoping/scripts/build_targets.py
"""
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent
DUR_PRIORITY = ["FM16", "U16", "W15", "KM19", "HZW19", "HC24"]


def merge_objects(d):
    """Group rows that share a transit time (|dt| < 1 d) on the same KIC into one candidate."""
    d = d.copy()
    d["t2_bkjd"] = np.where(d.n_kepler_transits == 2, d.t0_bkjd + d.P_cat, np.nan)
    objs = []
    for kic, g in d.groupby("KIC"):
        groups = []
        for i, r in g.iterrows():
            times = {round(r.t0_bkjd, 1)} | ({round(r.t2_bkjd, 1)} if np.isfinite(r.t2_bkjd) else set())
            hit = None
            for grp in groups:
                if any(abs(a - b) < 1.0 for a in times for b in grp["times"]):
                    hit = grp
                    break
            if hit is None:
                groups.append(dict(times=set(times), rows=[i]))
            else:
                hit["times"] |= times
                hit["rows"].append(i)
        for grp in groups:
            rows = d.loc[grp["rows"]]
            duo = rows[rows.n_kepler_transits == 2]
            ts = sorted(grp["times"])
            # collapse near-duplicate times
            tt = []
            for t in ts:
                if not tt or abs(t - tt[-1]) > 1.0:
                    tt.append(t)
            if len(duo):
                t1 = float(duo.t0_bkjd.iloc[0])
                dT = float(duo.P_cat.iloc[0])
                # a W15 'P' for a 2-transit system may be dT/2 (e.g. 10525077 427 vs 854): keep the largest
                dT = float(duo.P_cat.max())
                ntr = 2
            elif len(tt) >= 2:  # two separate single rows (e.g. U16 847 events) -> a duo
                t1, dT, ntr = tt[0], tt[-1] - tt[0], 2
            else:
                t1, dT, ntr = float(rows.t0_bkjd.median()), np.nan, 1
            durrow = None
            for s in DUR_PRIORITY:
                x = rows[(rows.source == s) & rows.dur_h.notna()]
                if len(x):
                    durrow = x.iloc[0]
                    break
            Pcat = rows[rows.n_kepler_transits == 1].P_cat.median()
            objs.append(dict(KIC=kic, sources=" ".join(sorted(set(rows.source))), n_sources=rows.source.nunique(),
                             KOI=" ".join(sorted({str(x) for x in rows.KOI if isinstance(x, str) and x.strip()})),
                             t1_bkjd=t1, dT_duo=dT, n_kepler_transits=ntr,
                             depth_ppm=float(rows.depth_ppm.median()),
                             dur_h=float(durrow.dur_h) if durrow is not None else np.nan,
                             dur_src=durrow.source if durrow is not None else "",
                             Rp_re=float(rows.Rp_re.median()), Rstar=float(rows.Rstar.median()),
                             Kp_cat=float(rows.Kp.median()), P_cat_single=Pcat,
                             km19_table=";".join(sorted(set(rows.table.dropna()) & {"clean", "flagged"})),
                             flags=" | ".join(f"{a}:{b}" for a, b in zip(rows.source, rows["flags"].fillna(""))
                                              if str(b).strip())))
    o = pd.DataFrame(objs)
    o["t1_btjd"] = o.t1_bkjd - 2167.0
    return o


def tic_by_kic(kics):
    from astroquery.mast import Catalogs
    out = []
    kics = list(kics)
    for i in range(0, len(kics), 50):
        for attempt in range(3):
            try:
                t = Catalogs.query_criteria(catalog="Tic", KIC=kics[i:i + 50]).to_pandas()
                break
            except Exception:  # noqa: BLE001
                time.sleep(5)
        out.append(t[["ID", "KIC", "ra", "dec", "Tmag", "Teff", "rad", "GAIA", "disposition"]])
    t = pd.concat(out)
    t["_bad"] = ~(t.disposition.isna() | (t.disposition == ""))  # prefer non-artifact/duplicate TIC entries
    t = t.sort_values("_bad").drop(columns="_bad")
    t["KIC"] = t.KIC.astype(int)
    t = t.rename(columns={"ID": "TIC", "rad": "Rstar_tic", "Teff": "Teff_tic"})
    return t.drop_duplicates("KIC")


def kepler_quarters(kic):
    import lightkurve as lk
    for attempt in range(3):
        try:
            r = lk.search_lightcurve(f"KIC {kic}", author="Kepler", exptime=1800)
            qs = sorted({int(m.split()[-1]) for m in r.table["mission"] if "Quarter" in m})
            return " ".join(map(str, qs))
        except Exception:  # noqa: BLE001
            time.sleep(5)
    return ""


def tess_sectors(d):
    import tess_stars2px as tp
    ends = pd.read_csv(ROOT / "singles" / "data" / "tess_orbit_times_latest.csv", comment="#")
    ends["end"] = pd.to_datetime(ends["End of Orbit"], utc=True)
    sec_end = ends.groupby("Sector").end.max()
    today = pd.Timestamp(date.today(), tz="UTC")
    done = {int(s) for s, e in sec_end.items() if e < today}
    alls, obs = [], []
    for r in d.itertuples():
        o = tp.tess_stars2px_function_entry(int(r.TIC), r.ra, r.dec)
        s = sorted({int(x) for x in o[3] if x > 0})
        alls.append(" ".join(map(str, s)))
        obs.append(" ".join(str(x) for x in s if x in done))
    d["tess_sectors_all"] = alls
    d["tess_sectors_observed"] = obs
    d["n_tess_observed"] = [len(x.split()) for x in obs]
    return d


def main():
    ev = pd.read_csv(HERE / "data" / "events_by_source.csv")
    o = merge_objects(ev)
    print(f"{len(ev)} source rows -> {len(o)} unique candidates on {o.KIC.nunique()} stars")
    t = tic_by_kic(o.KIC.unique())
    o = o.merge(t, on="KIC", how="left")
    print(f"TIC matched by KIC column: {o.TIC.notna().sum()}/{len(o)}")
    cache = HERE / "data" / "kepler_quarters_cache.csv"
    qc = pd.read_csv(cache) if cache.exists() else pd.DataFrame(columns=["KIC", "quarters"])
    have = set(qc.KIC)
    new = []
    for k in o.KIC.unique():
        if k not in have:
            new.append(dict(KIC=k, quarters=kepler_quarters(k)))
    if new:
        qc = pd.concat([qc, pd.DataFrame(new)])
        qc.to_csv(cache, index=False)
    o = o.merge(qc, on="KIC", how="left")
    o = tess_sectors(o)
    o.to_csv(HERE / "data" / "kepler_singles_targets.csv", index=False)
    print(o.n_tess_observed.describe())
    print(o[["KIC", "TIC", "sources", "n_kepler_transits", "t1_bkjd", "dT_duo", "depth_ppm", "dur_h", "Tmag",
             "n_tess_observed"]].round(2).to_string())


if __name__ == "__main__":
    main()
