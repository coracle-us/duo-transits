"""
Download the light curves, cut the windows and compute features (RULES.md section 2 and amendment A).

Input:  data/windows.csv, data/sample.csv (TOIs), data/eb_windows.csv, data/eb_sample.csv (EBs), data/toi_2026-11-24.csv
Cache:  ~/Documents/GitHub/coracle/duo-transits-cache/blend-ceiling/ (MAST FITS and binned npz; not in git)
Output: data/features.csv  one row per valid window: obj_id, tic, kind (pos, fpfa, eb, tce, sanity), label, sector,
                           author, t_c, the features in lcwin.FEATURES, and the 200-point window (for the CNN) in
        data/windows_flux.npz
        data/window_log.csv  every attempted window and why it was dropped
Usage:  .venv/bin/python -I make_windows.py [--workers 4]
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir; our modules live here
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

import lcwin

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
CACHE = os.path.expanduser("~/Documents/GitHub/coracle/duo-transits-cache/blend-ceiling")
SEED = 20261124
_LOG = open(os.path.join(CACHE, "progress.log"), "a") if os.path.isdir(CACHE) else None


def say(*a):
    """Progress messages go to our own log handle: something in the download stack closes sys.stdout."""
    msg = " ".join(str(x) for x in a)
    if _LOG:
        _LOG.write(msg + "\n")
        _LOG.flush()
    try:
        print(msg, flush=True)
    except ValueError:
        pass


def npz_path(tic, sector, author):
    return os.path.join(CACHE, "binned", f"{int(tic)}_{int(sector):03d}_{author}.npz")


def fetch(job):
    tic, sector, author, exptime = job
    p = npz_path(tic, sector, author)
    if os.path.exists(p):
        return job, "cached"
    import lightkurve as lk
    err = ""
    for attempt in range(3):
        try:
            sr = lk.search_lightcurve(f"TIC {tic}", sector=int(sector), author=author)
            sr = sr[np.asarray(sr.exptime.value) == exptime] if len(sr) > 1 else sr
            if len(sr) == 0:
                return job, "not found"
            lc = sr[0].download(download_dir=os.path.join(CACHE, "mast"))
            t, f = lcwin.prep(lc)
            os.makedirs(os.path.dirname(p), exist_ok=True)
            np.savez_compressed(p, t=t, f=f)
            return job, "ok"
        except Exception as e:
            err = str(e)[:150]
    return job, "error: " + err


def load(tic, sector, author):
    p = npz_path(tic, sector, author)
    if not os.path.exists(p):
        return None, None
    z = np.load(p)
    return z["t"], z["f"]


def exclusion_times(tic, toi, ebs, lo, hi):
    """Predicted transit or eclipse times on a star inside [lo-5, hi+5] d: every TOI (any disposition) and the EB
    primary and secondary eclipses."""
    out = []
    for r in toi[toi.tic == tic].itertuples():
        E, P = r.epoch - 2457000.0, r.period
        if not np.isfinite(P) or P == 0:
            out.append(E)
        else:
            out += list(E + P * np.arange(np.ceil((lo - 5 - E) / P), np.floor((hi + 5 - E) / P) + 1))
    for r in ebs[ebs.tic == tic].itertuples():
        E, P = r.epoch_btjd, r.period_d
        if not np.isfinite(P) or P == 0:
            out.append(E)
            continue
        ph2 = r.phis if np.isfinite(r.phis) else 0.5
        for off in (0.0, ph2 * P):
            out += list(E + off + P * np.arange(np.ceil((lo - 5 - E - off) / P), np.floor((hi + 5 - E - off) / P) + 1))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    s = pd.read_csv(os.path.join(D, "sample.csv"))
    s = s[s.drop_reason.isna()]
    eb = pd.read_csv(os.path.join(D, "eb_sample.csv"))
    eb = eb[eb.drop_reason.isna()]
    prsa = pd.read_csv(os.path.join(HERE, "..", "..", "sidequests", "eb_timing", "data", "prsa2022_tess_ebs.csv"))
    phis = prsa.drop_duplicates("TIC").set_index("TIC")["Phis-pf"]
    eb = eb.assign(phis=eb.tic.map(phis))
    toi = pd.read_csv(os.path.join(D, "toi_2026-11-24.csv"), usecols=["TIC ID", "Epoch (BJD)", "Period (days)"])
    toi = toi.rename(columns={"TIC ID": "tic", "Epoch (BJD)": "epoch", "Period (days)": "period"})
    w = pd.read_csv(os.path.join(D, "windows.csv"))
    w = w[w.toi.isin(s.toi)].merge(s[["toi", "disposition", "duration_h"]], on="toi")
    w["kind"] = np.where(w.label == 1, "pos", "fpfa")
    w["obj_id"] = "TOI-" + w.toi.astype(str)
    ew = pd.read_csv(os.path.join(D, "eb_windows.csv"))
    ew = ew[ew.tic.isin(eb.tic)].merge(eb[["tic", "duration_h"]], on="tic")
    ew["kind"], ew["obj_id"], ew["disposition"] = "eb", ew.toi, "EB"
    allw = pd.concat([w, ew], ignore_index=True)

    jobs = sorted({(int(r.tic), int(r.sector), r.author, float(r.exptime)) for r in allw.itertuples()})
    say(f"light curves needed: {len(jobs)}")
    status = {}
    with ThreadPoolExecutor(a.workers) as ex:
        for i, (job, st) in enumerate(ex.map(fetch, jobs), 1):
            status[job] = st
            if i % 100 == 0:
                say(f"  fetched {i}/{len(jobs)}")
    bad = {k: v for k, v in status.items() if v not in ("ok", "cached")}
    say(f"light curves: {len(jobs) - len(bad)} ok, {len(bad)} failed")

    rows, fluxes, log = [], [], []

    def add(base, t, f, T, t14):
        if t is None:
            log.append(dict(base, t_c=T, status="no light curve"))
            return
        if not lcwin.valid(t, T, t14):
            log.append(dict(base, t_c=T, status="invalid window"))
            return
        x, filled = lcwin.window(t, f, T)
        ft = None if x is None else lcwin.features(x, filled)
        if ft is None:
            log.append(dict(base, t_c=T, status="no features"))
            return
        rows.append(dict(base, t_c=T, **ft))
        fluxes.append(np.where(filled, 0.0, x).astype(np.float32))
        log.append(dict(base, t_c=T, status="ok"))

    for r in allw.itertuples():
        t, f = load(r.tic, r.sector, r.author)
        add(dict(obj_id=r.obj_id, tic=r.tic, kind=r.kind, label=int(r.kind == "pos"), disposition=r.disposition,
                 sector=r.sector, author=r.author), t, f, r.t_pred_btjd, r.duration_h)

    # TCE-like windows (A2): positive and EB stars, highest-SNR box per light curve, up to 2 per star
    rng = np.random.default_rng(SEED)
    lcs = allw[allw.kind.isin(["pos", "eb"])][["tic", "sector", "author"]].drop_duplicates()
    for tic, g in lcs.groupby("tic"):
        cands = []
        for r in g.itertuples():
            t, f = load(r.tic, r.sector, r.author)
            if t is None or len(t) == 0:
                continue
            ex = exclusion_times(tic, toi, eb, t.min(), t.max())
            b = lcwin.sector_box_search(t, f, ex)
            if b:
                cands.append((b["snr"], r.sector, r.author, b["t_c"], b["dur_h"]))
        cands.sort(reverse=True)
        pick = [c for i, c in enumerate(cands) if c[1] not in {x[1] for x in cands[:i]}][:2]
        if len(pick) < 2:
            pick += [c for c in cands if c not in pick][: 2 - len(pick)]
        for snr, sec, au, tc, dur in pick:
            t, f = load(tic, sec, au)
            add(dict(obj_id=f"TCE{tic}_{sec}", tic=tic, kind="tce", label=0, disposition="TCE-like", sector=sec,
                     author=au), t, f, tc, dur)

    # sanity windows: 50 random out-of-transit windows from the positives' light curves (RULES section 4)
    pl = allw[allw.kind == "pos"][["tic", "sector", "author", "duration_h"]].drop_duplicates(["tic", "sector", "author"])
    n, tries = 0, 0
    while n < 50 and tries < 5000:
        tries += 1
        r = pl.iloc[rng.integers(len(pl))]
        t, f = load(r.tic, r.sector, r.author)
        if t is None or len(t) < lcwin.L:
            continue
        T = rng.uniform(t.min() + 2, t.max() - 2)
        ex = exclusion_times(r.tic, toi, eb, t.min(), t.max())
        if len(ex) and np.min(np.abs(ex - T)) < 3 or not lcwin.valid(t, T, r.duration_h):
            continue
        before = len(rows)
        add(dict(obj_id=f"SAN{n}", tic=int(r.tic), kind="sanity", label=0, disposition="sanity",
                 sector=int(r.sector), author=r.author), t, f, T, r.duration_h)
        n += len(rows) > before
    say(f"sanity windows: {n} in {tries} draws")

    F = pd.DataFrame(rows)
    F.to_csv(os.path.join(D, "features.csv"), index=False)
    np.savez_compressed(os.path.join(D, "windows_flux.npz"), x=np.stack(fluxes), t_rel_h=lcwin.T_REL_H)
    L = pd.DataFrame(log)
    L.to_csv(os.path.join(D, "window_log.csv"), index=False)
    say(L.groupby(["kind", "status"]).size().to_string())
    say("objects with >= 1 valid window:", F.groupby("kind").obj_id.nunique().to_dict())


if __name__ == "__main__":
    main()
