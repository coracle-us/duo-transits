"""
Easy negatives 1 (RULES.md amendment A2): long-period Prsa et al. (2022) EBs that are not TOIs, with the same
light-curve, coverage and timing rules as the TOI sample.

Input:  ../../sidequests/eb_timing/data/prsa2022_tess_ebs.csv (full VizieR J/ApJS/258/16 table)
        data/toi_2026-11-24.csv (to drop TOI hosts)
Output: data/eb_sample.csv, data/eb_windows.csv; MAST search results are added to data/lc_search.csv
Usage:  .venv/bin/python -I build_ebs.py [--workers 4]
"""
import argparse
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I drops the script dir; our modules live here

import numpy as np
import pandas as pd

from build_sample import HERE, ROOT, TOI, lc_search, orbits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    e = pd.read_csv(os.path.join(ROOT, "sidequests", "eb_timing", "data", "prsa2022_tess_ebs.csv"))
    e = e.drop_duplicates("TIC")
    tois = set(pd.read_csv(TOI, usecols=["TIC ID"])["TIC ID"])
    e = e[~e.TIC.isin(tois) & ((e.Per > 10) | e.Per.isna() | (e.Per == 0))].copy()
    e["single"] = e.Per.isna() | (e.Per == 0)
    e["w_phase"] = e["Wp-pf"].fillna(e["Wp-2g"])
    e["depth_ppm"] = 1e6 * e["Dp-pf"].fillna(e["Dp-2g"])
    e["duration_h"] = 24 * e.w_phase * e.Per
    srch = lc_search(e.TIC.unique(), a.workers)
    orb = orbits()
    s = srch[srch.sector > 0].copy()
    s["rank"] = np.where((s.author == "SPOC") & (s.exptime == 120), 0, np.where(s.author == "QLP", 1, 9))
    s = s[s["rank"] < 9].sort_values(["tic", "sector", "rank"]).drop_duplicates(["tic", "sector"])
    by_tic = {t: g for t, g in s.groupby("tic")}
    failed = set(srch.loc[srch.sector == -2, "tic"])
    rows, wins = [], []
    for r in e.itertuples():
        rec = dict(tic=r.TIC, eb_id=f"EB{r.TIC}", label=0, period_d=r.Per, epoch_btjd=r.BJD0,
                   depth_ppm=r.depth_ppm, duration_h=r.duration_h, tmag=r.Tmag, morph=r.Morph)
        if not np.isfinite(r.duration_h):
            rows.append(dict(rec, n_windows=0, sectors_used="", drop_reason="no primary width in the catalogue"))
            continue
        if r.TIC in failed or r.TIC not in by_tic:
            rows.append(dict(rec, n_windows=0, sectors_used="",
                             drop_reason="MAST search failed" if r.TIC in failed else "no SPOC 2-min or QLP light curve"))
            continue
        sE = r.e_BJD0 if np.isfinite(r.e_BJD0) else 0.0
        sP = r.e_Per if np.isfinite(r.e_Per) else 0.0
        n_any = 0
        used = set()
        for sec in by_tic[r.TIC].itertuples():
            ob = orb[orb.sector == sec.sector]
            if ob.empty:
                continue
            lo, hi = ob.start.min(), ob.end.max()
            ks = [0] if r.single else range(int(np.ceil((lo - r.BJD0) / r.Per)), int(np.floor((hi - r.BJD0) / r.Per)) + 1)
            for k in ks:
                T = r.BJD0 + (0 if r.single else k * r.Per)
                if r.single and not (lo <= T <= hi):
                    continue
                if not ((ob.start <= T) & (T <= ob.end)).any():
                    continue
                n_any += 1
                sig_h = 24 * np.hypot(sE, k * sP)
                if 3 * sig_h >= 25:
                    continue
                used.add(sec.sector)
                wins.append(dict(tic=r.TIC, toi=f"EB{r.TIC}", label=0, sector=sec.sector, author=sec.author,
                                 exptime=sec.exptime, k=k, t_pred_btjd=T, sigma_t_h=sig_h))
        n = sum(1 for w in wins if w["tic"] == r.TIC)
        reason = "" if n else ("no predicted eclipse inside a sector with a light curve" if n_any == 0
                               else "all predicted eclipses have 3 sigma_T >= 25 h")
        rows.append(dict(rec, n_windows=n, sectors_used=" ".join(map(str, sorted(used))), drop_reason=reason))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(HERE, "data", "eb_sample.csv"), index=False)
    pd.DataFrame(wins).to_csv(os.path.join(HERE, "data", "eb_windows.csv"), index=False)
    k = out[out.drop_reason == ""]
    print(f"EBs: {len(out)} long-period non-TOI, kept {len(k)} with {int(k.n_windows.sum())} windows")
    print(out[out.drop_reason != ""].drop_reason.value_counts().to_string())


if __name__ == "__main__":
    main()
