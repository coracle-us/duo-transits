"""Step 1+2: short-cadence coverage per host and flare detection.

Writes (in sidequests/flares/):
  coverage.csv        - one row per host (all 85): 2-min/20-s sectors, cameras, observing time
  sectors.csv         - one row per host-sector with 2-min data: noise, obs time, window used
  candidates_raw.csv  - every event passing N1/N2/N3, with automatic flags
  cache/<tic>_s<sec>.npz - detrended residuals (used by the cross-star check)
"""
import numpy as np
import pandas as pd
from tess_stars2px import tess_stars2px_function_entry

import flarelib as F


def cameras(row) -> dict:
    out = tess_stars2px_function_entry(int(row.TIC), float(row.RA), float(row.Dec))
    return {int(s): int(c) for s, c in zip(out[3], out[4])}


def main():
    cat = F.load_catalog()
    F.CACHE.mkdir(exist_ok=True)
    cov_rows, sec_rows, cands = [], [], []
    for tic, row in cat.iterrows():
        f120 = F.spoc_files(tic, 120)
        f20 = F.spoc_files(tic, 20)
        cams = cameras(row)
        obs_tot, raw_tot = 0.0, 0.0
        for sec in f120:
            lc = F.detrend(F.load_sector(tic, sec, cat))
            obs = F.obs_time_days(lc)
            raw = float(lc["good"].sum() * F.CADENCE_D)
            obs_tot += obs
            raw_tot += raw
            sec_rows.append(dict(tic=tic, sector=sec, camera=cams.get(sec, -1), has_20s=sec in f20,
                                 good_days=round(raw, 3), obs_days=round(obs, 3),
                                 sigma_ppm=round(float(np.nanmedian(lc["sigma"][lc["good"]])) * 1e6, 1),
                                 n_segments=len(lc["segs"]),
                                 win_h=";".join(f"{w*24:.0f}" for w in lc["win_used"]),
                                 n_dumps=len(lc["dumps"])))
            np.savez_compressed(F.CACHE / f"{tic}_s{sec:03d}.npz", t=lc["t"], resid=lc["resid"],
                                sigma=lc["sigma"], good=lc["good"])
            fl = F.auto_flags(F.find_flares(lc), lc)
            if len(fl):
                fl["camera"] = cams.get(sec, -1)
                cands.append(fl)
        cov_rows.append(dict(TIC=tic, Teff=row.Teff, Rstar=row.Rstar, Tmag=row.Tmag,
                             has_2min=len(f120) > 0, has_20s=len(f20) > 0,
                             sectors_2min=" ".join(str(s) for s in f120),
                             sectors_20s=" ".join(str(s) for s in f20),
                             n_sectors_2min=len(f120), n_sectors_20s=len(f20),
                             good_days_2min=round(raw_tot, 2), obs_days_2min=round(obs_tot, 2),
                             note="" if f120 else "no SPOC 2-min/20-s data on disk; skipped"))
        print(tic, len(f120), "sectors", f"{obs_tot:.1f} d",
              sum(len(c) for c in cands if len(c) and c.tic.iloc[0] == tic), "cands")
    pd.DataFrame(cov_rows).to_csv(F.SQ / "coverage.csv", index=False)
    pd.DataFrame(sec_rows).to_csv(F.SQ / "sectors.csv", index=False)
    cand = pd.concat(cands, ignore_index=True) if cands else pd.DataFrame()
    cand.insert(0, "cand_id", [f"C{k:03d}" for k in range(len(cand))])
    cand.to_csv(F.SQ / "candidates_raw.csv", index=False)
    print("hosts with 2-min:", sum(r["has_2min"] for r in cov_rows), "candidates:", len(cand))


if __name__ == "__main__":
    main()
