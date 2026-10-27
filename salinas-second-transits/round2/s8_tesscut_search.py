#!/usr/bin/env python3
"""Round 2, section 8: TESScut photometry of a sector with no light curve product,
searched with the pipeline's own detrending, box search and cuts.

  python round2/s8_tesscut_search.py 419959752 97
"""
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import st2_pipeline as pl  # noqa: E402
import st2core as core  # noqa: E402
import tesscut_phot as tp  # noqa: E402

CUTS = SimpleNamespace(snr_min=7.0, depth_lo=0.5, depth_hi=2.0, dur_lo=0.5, dur_hi=2.0, anti_ratio=1.5)


def main(tic, sector):
    work = HERE / "work"
    out = HERE / "results" / "round2"
    out.mkdir(parents=True, exist_ok=True)
    cov = pd.read_csv(work / "coverage.csv").set_index("tic").loc[tic]
    rows = []
    for p in sorted((work / "per_cand").glob(f"{tic}*.json")):
        r = json.loads(p.read_text())
        try:
            t, f, e, info = tp.sector_lightcurve(tic, cov.ra, cov.dec, sector, work / "tesscut")
        except RuntimeError as err:
            rows.append(dict(cand_id=r["cand_id"], tic=tic, sector=sector, error=str(err)))
            continue
        cad = core.cadence_days(t)
        if cad < 9 / 1440:  # same 10-min binning as load_sector_lc
            edges = np.arange(t.min(), t.max() + 10 / 1440, 10 / 1440)
            idx = np.digitize(t, edges)
            ok = np.bincount(idx) > 0
            t, f = (np.bincount(idx, t) / np.bincount(idx))[ok], (np.bincount(idx, f) / np.bincount(idx))[ok]
        dur_ref, depth_ref = r["dur_ref_d"], r["depth_ref"]
        t, f = pl.detrend(t, f, dur_ref)
        np.savez_compressed(work / "lc" / f"{tic}_s{sector:03d}_tesscut.npz", t=t, f=f)
        cen, dep, snr, frac, sig = core.box_search(t, f, dur_ref)
        anti = -np.nanmin(snr)
        for i in core.pick_events(cen, snr, dur_ref, n_max=3):
            ev = core.refine_event(t, f, cen[i], dur_ref)
            if not ev:
                continue
            ev.update(sector=sector, author="TESScut", noise_ppm=sig * 1e6, anti_snr=float(anti),
                      depth_ratio=ev["depth"] / depth_ref, dur_ratio=ev["dur_d"] / dur_ref,
                      edge_dist_d=core.gap_edge_distance(t, ev["tc"]))
            rows.append(dict(cand_id=r["cand_id"], tic=tic, sector=sector, n_ap_pix=info["n_ap"],
                             target_to_detector_edge_px=info["edge_px"],
                             n_high_bkg_dropped=info["n_high_bkg_dropped"],
                             n_points=t.size, span_d=t.max() - t.min(), t2=ev["tc"],
                             depth_ppm=ev["depth"] * 1e6, depth_err_ppm=ev["depth_err"] * 1e6,
                             depth_ref_ppm=depth_ref * 1e6, dur_h=ev["dur_d"] * 24, dur_ref_h=dur_ref * 24,
                             snr=ev["snr"], anti_snr=anti, depth_ratio=ev["depth_ratio"],
                             dur_ratio=ev["dur_ratio"], edge_dist_d=ev["edge_dist_d"], frac=ev["frac"],
                             passes=pl.passes(ev, CUTS)))
        plot(t, f, rows, tic, sector, depth_ref, out, info["edge_px"])
    df = pd.DataFrame(rows)
    if "error" in df:
        print("\n".join(df["error"].dropna().unique()))
    df.to_csv(out / f"tic{tic}_s{sector}_tesscut_search.csv", index=False)
    print(df.to_string(index=False))


def plot(t, f, rows, tic, sector, depth_ref, out, edge_px):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(t, (f - 1) * 1e6, ".", ms=1, color="0.5")
    for r in rows:
        ax.axvline(r["t2"], color="C3" if r["passes"] else "C0", alpha=0.6)
    ax.axhline(-depth_ref * 1e6, color="C3", ls="--", lw=0.8)
    ax.set_title(f"TIC {tic} S{sector} TESScut aperture photometry (detrended); target {edge_px:.1f} px from "
                 f"detector edge; lines: top box-search dips, red = passes", fontsize=9)
    ax.set_xlabel("BTJD")
    ax.set_ylabel("ppm")
    fig.tight_layout()
    fig.savefig(out / f"tic{tic}_s{sector}_tesscut.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
