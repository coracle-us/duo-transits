"""
Epoch-offset check: re-run the pipeline's own event test (alias_pipeline.test_window on the
pipeline-loaded light curve, same products and quality cuts as pilot.py) at the catalogued epoch
shifted by -0.5 d, 0 and +0.5 d, and find the best dip within +-0.15 d of each shifted epoch.

Motivation: several UAN-working-group ("ExoMAST") CTOIs show their dip ~12 h after the
catalogued epoch, as if a JD/MJD half-day were mishandled. Run on ALL analysed stars by the same
uploader so the recovered ones serve as a control.

Usage:  python singles/revision/offset_check.py [user]   (default: all 34 + all amarino CTOIs)
Output: singles/revision/offset_check.csv
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
REV = Path(__file__).resolve().parent
sys.path.insert(0, str(REV))
import revlib as rl  # noqa: E402

ap = rl.ap
ROOT = rl.ROOT
PRODUCTS = ap.PRODUCT_RANK + ["OWN"]


def best_near(t, f, sg, T0, D, delta, half=0.15):
    best = None
    for T in np.arange(T0 - half, T0 + half + 1e-9, D / 20):
        r = ap.test_window(t, f, T, 0.0, D, delta, sg)
        if "snr" in r and r.get("cover", 0) >= 0.5 and (best is None or r["snr"] > best["snr"]):
            best = dict(r, T=T)
    return best


def main():
    cl = pd.read_csv(ROOT / "singles" / "results" / "final_classification_2026-10-11_rnaas.csv")
    cs = pd.read_csv(ROOT / "singles" / "data" / "ctoi_singles.csv")
    m = cl.merge(cs[["TIC", "User", "epoch_btjd"]].drop_duplicates("TIC"), on="TIC")
    sel = m[(m.cls == "event-not-recovered") | (m.User == "amarino")]
    rows = []
    for _, r in sel.iterrows():
        tic = int(r.TIC)
        D, delta, ep = r.dur_h / 24, r.depth_ppm * 1e-6, float(r.epoch_btjd)
        lc = ap.load_lightcurves(tic, products=PRODUCTS)
        if not len(lc):
            continue
        t, f, sg = lc.time.to_numpy(), lc.flux.to_numpy(), lc.segment.to_numpy()
        out = dict(TIC=tic, user=r.User, cls=r.cls, epoch=ep, depth_ppm=r.depth_ppm, dur_h=r.dur_h)
        for lab, off in (("m05", -0.5), ("0", 0.0), ("p05", 0.5)):
            e = ap.test_window(t, f, ep + off, rl.EPOCH_ERR, D, delta, sg)
            out[f"v_{lab}"] = e["verdict"]
            out[f"d_{lab}"] = round(e.get("depth", np.nan) * 1e6, 0) if "depth" in e else np.nan
            out[f"snr_{lab}"] = round(e.get("snr", np.nan), 1) if "snr" in e else np.nan
            b = best_near(t, f, sg, ep + off, D, delta)
            if b:
                out[f"bestdt_h_{lab}"] = round((b["T"] - ep) * 24, 2)
                out[f"bestd_{lab}"] = round(b["depth"] * 1e6, 0)
                out[f"bestsnr_{lab}"] = round(b["snr"], 1)
                j = min(np.searchsorted(t, b["T"]), len(t) - 1)
                out[f"prod_{lab}"] = f"S{lc.sector.iloc[j]} {lc['product'].iloc[j]}"
        rows.append(out)
        print(out, flush=True)
    pd.DataFrame(rows).to_csv(REV / "offset_check.csv", index=False)


if __name__ == "__main__":
    main()
