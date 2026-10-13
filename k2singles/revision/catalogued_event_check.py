"""
Is each catalogued K2 single event present in the K2 light curve at its catalogued epoch? Pre-empts the AAS82207
editor's point 6 ("do you detect the known transit?") for AAS82225. Uses the same window test as the repeat search
(k2lib.test_window_k2, K2 PDCSAP only) at the LaCourse & Jacobs (2018) epoch, depth and duration, and also the best
window within +-0.1 d (to separate epoch rounding from a missing event).
Output: k2singles/revision/catalogued_event_check.csv and a summary on stdout.
Usage:  python k2singles/revision/catalogued_event_check.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "k2singles" / "scripts"))
import k2lib  # noqa: E402
import run_full3  # noqa: E402


def main():
    rows = []
    for _, r in run_full3.select().iterrows():
        D, delta, T = r.dur_h / 24, r.depth_ppm * 1e-6, float(r.epoch_btjd)
        out = dict(EPIC=int(r.EPIC), event_idx=int(r.event_idx), depth_ppm=r.depth_ppm, dur_h=r.dur_h, epoch=T)
        k2 = k2lib.k2_lightcurves(int(r.EPIC))
        if k2 is None:
            rows.append({**out, "verdict": "no-k2-data"})
            continue
        k2 = k2.sort_values("time")
        t, f, sg = k2.time.to_numpy(), k2.flux.to_numpy(), k2.segment.to_numpy()
        w = k2lib.test_window_k2(t, f, T, D, delta, sg)
        best = max((k2lib.test_window_k2(t, f, T + dt, D, delta, sg) for dt in np.arange(-0.1, 0.1001, D / 8)),
                   key=lambda x: x.get("snr", -np.inf))
        rows.append({**out, "verdict": w["verdict"], "depth_ratio": w.get("depth", np.nan) / delta, "snr": w.get("snr"),
                     "best_dt_h": (best["T_pred"] - T) * 24, "best_verdict": best["verdict"],
                     "best_depth_ratio": best.get("depth", np.nan) / delta, "best_snr": best.get("snr")})
    d = pd.DataFrame(rows)
    d.to_csv(HERE / "catalogued_event_check.csv", index=False)
    print(f"{len(d)} events; at the catalogued epoch: {d.verdict.value_counts().to_dict()}")
    print(f"best within +-0.1 d: {d.best_verdict.value_counts().to_dict()}")
    print(d[d.best_verdict != "supported"].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
