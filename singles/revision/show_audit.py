"""Print a compact per-star view of the audit outputs (for hand classification)."""
import json
import sys
from pathlib import Path

import pandas as pd

REV = Path(__file__).resolve().parent
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
meta = json.loads((REV / "audit_meta.json").read_text())
m = pd.read_csv(REV / "audit_measurements.csv")
s = pd.read_csv(REV / "audit_sector_scan.csv")
tics = [int(x) for x in sys.argv[1:]] or [int(k) for k in meta]
for tic in tics:
    d = meta[str(tic)]
    print("=" * 100)
    print(f"TIC {tic} CTOI {d['CTOI']} {d['user']} T={d['Tmag']:.1f} epoch {d['epoch_btjd']:.4f} depth {d['depth_ppm']:.0f} "
          f"D {d['dur_h']} h | epoch sector {d['epoch_sector']} observed {d['observed_epoch_sectors']} "
          f"pipeline {d['pipeline_product']}")
    print("   products at epoch sector:", d["products_at_epoch_sector"], "| tesscut:", d["tesscut_sectors"])
    for c in d.get("exofop_params", []):
        print("   exofop:", c)
    x = m[m.TIC == tic].copy()
    if len(x):
        x["pd"] = (x.pipe_depth * 1e6).round(0)
        x["ad"] = (x.all_depth * 1e6).round(0)
        x["sd"] = (x.scan_depth * 1e6).round(0)
        print(x[["sector", "product", "exptime", "column", "pipeline_used", "pipe_n_in", "n_flagged_in", "pipe_verdict",
                 "pd", "pipe_snr", "all_verdict", "ad", "all_snr", "scan_dt_h", "sd", "scan_snr"]]
              .round(2).to_string(index=False))
    y = s[s.TIC == tic]
    if len(y):
        y = y.copy()
        y["top_depth_ppm"] = (y.sec_top_depth * 1e6).round(0)
        cols = [c for c in ["sector", "product", "note", "sec_top_T", "top_depth_ppm", "sec_top_snr", "sec_n_dips", "sec_dips"]
                if c in y]
        print(y[cols].round(3).to_string(index=False))
