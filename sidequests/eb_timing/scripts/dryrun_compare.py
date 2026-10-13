"""
Dry run check for the frozen southern-run code: compare a fresh run (results/oc_dryrun/,
results/noise_dryrun.csv) with the committed pilot results (results/oc/, results/noise_fit.csv).
Reports, per system: number of eclipses in each, max |T_new - T_old| over matched eclipses (s),
and the differences in P3, A, e, their sigmas, class.
Usage:  python sidequests/eb_timing/scripts/dryrun_compare.py TIC [TIC ...]
Output: results/dryrun_compare.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]


def main():
    old_n = pd.read_csv(HERE / "results" / "noise_fit.csv").set_index("tic")
    new_n = pd.read_csv(HERE / "results" / "noise_dryrun.csv").set_index("tic")
    rows = []
    for tic in [int(x) for x in sys.argv[1:]]:
        o = pd.read_csv(HERE / "results" / "oc" / f"{tic}_oc.csv")
        n = pd.read_csv(HERE / "results" / "oc_dryrun" / f"{tic}_oc.csv")
        m = o.merge(n, on=["kind", "epoch"], suffixes=("_o", "_n"))
        r = dict(tic=tic, n_old=len(o), n_new=len(n), n_matched=len(m),
                 max_dT_s=float(np.max(np.abs(m.T_o - m.T_n)) * 86400) if len(m) else np.nan,
                 max_dsig_s=float(np.max(np.abs(m.sig_s_o - m.sig_s_n)) * 86400) if len(m) else np.nan,
                 identical_times=bool(len(o) == len(n) == len(m) and np.allclose(m.T_o, m.T_n, atol=1e-9, rtol=0)))
        for k in ("P3", "P3_sig", "A_s", "A_sig", "e3", "dBIC"):
            r[f"d_{k}"] = float(new_n.loc[tic, k] - old_n.loc[tic, k])
        r["cls_old"], r["cls_new"] = old_n.loc[tic, "cls"], new_n.loc[tic, "cls"]
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "results" / "dryrun_compare.csv", index=False)
    print(df.to_string())


if __name__ == "__main__":
    main()
