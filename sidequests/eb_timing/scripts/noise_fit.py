"""
Run the O-C noise model + bootstrap (ocfit.analyse) on saved eclipse timings.

Usage:  python sidequests/eb_timing/scripts/noise_fit.py [--oc=DIR] [--out=CSV] [--cal=F] [--nospot] TIC [TIC ...]
  --oc   directory with <TIC>_oc.csv (default results/oc)
  --out  output csv (default results/noise_fit.csv); rows are updated per TIC
  --cal  global sigma calibration factor (default ocfit.CAL)
Reads only rows without a flag (outliers removed by run_etv) and uses sig_s as sig.
"""
import json
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ocfit  # noqa: E402

HERE = Path(__file__).resolve().parents[1]


def load_oc(path):
    oc = pd.read_csv(path)
    oc = oc[oc.flag.isna() | (oc.flag == "")]
    return pd.DataFrame({"T": oc["T"], "kind": oc.kind, "sector": oc.sector, "sig": oc.sig_s, "oc": oc.oc})


def main():
    args = sys.argv[1:]
    ocdir, out, spot = HERE / "results" / "oc", HERE / "results" / "noise_fit.csv", True
    while args and args[0].startswith("--"):
        a = args.pop(0)
        if a == "--nospot":
            spot = False
            continue
        k, v = a[2:].split("=", 1)
        if k == "oc":
            ocdir = HERE / v
        elif k == "out":
            out = HERE / v
        elif k == "cal":
            ocfit.CAL = float(v)
    old = pd.read_csv(out) if out.exists() else pd.DataFrame(columns=["tic"])
    rows = {int(r["tic"]): r for r in old.to_dict("records")}
    for tic in [int(x) for x in args]:
        f = ocdir / f"{tic}_oc.csv"
        if not f.exists():
            print(f"TIC {tic}: no O-C file")
            continue
        t0 = time.time()
        fj = ocdir / f"{tic}_fit.json"
        P1 = json.loads(fj.read_text())["P"] if fj.exists() else None
        r = ocfit.analyse(load_oc(f), spot=spot, P1=P1)
        r = {k: v for k, v in r.items() if k != "bics"}
        r["tic"] = tic
        r["spot_model"] = spot
        rows[tic] = r
        print(f"TIC {tic}: {r['cls']:9s} (v1 {r['cls_v1']}, veto {r['v2_veto'] or '-'}) P3={r['P3']:.1f}±{r['P3_sig']:.1f} A={r['A_s']:.1f}±{r['A_sig']:.1f}s "
              f"e={r['e3']:.2f}±{r['e3_sig']:.2f} dBIC={r['dBIC']:.1f} snr={r['snr']:.1f} rho={r['rho']:.2f} "
              f"tau={r['tau_s']:.1f}s sw={r['white_scale']:.2f} ({time.time() - t0:.0f}s)", flush=True)
        pd.DataFrame(list(rows.values())).to_csv(out, index=False)


if __name__ == "__main__":
    main()
