"""
Run mtphot + period on a list of (number, sector[, known period]) and collect a summary table.

Usage: python run_batch.py LIST.csv OUT.csv [NTHREADS_DOWNLOAD]
LIST.csv columns: number, sector, P_known_h (optional), label (optional), set (optional)
Downloads run in parallel threads (TESScut is slow, ~1-2 min per 27-d 200-s cutout); photometry
and the period search run serially afterwards.
"""
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mtphot  # noqa: E402
import period  # noqa: E402


def dl(row):
    try:
        mtphot.fetch_cutout(f"{int(row.number)};", int(row.sector))
        return None
    except Exception as e:  # noqa: BLE001
        return f"{row.number} S{row.sector}: {e}"


def main():
    lst = pd.read_csv(sys.argv[1])
    out = Path(sys.argv[2])
    nth = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    with ThreadPoolExecutor(nth) as ex:
        errs = [e for e in ex.map(dl, lst.itertuples()) if e]
    for e in errs:
        print("DOWNLOAD FAIL", e, flush=True)
    rows = []
    for r in lst.itertuples():
        obj = f"{int(r.number)};"
        Pk = getattr(r, "P_known_h", np.nan)
        Pk = None if pd.isna(Pk) else float(Pk)
        label = getattr(r, "label", None)
        try:
            _, meta = mtphot.run(obj, int(r.sector))
            res = period.analyse(str(int(r.number)), int(r.sector), Pk,
                                 label=f"({int(r.number)}) {label}" if isinstance(label, str) else None)
            res.update({f"m_{k}": v for k, v in meta.items() if k not in ("obj", "sector")})
        except Exception as e:  # noqa: BLE001
            traceback.print_exc()
            res = dict(id=str(int(r.number)), sector=int(r.sector), status=f"ERR {e}"[:120])
        res["label"] = label
        res["set"] = getattr(r, "set", None)
        rows.append(res)
        print({k: res.get(k) for k in ["id", "sector", "V_med", "n_keep", "P_LS_h", "P2_h", "amp2",
                                        "p_odd_red", "P_known_h", "status"]}, flush=True)
        pd.DataFrame(rows).to_csv(out, index=False)


if __name__ == "__main__":
    main()
