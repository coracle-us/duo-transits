"""Fresh MAST searches (light curves incl. all HLSPs, and TESScut sectors) for the 34 stars whose
catalogued event was not recovered. Cached to singles/revision/data/search/.
Usage: python singles/revision/search_products.py"""
import sys
import time
from pathlib import Path

import pandas as pd

REV = Path(__file__).resolve().parent
ROOT = REV.parents[1]
OUT = REV / "data" / "search"
OUT.mkdir(parents=True, exist_ok=True)


def main(tics=None):
    import lightkurve as lk
    cl = pd.read_csv(ROOT / "singles" / "results" / "final_classification_2026-10-11_rnaas.csv")
    if tics is None:
        tics = cl[cl.cls == "event-not-recovered"].TIC.tolist()
    for tic in tics:
        f1, f2 = OUT / f"{tic}_lc.csv", OUT / f"{tic}_tesscut.csv"
        if f1.exists() and f2.exists():
            continue
        for attempt in range(3):
            try:
                sr = lk.search_lightcurve(f"TIC {tic}")
                tab = sr.table.to_pandas()
                cols = [c for c in ["mission", "author", "exptime", "target_name", "distance", "dataURI",
                                    "productFilename", "provenance_name"] if c in tab.columns]
                tab[cols].to_csv(f1, index=False)
                sc = lk.search_tesscut(f"TIC {tic}")
                sc.table.to_pandas()[["mission"]].to_csv(f2, index=False)
                break
            except Exception as e:  # noqa: BLE001
                print(f"TIC {tic} retry {attempt + 1}: {e}")
                time.sleep(5)
        print(f"TIC {tic}: {len(tab)} LC products", flush=True)


if __name__ == "__main__":
    main([int(x) for x in sys.argv[1:]] or None)
