"""
Download the Prsa et al. 2022 TESS EB catalog (VizieR J/ApJS/258/16, 4584 EBs from 2-min
data of Sectors 1-26) to sidequests/eb_timing/data/prsa2022_tess_ebs.csv.
Usage:  python sidequests/eb_timing/scripts/get_catalog.py
"""
from pathlib import Path

from astroquery.vizier import Vizier

OUT = Path(__file__).resolve().parents[1] / "data" / "prsa2022_tess_ebs.csv"


def main():
    v = Vizier(row_limit=-1, columns=["**"])
    t = v.get_catalogs("J/ApJS/258/16")[0]
    df = t.to_pandas()
    df = df.drop(columns=[c for c in ("TESSebs", "Simbad") if c in df.columns])
    df.to_csv(OUT, index=False)
    print(f"{len(df)} rows ({df.TIC.nunique()} unique TIC) -> {OUT}")


if __name__ == "__main__":
    main()
