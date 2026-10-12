"""
VizieR census for K2-singles leads: every VizieR catalogue with an entry within 3" of the star,
plus any column that looks like a period. Used to check whether an eclipse period is already
catalogued before we claim it (TOI-5560.01 / K2-158 lesson).
Output: k2singles/data/vizier_census.csv
Usage:  python k2singles/scripts/vizier_census.py EPIC [EPIC ...]
"""
import sys
from pathlib import Path

import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "k2singles" / "data" / "vizier_census.csv"
SKIP = ("I/", "II/246", "II/328", "II/336", "II/349", "II/363", "IV/38", "IV/39", "V/154", "J/AJ/156/58")


def main():
    cat = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv").drop_duplicates("EPIC").set_index("EPIC")
    rows = pd.read_csv(OUT).to_dict("records") if OUT.exists() else []
    v = Vizier(columns=["**"], row_limit=5)
    for e in map(int, sys.argv[1:]):
        r = cat.loc[e]
        res = v.query_region(SkyCoord(r.ra, r.dec, unit="deg"), radius=3 * u.arcsec)
        for name in res.keys():
            t = res[name]
            if name.startswith(SKIP):
                continue
            pcols = [c for c in t.colnames if c.lower() in ("per", "period", "p", "porb", "plc", "pdays", "p_orb")
                     or c.lower().startswith("per")]
            per = "; ".join(f"{c}={t[c][0]}" for c in pcols)
            rows.append(dict(EPIC=e, TIC=int(r.TIC), catalog=name, title=t.meta.get("description", "")[:90],
                             periods=per))
        print(e, len(res), "catalogues", flush=True)
    pd.DataFrame(rows).drop_duplicates(["EPIC", "catalog"]).to_csv(OUT, index=False)


if __name__ == "__main__":
    main()
