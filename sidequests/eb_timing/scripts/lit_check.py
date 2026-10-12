"""
Literature/identity check for EB timing systems: SIMBAD (main id, object type, number of
references) and AAVSO VSX (VizieR B/vsx: name, type, period) within 10 arcsec of the
Prsa et al. 2022 position, plus membership in data/known_multiples.csv.

SIMBAD bibcodes of the matched object are listed so the reference list can be checked
for prior ETV / third-body work.

Usage:  python sidequests/eb_timing/scripts/lit_check.py TIC [TIC ...]
Output: sidequests/eb_timing/data/lit_check.csv (rows updated per TIC)
"""
import sys
from pathlib import Path

import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "data" / "lit_check.csv"


def main():
    cat = pd.read_csv(HERE / "data" / "prsa2022_tess_ebs.csv").drop_duplicates("TIC").set_index("TIC")
    known = pd.read_csv(HERE / "data" / "known_multiples.csv").set_index("TIC").papers
    old = pd.read_csv(OUT) if OUT.exists() else pd.DataFrame(columns=["tic"])
    rows = {int(r["tic"]): r for r in old.to_dict("records")}
    sim = Simbad()
    sim.add_votable_fields("otype", "nbref")
    viz = Vizier(columns=["Name", "Type", "Period", "max", "min"], row_limit=5)
    for tic in [int(x) for x in sys.argv[1:]]:
        c = SkyCoord(cat.loc[tic, "RAJ2000"] * u.deg, cat.loc[tic, "DEJ2000"] * u.deg)
        row = dict(tic=tic, known_multiple=known.get(tic, ""))
        try:
            s = sim.query_region(c, radius=10 * u.arcsec)
            if s is not None and len(s):
                s = s.to_pandas()
                s.columns = [x.lower() for x in s.columns]
                row.update(simbad_id=s.iloc[0].get("main_id"), simbad_otype=s.iloc[0].get("otype"),
                           simbad_nbref=s.iloc[0].get("nbref"), simbad_n=len(s))
                mid = str(s.iloc[0].get("main_id")).replace("'", "''")
                q = ("SELECT b.bibcode FROM basic AS o JOIN has_ref AS h ON h.oidref=o.oid "
                     f"JOIN ref AS b ON b.oidbib=h.oidbibref WHERE o.main_id='{mid}'")
                refs = Simbad.query_tap(q).to_pandas()
                row["simbad_bibcodes"] = ";".join(sorted(refs.bibcode.astype(str)))
        except Exception as e:  # noqa: BLE001
            row["simbad_err"] = str(e)[:80]
        try:
            v = viz.query_region(c, radius=10 * u.arcsec, catalog="B/vsx/vsx")
            if len(v):
                v = v[0].to_pandas()
                row.update(vsx_name=v.iloc[0]["Name"], vsx_type=v.iloc[0]["Type"], vsx_period=v.iloc[0]["Period"])
        except Exception as e:  # noqa: BLE001
            row["vsx_err"] = str(e)[:80]
        rows[tic] = row
        print(row, flush=True)
    pd.DataFrame(list(rows.values())).to_csv(OUT, index=False)


if __name__ == "__main__":
    main()
