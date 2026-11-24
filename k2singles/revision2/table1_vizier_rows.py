"""
AAS82225 revision: full VizieR rows for the Table 1 stars' fixed-list matches (vizier_targeted.csv), for classification
under PLAN.md "Literature classification: rules" (a VizieR row with a period column must be read). One query per star,
all columns, same catalogues and 3" radius as vizier_targeted.py; every non-empty field is written out.
Output: data/literature/table1_vizier_rows.txt and table1_vizier_rows.csv (one row per star and table)
Usage:  python k2singles/revision2/table1_vizier_rows.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import vizier_targeted as vt  # noqa: E402

LIT = HERE / "data" / "literature"


def main():
    t1 = pd.read_csv(HERE.parents[0] / "rnaas" / "table.csv")
    stars = sorted(set(t1.EPIC.astype(str).str.extract(r"(\d{9})")[0].astype(int)))
    v = pd.read_csv(LIT / "vizier_targeted.csv")
    cat = pd.read_csv(HERE.parents[0] / "data" / "k2_singles.csv").drop_duplicates("EPIC").set_index("EPIC")
    lines, rows = [], []
    for e in stars:
        tabs = sorted(v[v.EPIC == e].catalog.unique())
        if not tabs:
            continue
        c = cat.loc[e]
        res = vt.patient({"-source": " ".join(tabs + [vt.GAIA, vt.TMASS]), "-c": f"{c.ra:.6f} {c.dec:+.6f}", "-c.rs": 3,
                          "-out.max": 5, "-out.all": 1, "-out.add": "_r"}, f"EPIC {e}",
                         # every table the main run matched must come back: the mirror also drops single catalogues
                         # from otherwise complete replies (I/357/tbosb1 for 211892898 on the first try)
                         ok=lambda x: set(tabs) <= {t.name for t, _ in vt.tables(x) if len(t.to_table())})
        lines.append(f"\n==================== EPIC {e}")
        got = set()
        for t, desc in vt.tables(res):
            if t.name in (vt.GAIA, vt.TMASS) or t.name not in tabs:
                continue
            tab = t.to_table(use_names_over_ids=True)
            got.add(t.name)
            for r in tab:
                f = {k: str(r[k]) for k in tab.colnames
                     if str(r[k]).strip() not in ("", "--", "nan", "None", "-99.999") and not (np.ma.is_masked(r[k]))}
                txt = "; ".join(f"{k}={x}" for k, x in f.items())
                lines.append(f"--- {t.name}  ({desc[:100]})\n    {txt}")
                rows.append(dict(EPIC=e, catalog=t.name, title=desc[:160], fields=txt))
        missing = sorted(set(tabs) - got)
        if missing:
            lines.append(f"(not returned on this query: {' '.join(missing)})")
    (LIT / "table1_vizier_rows.txt").write_text("\n".join(lines) + "\n")
    pd.DataFrame(rows).to_csv(LIT / "table1_vizier_rows.csv", index=False)
    print(f"{len(rows)} rows on {len({r['EPIC'] for r in rows})} stars")


if __name__ == "__main__":
    main()
