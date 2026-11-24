"""
Leads triage, step 1: catalogue rows for the 47 lead stars (VizieR row queries only; no files are downloaded).

Per star, around the K2 catalogue position (k2singles/data/k2_singles.csv):
  B/vsx/vsx          within 5"  : VSX variability type and period
  I/358/vclassre     within 5"  : Gaia DR3 variability classifier class (best class)
  I/355/gaiadr3      within 60" : the nearest source within 5" is taken as the target (RUWE, G, NSS flag,
                                  ipd_frac_multi_peak); the rest are neighbours (count within 30" and 60" with
                                  G < G_target + 5, the brightest one)
The existing literature products (results/literature/vizier_targeted.csv, gaia_nss.csv) are read by triage.py;
this adds the VSX type and Gaia RUWE they do not keep.

Output: k2singles/leads_triage/catalog_lookup.csv (re-running skips stars already in it)
Usage:  python k2singles/leads_triage/catalog_lookup.py
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "revision2"))
import vizier_targeted as vt  # noqa: E402  (same CfA mirror and reply checks as the literature search)
from leads import lead_epics  # noqa: E402

ROOT = HERE.parents[1]
OUT = HERE / "catalog_lookup.csv"


def q(cat, ra, dec, rad):
    """Rows of one catalogue within rad arcsec (pandas, with _r), [] -> empty frame, None if the mirror failed."""
    for attempt in range(4):
        v = vt.get({"-source": cat, "-c": f"{ra:.6f} {dec:+.6f}", "-c.rs": rad, "-out.max": 500, "-out.all": 1,
                    "-out.add": "_r"})
        if v is not None:
            tabs = [t.to_table(use_names_over_ids=True) for t, _ in vt.tables(v)]
            tabs = [t for t in tabs if len(t)]
            return tabs[0].to_pandas() if tabs else pd.DataFrame()
        time.sleep(20 * (attempt + 1))
    return None


def main():
    cat = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv", dtype={"EPIC": str}).drop_duplicates("EPIC")
    cat = cat.set_index("EPIC")
    done = pd.read_csv(OUT, dtype={"EPIC": str}) if OUT.exists() else pd.DataFrame(columns=["EPIC"])
    rows = done.to_dict("records")
    for e in lead_epics():
        if e in set(done.EPIC):
            continue
        r = cat.loc[e]
        row = dict(EPIC=e, TIC=int(r.TIC), status="ok")
        vsx = q("B/vsx/vsx", r.ra, r.dec, 5)
        vcl = q("I/358/vclassre", r.ra, r.dec, 5)
        gai = q("I/355/gaiadr3", r.ra, r.dec, 60)
        if vsx is None or vcl is None or gai is None:
            row["status"] = "query-failed"
        if vsx is not None and len(vsx):
            v = vsx.sort_values("_r").iloc[0]
            row.update(vsx_name=v.get("Name"), vsx_type=v.get("Type"), vsx_period=v.get("Period"))
        if vcl is not None and len(vcl):
            v = vcl.sort_values("_r").iloc[0]
            row.update(gaia_vclass=v.get("Class"), gaia_vclass_score=v.get("ClassSc"))
        if gai is not None and len(gai):
            g = gai.sort_values("_r")
            tgt = g[g._r <= 5]
            if len(tgt):
                t = tgt.iloc[0]
                row.update(gaia_source=str(t.Source), gaia_sep=round(float(t._r), 2), G=t.Gmag, RUWE=t.RUWE,
                           NSS=t.NSS, ipd_frac_multi_peak=t.IPDfmp)
                nb = g[(g.Source != t.Source) & (g.Gmag < t.Gmag + 5)]
                row.update(n_nb30=int((nb._r <= 30).sum()), n_nb60=int(len(nb)))
                if len(nb):
                    b = nb.sort_values("Gmag").iloc[0]
                    row.update(brightest_nb_sep=round(float(b._r), 1), brightest_nb_dG=round(float(b.Gmag - t.Gmag), 2))
                close = nb[nb._r <= 30].sort_values("_r")
                if len(close):
                    row.update(nearest_nb_sep=round(float(close.iloc[0]._r), 1),
                               nearest_nb_dG=round(float(close.iloc[0].Gmag - t.Gmag), 2))
        print(e, row["status"], {k: row.get(k) for k in ("vsx_type", "vsx_period", "gaia_vclass", "RUWE", "n_nb30")})
        rows.append(row)
        pd.DataFrame(rows).to_csv(OUT, index=False)


if __name__ == "__main__":
    main()
