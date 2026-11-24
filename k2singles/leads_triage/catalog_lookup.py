"""
Leads triage, step 1: catalogue rows for the 47 lead stars (row queries only; no files are downloaded).

  Gaia DR3   ESA archive TAP, one cone per star: every gaia_source row within 60" of each star (K2 catalogue position,
             k2singles/data/k2_singles.csv). The nearest within 5" is the target (RUWE, G, non_single_star,
             ipd_frac_multi_peak); the rest are neighbours (counts within 30" and 60" with G < G_target + 5,
             the nearest and the brightest).
  VSX type   VizieR B/vsx/vsx within 5", only for stars whose literature census (results/literature/
             vizier_targeted.csv) already has a VSX row; the census keeps the period but not the type.
  Gaia class VizieR I/358/vclassre within 5", only for stars with a vclassre row in that census.
The VizieR CfA mirror took ~10 min per star for a full 60" Gaia cone on 2026-11-23, hence the split.

Output: k2singles/leads_triage/catalog_lookup.csv
Usage:  python k2singles/leads_triage/catalog_lookup.py
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astroquery.gaia import Gaia

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "revision2"))
import vizier_targeted as vt  # noqa: E402  (same CfA mirror and reply checks as the literature search)
from leads import REV, lead_epics  # noqa: E402

ROOT = HERE.parents[1]
OUT = HERE / "catalog_lookup.csv"
R_NB = 60.0


def q(cat, ra, dec, rad):
    """Rows of one VizieR catalogue within rad arcsec (pandas, with _r); None if the mirror failed."""
    for attempt in range(4):
        v = vt.get({"-source": cat, "-c": f"{ra:.6f} {dec:+.6f}", "-c.rs": rad, "-out.max": 50, "-out.all": 1,
                    "-out.add": "_r"})
        if v is not None:
            tabs = [t.to_table(use_names_over_ids=True) for t, _ in vt.tables(v)]
            tabs = [t for t in tabs if len(t)]
            return tabs[0].to_pandas() if tabs else pd.DataFrame()
        time.sleep(20 * (attempt + 1))
    return None


def gaia(cat):
    """One cone per star (an OR of cones defeats the archive's spatial index)."""
    parts = []
    for r in cat.itertuples():
        for attempt in range(3):
            try:
                j = Gaia.launch_job("SELECT source_id, ra, dec, phot_g_mean_mag, ruwe, non_single_star, ipd_frac_multi_peak "
                                    "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec),"
                                    f"CIRCLE({r.ra:.6f},{r.dec:.6f},{R_NB / 3600:.6f}))")
                parts.append(j.get_results().to_pandas())
                break
            except Exception as e:  # noqa: BLE001
                print("  gaia retry", r.EPIC, e, flush=True)
                time.sleep(10)
    return pd.concat(parts, ignore_index=True).drop_duplicates("source_id")


def main():
    epics = lead_epics()
    cat = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv", dtype={"EPIC": str}).drop_duplicates("EPIC")
    cat = cat.set_index("EPIC").loc[epics]
    census = pd.read_csv(REV / "results" / "literature" / "vizier_targeted.csv", dtype={"EPIC": str})
    ids = pd.read_csv(REV / "results" / "literature" / "gaia_ids.csv", dtype={"EPIC": str}).drop_duplicates("EPIC")
    ids = ids.set_index("EPIC").source_id
    G = gaia(cat.reset_index())
    rows = []
    for e, r in cat.iterrows():
        row = dict(EPIC=e, TIC=int(r.TIC), status="ok")
        dra = (G.ra - r.ra) * np.cos(np.radians(r.dec)) * 3600
        g = G.assign(_r=np.hypot(dra, (G.dec - r.dec) * 3600))
        g = g[g._r <= R_NB].sort_values("_r")
        tgt = g[g._r <= 5]
        if not len(tgt) and e in ids.index:  # high proper motion: match the census's Gaia DR3 id (stored as a float)
            tgt = g[np.abs(g.source_id.astype(float) - float(ids.loc[e])) < 4096]
        if len(tgt):
            t = tgt.iloc[0]
            row.update(gaia_source=str(int(t.source_id)), gaia_sep=round(float(t._r), 2), G=round(float(t.phot_g_mean_mag), 3),
                       RUWE=t.ruwe, NSS=int(t.non_single_star), ipd_frac_multi_peak=t.ipd_frac_multi_peak)
            nb = g[(g.source_id != t.source_id) & (g.phot_g_mean_mag < t.phot_g_mean_mag + 5)]
            row.update(n_nb30=int((nb._r <= 30).sum()), n_nb60=int(len(nb)))
            if len(nb):
                b = nb.sort_values("phot_g_mean_mag").iloc[0]
                n = nb.iloc[0]
                row.update(nearest_nb_sep=round(float(n._r), 1), nearest_nb_dG=round(float(n.phot_g_mean_mag - t.phot_g_mean_mag), 2),
                           brightest_nb_sep=round(float(b._r), 1), brightest_nb_dG=round(float(b.phot_g_mean_mag - t.phot_g_mean_mag), 2))
        else:
            row["status"] = "no Gaia source within 5 arcsec"
        c = census[census.EPIC == e].catalog
        if c.str.startswith("B/vsx").any():
            v = q("B/vsx/vsx", r.ra, r.dec, 5)
            if v is None:
                row["status"] = "vsx query failed"
            elif len(v):
                v = v.sort_values("_r").iloc[0]
                row.update(vsx_name=v.get("Name"), vsx_type=v.get("Type"), vsx_period=v.get("Period"))
        if c.str.startswith("I/358/vclassre").any():
            v = q("I/358/vclassre", r.ra, r.dec, 5)
            if v is None:
                row["status"] = "vclassre query failed"
            elif len(v):
                v = v.sort_values("_r").iloc[0]
                row.update(gaia_vclass=v.get("Class"), gaia_vclass_score=v.get("ClassSc"))
        print(e, row["status"], {k: row.get(k) for k in ("vsx_type", "vsx_period", "gaia_vclass", "RUWE", "n_nb30")}, flush=True)
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT, index=False)


if __name__ == "__main__":
    main()
