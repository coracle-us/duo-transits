"""
Gaia DR3 binarity flags + TIC mass for EB timing systems, reusing singles/scripts/gaia_flags.py
(query(), rv_variable()). Brightest Gaia source within 5.4 arcsec of the catalogue position.
Adds TIC v8.2 mass/Teff (MAST) for the third-body mass estimate.

Usage:  python sidequests/eb_timing/scripts/gaia_check.py TIC [TIC ...]
Output: sidequests/eb_timing/data/gaia_tic_flags.csv (appends/updates rows)
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
import gaia_flags as gf  # noqa: E402

OUT = HERE / "data" / "gaia_tic_flags.csv"


def tic_info(tic):
    from astroquery.mast import Catalogs
    try:
        r = Catalogs.query_criteria(catalog="TIC", ID=int(tic)).to_pandas()
        if len(r):
            return dict(tic_mass=r["mass"].iloc[0], tic_teff=r["Teff"].iloc[0], tic_rad=r["rad"].iloc[0],
                        tic_dist=r["d"].iloc[0])
    except Exception as e:  # noqa: BLE001
        print(f"TIC query failed {tic}: {e}")
    return {}


def main():
    cat = pd.read_csv(HERE / "data" / "prsa2022_tess_ebs.csv").drop_duplicates("TIC").set_index("TIC")
    old = pd.read_csv(OUT) if OUT.exists() else pd.DataFrame(columns=["tic"])
    rows = {int(r["tic"]): r for r in old.to_dict("records")}
    for tic in [int(x) for x in sys.argv[1:]]:
        ra, dec = cat.loc[tic, "RAJ2000"], cat.loc[tic, "DEJ2000"]
        g = gf.query(f"""SELECT source_id, phot_g_mean_mag, ruwe, non_single_star, rv_chisq_pvalue,
                     rv_renormalised_gof, rv_amplitude_robust, rv_nb_transits,
                     DISTANCE(POINT({ra},{dec}),POINT(ra,dec))*3600 AS sep
                     FROM gaiadr3.gaia_source
                     WHERE 1=CONTAINS(POINT(ra,dec),CIRCLE({ra},{dec},0.0015))
                     ORDER BY phot_g_mean_mag""")
        s = dict(tic=tic)
        if g is not None and len(g):
            s.update(g.iloc[0].to_dict())
            sid = int(g["source_id"].iloc[0])
            s["source_id"] = sid
            s["n_gaia_5arcsec"] = len(g)
            n = gf.query(f"SELECT nss_solution_type, period, eccentricity, semi_amplitude_primary "
                         f"FROM gaiadr3.nss_two_body_orbit WHERE source_id = {sid}")
            if n is not None and len(n):
                s.update({f"nss_{k}": v for k, v in n.iloc[0].to_dict().items()})
            vs = gf.query(f"SELECT source_id FROM gaiadr3.vari_eclipsing_binary WHERE source_id = {sid}")
            s["gaia_vari_eb"] = bool(vs is not None and len(vs))
            s["rv_variable"] = gf.rv_variable(s.get("rv_chisq_pvalue"), s.get("rv_renormalised_gof"),
                                              s.get("rv_nb_transits"))
        s.update(tic_info(tic))
        rows[tic] = s
        print(tic, {k: s.get(k) for k in ("ruwe", "non_single_star", "nss_nss_solution_type", "nss_period",
                                           "tic_mass", "tic_teff")}, flush=True)
    pd.DataFrame(list(rows.values())).to_csv(OUT, index=False)


if __name__ == "__main__":
    main()
