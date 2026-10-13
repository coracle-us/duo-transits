"""
Gaia DR3 binarity flags for the pilot hosts (same queries and RV-variability cut as
singles/scripts/gaia_flags.py): RUWE, non_single_star, RV variability, NSS two-body orbit.
Output: toi_singles_pilot/gaia_flags.csv
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
import gaia_flags as gf  # noqa: E402
import s01_run  # noqa: E402


def main():
    rows = []
    for r in s01_run.targets("pilot"):
        g = gf.query(f"""SELECT source_id, phot_g_mean_mag, ruwe, non_single_star, rv_chisq_pvalue,
                     rv_renormalised_gof, rv_amplitude_robust, rv_nb_transits,
                     DISTANCE(POINT({r['ra']},{r['dec']}),POINT(ra,dec))*3600 AS sep
                     FROM gaiadr3.gaia_source
                     WHERE 1=CONTAINS(POINT(ra,dec),CIRCLE({r['ra']},{r['dec']},0.0015))
                     ORDER BY phot_g_mean_mag""")
        s = dict(TOI=r["TOI"], TIC=int(r["TIC"]))
        if g is None or not len(g):
            s["note"] = "no Gaia result"
            rows.append(s)
            continue
        s.update({k: v for k, v in g.iloc[0].to_dict().items() if k != "source_id"})
        sid = int(g["source_id"].iloc[0])
        s["source_id"] = str(sid)
        s["n_gaia_within_5.4as"] = len(g)
        n = gf.query(f"SELECT nss_solution_type, period, eccentricity, semi_amplitude_primary "
                     f"FROM gaiadr3.nss_two_body_orbit WHERE source_id = {sid}")
        if n is not None and len(n):
            s.update({f"nss_{k}": v for k, v in n.iloc[0].to_dict().items()})
        s["rv_variable"] = gf.rv_variable(s.get("rv_chisq_pvalue"), s.get("rv_renormalised_gof"),
                                          s.get("rv_nb_transits"))
        s["binarity_flag"] = bool((s.get("non_single_star") or 0) > 0 or s["rv_variable"]
                                  or (s.get("ruwe") or 0) > 1.4)
        rows.append(s)
        print(s, flush=True)
    pd.DataFrame(rows).to_csv(HERE.parent / "gaia_flags.csv", index=False)


if __name__ == "__main__":
    main()
