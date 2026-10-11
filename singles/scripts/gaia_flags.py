"""
Gaia DR3 binarity flags for every target in singles/data/ctoi_singles.csv.

Why: in batch 2, two 'single transit' CTOIs (231969683, 230130556) turned out to be eccentric
eclipsing binaries whose TESS periods match Gaia DR3 SB1 orbits, and TIC 147035057 is
strongly RV-variable. Checking Gaia first is cheap and flags such systems before any
light-curve work.

Per target (brightest Gaia source within 5.4 arcsec): RUWE, non_single_star, RV-variability
statistics, and any nss_two_body_orbit solution (type, period, e, K1).
rv_variable := rv_chisq_pvalue <= 0.01 & rv_renormalised_gof > 4 & rv_nb_transits >= 10
(the cut used for Gaia DR3 RV variability, Katz et al. 2023).

One small synchronous query per star (the archive drops bulk upload jobs); resumable.
Output: singles/data/gaia_binarity_flags.csv
Usage:  python singles/scripts/gaia_flags.py
"""
import time
from pathlib import Path

import pandas as pd
from astroquery.gaia import Gaia

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "singles" / "data" / "gaia_binarity_flags.csv"


def query(sql):
    for attempt in range(4):
        try:
            return Gaia.launch_job(sql).get_results().to_pandas()
        except Exception as e:  # noqa: BLE001
            print(f"  retry {attempt + 1}: {str(e)[:80]}")
            time.sleep(10 * (attempt + 1))
    return None


def rv_variable(p, gof, n):
    """Katz et al. 2023 cut. NaN/None means 'not tested' (False). Note a p-value of exactly 0
    is the most variable case; an earlier `(p or 1)` idiom wrongly turned it into 1."""
    if any(v is None or pd.isna(v) for v in (p, gof, n)):
        return False
    return bool(p <= 0.01 and gof > 4 and n >= 10)


def main():
    c = pd.read_csv(ROOT / "singles" / "data" / "ctoi_singles.csv").drop_duplicates("TIC")
    done = pd.read_csv(OUT) if OUT.exists() else pd.DataFrame(columns=["tic"])
    rows = done.to_dict("records")
    have = set(done.tic.astype(int)) if len(done) else set()
    for i, r in enumerate(c.itertuples(), 1):
        if int(r.TIC) in have:
            continue
        g = query(f"""SELECT source_id, phot_g_mean_mag, ruwe, non_single_star, rv_chisq_pvalue,
                     rv_renormalised_gof, rv_amplitude_robust, rv_nb_transits,
                     DISTANCE(POINT({r.RA},{r.Dec}),POINT(ra,dec))*3600 AS sep
                     FROM gaiadr3.gaia_source
                     WHERE 1=CONTAINS(POINT(ra,dec),CIRCLE({r.RA},{r.Dec},0.0015))
                     ORDER BY phot_g_mean_mag""")
        if g is None:
            rows.append(dict(tic=int(r.TIC), note="query failed"))
            continue
        if not len(g):
            rows.append(dict(tic=int(r.TIC), note="no Gaia source within 5.4 arcsec"))
            continue
        s = g.iloc[0].to_dict()
        # source_id must stay an exact int64: a row Series of all-numeric columns is float64,
        # which corrupts 19-digit IDs (bug fixed 2026-10-11)
        sid = int(g["source_id"].iloc[0])
        s["source_id"] = sid
        s["tic"] = int(r.TIC)
        n = query(f"SELECT nss_solution_type, period, eccentricity, semi_amplitude_primary "
                  f"FROM gaiadr3.nss_two_body_orbit WHERE source_id = {sid}")
        if n is not None and len(n):
            s.update({f"nss_{k}": v for k, v in n.iloc[0].to_dict().items()})
        s["rv_variable"] = rv_variable(s.get("rv_chisq_pvalue"), s.get("rv_renormalised_gof"),
                                       s.get("rv_nb_transits"))
        rows.append(s)
        if i % 20 == 0:
            pd.DataFrame(rows).to_csv(OUT, index=False)
            print(f"{i}/{len(c)}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    flag = (df.non_single_star.fillna(0) > 0) | df.rv_variable.fillna(False) | (df.ruwe.fillna(0) > 1.4)
    print(f"done: {len(df)} stars; any binarity flag: {int(flag.sum())}; NSS orbits: "
          f"{df.get('nss_nss_solution_type', pd.Series(dtype=str)).notna().sum()}; "
          f"RV-variable: {int(df.rv_variable.fillna(False).sum())}")


if __name__ == "__main__":
    main()
