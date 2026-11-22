"""
AAS82225 revision, literature search part 3: Gaia DR3 orbits and eclipsing-binary periods for all 163 stars.
The VizieR census (vizier_all.py) skips the "I/" section and keeps only I/358 (variability), so it misses I/357, the
Gaia DR3 non-single-star catalogue, which PLAN.md says must be read (it carries orbital periods). This queries the
Gaia archive directly by source_id (from SIMBAD's identifiers, identifiers.csv):
  gaiadr3.nss_two_body_orbit    SB1 / SB2 / EclipsingBinary / AstroSpectroSB1 / Orbital solutions, with period
  gaiadr3.vari_eclipsing_binary the variability pipeline's eclipsing-binary model, with frequency (1/period)
  gaiadr3.nss_non_linear_spectro, nss_acceleration_astro, nss_vim_fl  solution types without a period (listed only)
Server: Gaia archive TAP (gea.esac.esa.int). Output: k2singles/revision2/data/literature/gaia_nss.csv
Stars with no SIMBAD Gaia id: the brightest Gaia DR3 source within 2" (gaia_ids.csv says which).
Usage:  python k2singles/revision2/gaia_nss.py
"""
import io
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
LIT = HERE / "data" / "literature"
TAP = "https://gea.esac.esa.int/tap-server/tap/sync"


def tap(q):
    data = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q)).encode()
    for k in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(TAP, data=data), timeout=300) as r:
                return pd.read_csv(io.BytesIO(r.read()))
        except Exception as e:  # noqa: BLE001
            print("retry", k, str(e)[:120], flush=True)
            time.sleep(10 * (k + 1))
    raise RuntimeError("Gaia TAP failed")


def main():
    ids = pd.read_csv(LIT / "identifiers.csv")
    rows = []
    for r in ids.itertuples():
        m = re.search(r"Gaia DR3 (\d+)", r.identifiers)
        rows.append(dict(EPIC=r.EPIC, TIC=r.TIC, source_id=int(m.group(1)) if m else None))
    s = pd.DataFrame(rows)
    s["id_from"] = s.source_id.notna().map({True: "simbad", False: ""})
    print(f"{s.source_id.notna().sum()} of {len(s)} stars have a Gaia DR3 source_id in SIMBAD")
    # the rest (no SIMBAD entry): the brightest Gaia DR3 source within 2" of the catalogue position
    cat = pd.read_csv(HERE.parents[1] / "k2singles" / "data" / "k2_singles.csv").drop_duplicates("EPIC").set_index("EPIC")
    for i in s.index[s.source_id.isna()]:
        c = cat.loc[s.at[i, "EPIC"]]
        x = tap(f"SELECT TOP 1 source_id, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), "
                f"CIRCLE('ICRS', {c.ra}, {c.dec}, {2 / 3600})) ORDER BY phot_g_mean_mag ASC")
        if len(x):
            s.at[i, "source_id"], s.at[i, "id_from"] = int(x.source_id[0]), "position 2 arcsec"
    print(f"after the position match: {s.source_id.notna().sum()} of {len(s)}")
    sid = ",".join(str(int(x)) for x in s.source_id.dropna())
    orb = tap(f"SELECT source_id, nss_solution_type, period, period_error, t_periastron, eccentricity, "
              f"significance, goodness_of_fit FROM gaiadr3.nss_two_body_orbit WHERE source_id IN ({sid})")
    eb = tap(f"SELECT source_id, frequency, frequency_error, global_ranking, model_type "
             f"FROM gaiadr3.vari_eclipsing_binary WHERE source_id IN ({sid})")
    other = []
    for t in ("nss_non_linear_spectro", "nss_acceleration_astro", "nss_vim_fl"):
        x = tap(f"SELECT source_id, nss_solution_type FROM gaiadr3.{t} WHERE source_id IN ({sid})")
        other.append(x)
    other = pd.concat(other) if other else pd.DataFrame(columns=["source_id", "nss_solution_type"])
    out = []
    for r in s.itertuples():
        if pd.isna(r.source_id):
            out.append(dict(EPIC=r.EPIC, TIC=r.TIC, source_id=None, table="(no Gaia DR3 source)"))
            continue
        for x in orb[orb.source_id == r.source_id].itertuples():
            out.append(dict(EPIC=r.EPIC, TIC=r.TIC, source_id=r.source_id, table="nss_two_body_orbit",
                            solution=x.nss_solution_type, period_d=x.period, period_err_d=x.period_error,
                            eccentricity=x.eccentricity, significance=x.significance))
        for x in eb[eb.source_id == r.source_id].itertuples():
            out.append(dict(EPIC=r.EPIC, TIC=r.TIC, source_id=r.source_id, table="vari_eclipsing_binary",
                            solution=x.model_type, period_d=1 / x.frequency,
                            period_err_d=x.frequency_error / x.frequency ** 2, significance=x.global_ranking))
        for x in other[other.source_id == r.source_id].itertuples():
            out.append(dict(EPIC=r.EPIC, TIC=r.TIC, source_id=r.source_id, table="nss (no period)", solution=x.nss_solution_type))
    o = pd.DataFrame(out)
    o.to_csv(LIT / "gaia_nss.csv", index=False)
    s.to_csv(LIT / "gaia_ids.csv", index=False)
    print(o[o.table != "(no Gaia DR3 source)"].to_string(index=False))


if __name__ == "__main__":
    main()
