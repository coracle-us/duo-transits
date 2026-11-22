"""
AAS82225 revision, literature search part 2: VizieR census of all 163 stars (as k2singles/scripts/vizier_census.py,
which covered only the leads). Every VizieR catalogue with an entry within 3" of the star, with any period-like
column. Catalogue papers often name a star only in a machine-readable table, which neither ADS full text nor
SIMBAD's bibliography finds (e.g. IJspeert et al. 2024, Kovalev et al. 2024 for our C0 binaries).
Server: the CfA mirror (vizier.cfa.harvard.edu); the CDS host returned 404 on 2026-11-15.
Output: k2singles/revision2/data/literature/vizier.csv (one row per star and catalogue). Resumable.
A response without Gaia DR3 (I/355) counts as a failure and is not cached.
Usage:  python k2singles/revision2/vizier_all.py [n_threads]
"""
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "data" / "literature"
CACHE = OUT / "vizier_cache"
SKIP = ("I/", "II/246", "II/328", "II/336", "II/349", "II/363", "IV/38", "IV/39", "V/154", "J/AJ/156/58")
KEEP = ("I/358",)  # Gaia DR3 variability (eclipsing binaries with periods), inside the skipped "I/" astrometric section
PCOLS = ("per", "period", "p", "porb", "plc", "pdays", "p_orb", "p1", "pbls", "ptess")


def one(epic, tic, ra, dec):
    f = CACHE / f"{epic}.json"
    if f.exists():
        return json.loads(f.read_text())
    from astroquery.vizier import Vizier
    v = Vizier(columns=["**"], row_limit=5, timeout=300)
    v.VIZIER_SERVER = "vizier.cfa.harvard.edu"
    for k in range(4):
        try:
            res = v.query_region(SkyCoord(ra, dec, unit="deg"), radius=3 * u.arcsec)
        except Exception:  # noqa: BLE001
            time.sleep(20 * (k + 1))
            continue
        # every one of these stars is in Gaia DR3 (I/355): a result without it is a failed or truncated response,
        # not an empty sky (on 2026-11-22 the CfA mirror returned 66 such empty results, which were cached as real)
        if any(name.startswith("I/355") for name in res.keys()):
            break
        time.sleep(20 * (k + 1))
    else:
        raise RuntimeError(f"VizieR failed for EPIC {epic}")
    rows = []
    for name in res.keys():
        if name.startswith(SKIP) and not name.startswith(KEEP):
            continue
        t = res[name]
        pc = [c for c in t.colnames if c.lower() in PCOLS or c.lower().startswith("per")]
        rows.append(dict(EPIC=epic, TIC=tic, catalog=name, title=str(t.meta.get("description", ""))[:120],
                         periods="; ".join(f"{c}={t[c][0]}" for c in pc)))
    f.write_text(json.dumps(rows))
    return rows


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    CACHE.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv").drop_duplicates("EPIC")
    allrows = []
    with ThreadPoolExecutor(n) as ex:
        futs = {ex.submit(one, int(r.EPIC), int(r.TIC), float(r.ra), float(r.dec)): int(r.EPIC) for r in d.itertuples()}
        for k, fu in enumerate(as_completed(futs), 1):
            try:
                rows = fu.result()
            except Exception as e:  # noqa: BLE001
                print(f"[{k}/{len(futs)}] EPIC {futs[fu]} FAILED {e}", flush=True)
                continue
            allrows += rows
            print(f"[{k}/{len(futs)}] EPIC {futs[fu]}: {len(rows)} catalogues, "
                  f"{sum(bool(r['periods']) for r in rows)} with a period column", flush=True)
    pd.DataFrame(allrows).to_csv(OUT / "vizier.csv", index=False)


if __name__ == "__main__":
    main()
