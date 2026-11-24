"""
Fetch the raw inputs for the host-rotation retest into data/ (git-ignored), with the download date in each file name.

  VizieR (IDIA mirror, asu-tsv; files <stem>_<table>_<date>.csv):
    J/ApJ/775/L11   McQuillan, Mazeh & Aigrain 2013, KOI rotation periods   -> mma2013_<date>.csv
    J/ApJS/211/24   McQuillan, Mazeh & Aigrain 2014, table1 (periods) + table2 (no period)
    J/ApJS/244/21   Santos et al. 2019 (M and K dwarfs)                     -> santos2019_<date>.csv
    J/ApJS/255/17   Santos et al. 2021 (F and G stars)                      -> santos2021_<date>.csv
  NASA Exoplanet Archive TAP:
    cumulative (KOI table: kepid, koi disposition)                          -> koi_cumulative_<date>.csv
    keplerstellar (KIC -> 2MASS designation)                                -> kic_2mass_<date>.csv
  Gaia DR3, Kepler field, G < 17.5 (ESA Gaia archive TAP, RA strips)       -> gaia_field_<date>.csv.gz
Usage:  .venv/bin/python fetch.py [step ...]      steps: vizier koi kicstellar gaia (default: all; a step is skipped if its file exists)
"""
import datetime as dt
import io
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TODAY = dt.date.today().isoformat()
CATS = {"mma2013": "J/ApJ/775/L11", "mma2014": "J/ApJS/211/24", "santos2019": "J/ApJS/244/21",
        "santos2021": "J/ApJS/255/17"}


def latest(stem):
    """Newest file for a stem (tables are stem_table_date.csv)."""
    fs = sorted(DATA.glob(f"{stem}_*.csv*"))
    return fs[-1] if fs else None


# On 2026-11-24 the CDS host answered 404 and the CfA mirror's database was down ("Postgres connect error"); the
# IDIA mirror (South Africa) answered. Plain asu-tsv, every column, all rows, positions added.
VIZIER_HOST = "vizier.idia.ac.za"
TABLES = {"mma2013": ["table1"], "mma2014": ["table1", "table2"], "santos2019": ["table3", "table4"],
          "santos2021": ["table1", "table2"]}


def parse_tsv(text):
    """VizieR asu-tsv -> DataFrame (header line, units line, dashes line, then data)."""
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    if "QUERY_STATUS=ERROR" in text or len(lines) < 3:
        raise RuntimeError("VizieR error or empty: " + text[:500])
    head = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[3:]]
    d = pd.DataFrame(rows, columns=head).apply(lambda c: c.str.strip())
    for c in d.columns:
        conv = pd.to_numeric(d[c], errors="coerce")
        if conv.notna().sum() >= 0.9 * (d[c] != "").sum():
            d[c] = conv
    return d


def vizier():
    for stem, cat in CATS.items():
        for tab in TABLES[stem]:
            out = DATA / f"{stem}_{tab}_{TODAY}.csv"
            if list(DATA.glob(f"{stem}_{tab}_*.csv")):
                print(stem, tab, "exists")
                continue
            url = (f"https://{VIZIER_HOST}/viz-bin/asu-tsv?" +
                   urllib.parse.urlencode({"-source": f"{cat}/{tab}", "-out.max": "unlimited", "-out.all": "",
                                           "-out.add": "_RAJ,_DEJ"}))
            with urllib.request.urlopen(url, timeout=900) as r:
                d = parse_tsv(r.read().decode("utf-8", "replace"))
            d.to_csv(out, index=False)
            print("wrote", out.name, len(d), "rows; columns:", list(d.columns))


def koi():
    if latest("koi_cumulative"):
        print("koi exists:", latest("koi_cumulative").name)
        return
    q = ("select kepid,kepoi_name,kepler_name,koi_disposition,koi_pdisposition,koi_period,ra,dec "
         "from cumulative")
    url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode(dict(query=q, format="csv"))
    with urllib.request.urlopen(url, timeout=300) as r:
        d = pd.read_csv(io.BytesIO(r.read()))
    out = DATA / f"koi_cumulative_{TODAY}.csv"
    d.to_csv(out, index=False)
    print("wrote", out.name, len(d), "rows")


def kicstellar():
    """KIC -> 2MASS designation, from the NASA Exoplanet Archive keplerstellar table."""
    if latest("kic_2mass"):
        print("kic_2mass exists:", latest("kic_2mass").name)
        return
    q = "select distinct kepid, tm_designation from keplerstellar"
    url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode(dict(query=q, format="csv"))
    with urllib.request.urlopen(url, timeout=900) as r:
        d = pd.read_csv(io.BytesIO(r.read()))
    out = DATA / f"kic_2mass_{TODAY}.csv"
    d.to_csv(out, index=False)
    print("wrote", out.name, len(d), "rows,", d.kepid.nunique(), "KIC")


def gaia():
    """Gaia DR3 for the whole Kepler field (G < 17.5), downloaded in RA strips; the KIC match is done locally
    (analysis.gaia_match: nearest source within 1 arcsec after moving Gaia back to J2000 with its proper motion).
    Why not a server-side match (2026-11-24): gaia-kepler.fun failed its TLS check, CDS XMatch was down, and the
    ESA archive timed out on joins with uploaded tables (positions and 2MASS IDs alike). Anonymous queries return at
    most 50,000 rows, so a strip that hits the cap is split in two."""
    if latest("gaia_field"):
        print("gaia exists:", latest("gaia_field").name)
        return
    from astroquery.gaia import Gaia
    Gaia.ROW_LIMIT = -1
    cols = ("source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, bp_rp, parallax, parallax_error, ruwe, teff_gspphot, "
            "logg_gspphot, mh_gspphot, distance_gspphot, ag_gspphot, ebpminrp_gspphot")
    parts_dir = DATA / "_gaia_strips"
    parts_dir.mkdir(exist_ok=True)

    def strip(r0, r1):
        f = parts_dir / f"{r0:.4f}_{r1:.4f}.csv"
        if f.exists():
            return [f]
        q = (f"SELECT {cols} FROM gaiadr3.gaia_source WHERE ra >= {r0} AND ra < {r1} "
             "AND dec BETWEEN 36.4 AND 52.7 AND phot_g_mean_mag < 17.5")
        for k in range(4):
            try:
                t = Gaia.launch_job_async(q).get_results().to_pandas()
                break
            except Exception as e:  # noqa: BLE001
                print("   retry", r0, str(e)[:120], flush=True)
        else:
            raise RuntimeError(f"Gaia strip {r0}-{r1} failed")
        if len(t) >= 50000:
            m = (r0 + r1) / 2
            return strip(r0, m) + strip(m, r1)
        t.to_csv(f, index=False)
        print(f"  RA {r0:.3f}-{r1:.3f}: {len(t)}", flush=True)
        return [f]

    files = []
    r = 279.5
    while r < 301.9:
        files += strip(round(r, 4), round(r + 0.25, 4))
        r += 0.25
    d = pd.concat([pd.read_csv(f) for f in files]).drop_duplicates("source_id")
    o = DATA / f"gaia_field_{TODAY}.csv.gz"
    d.to_csv(o, index=False)
    print("wrote", o.name, len(d), "sources")


if __name__ == "__main__":
    DATA.mkdir(parents=True, exist_ok=True)
    steps = sys.argv[1:] or ["vizier", "koi", "kicstellar", "gaia"]
    for s in steps:
        globals()[s]()
