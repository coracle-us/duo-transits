"""
Fetch the raw inputs for the host-rotation retest into data/ (git-ignored), with the download date in each file name.

  VizieR (IDIA mirror, asu-tsv; files <stem>_<table>_<date>.csv):
    J/ApJ/775/L11   McQuillan, Mazeh & Aigrain 2013, KOI rotation periods   -> mma2013_<date>.csv
    J/ApJS/211/24   McQuillan, Mazeh & Aigrain 2014, table1 (periods) + table2 (no period)
    J/ApJS/244/21   Santos et al. 2019 (M and K dwarfs)                     -> santos2019_<date>.csv
    J/ApJS/255/17   Santos et al. 2021 (F and G stars)                      -> santos2021_<date>.csv
  NASA Exoplanet Archive TAP:
    cumulative (KOI table: kepid, koi disposition)                          -> koi_cumulative_<date>.csv
  Gaia DR3 (ESA Gaia archive TAP, every source within 2 arcsec of each KIC position)
                                                                            -> gaia_xmatch_<date>.csv
Usage:  .venv/bin/python fetch.py [step ...]      steps: vizier koi gaia (default: all; a step is skipped if its file exists)
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
    fs = sorted(DATA.glob(f"{stem}_*.csv"))
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


def gaia():
    if latest("gaia_xmatch"):
        print("gaia exists:", latest("gaia_xmatch").name)
        return
    # CDS XMatch was down on 2026-11-24 (as was gaia-kepler.fun's TLS certificate), so the match is a 2" cone search
    # in the ESA Gaia archive (TAP, uploaded KIC positions, anonymous). Every Gaia source within 2" is kept;
    # match.py picks the counterpart after moving Gaia positions back to J2000 with the proper motion.
    from astropy.table import Table
    from astroquery.gaia import Gaia
    pos = []
    for f in DATA.glob("*_table*_*.csv"):
        d = pd.read_csv(f, usecols=lambda c: c in ("KIC", "_RAJ2000", "_DEJ2000"))
        pos.append(d)
    k = pd.read_csv(latest("koi_cumulative"))[["kepid", "ra", "dec"]].rename(
        columns={"kepid": "KIC", "ra": "_RAJ2000", "dec": "_DEJ2000"})
    pos = pd.concat(pos + [k]).dropna().drop_duplicates("KIC")
    pos.columns = ["kic", "ra", "dec"]
    print("positions:", len(pos))
    cols = ("g.source_id, g.ra AS gra, g.dec AS gdec, g.pmra, g.pmdec, g.phot_g_mean_mag, g.bp_rp, g.parallax, "
            "g.parallax_error, g.ruwe, g.teff_gspphot, g.logg_gspphot, g.mh_gspphot, g.distance_gspphot, "
            "g.ag_gspphot, g.ebpminrp_gspphot")
    q = (f"SELECT u.kic, u.ra, u.dec, {cols} FROM tap_upload.kic AS u JOIN gaiadr3.gaia_source AS g "
         "ON 1 = CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', u.ra, u.dec, 2.0/3600.0))")
    parts = []
    tmp = DATA / "_upload.xml"
    step = 20000
    for i in range(0, len(pos), step):
        part = DATA / f"_gaia_part_{i:06d}.csv"
        if part.exists():
            parts.append(pd.read_csv(part))
            continue
        Table.from_pandas(pos.iloc[i:i + step].reset_index(drop=True)).write(tmp, format="votable", overwrite=True)
        job = Gaia.launch_job_async(q, upload_resource=str(tmp), upload_table_name="kic")
        r = job.get_results().to_pandas()
        r.to_csv(part, index=False)
        print(f"  rows {i}-{i + step}: {len(r)} Gaia sources within 2 arcsec", flush=True)
        parts.append(r)
    o = DATA / f"gaia_xmatch_{TODAY}.csv"
    pd.concat(parts).to_csv(o, index=False)
    for p in DATA.glob("_gaia_part_*.csv"):
        p.unlink()
    tmp.unlink(missing_ok=True)
    print("wrote", o.name)


if __name__ == "__main__":
    DATA.mkdir(parents=True, exist_ok=True)
    steps = sys.argv[1:] or ["vizier", "koi", "gaia"]
    for s in steps:
        globals()[s]()
