"""Download every TESS light curve of HIP 7666 used by RULES.md section 1 into data/lc/."""
import warnings
from pathlib import Path
import lightkurve as lk

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "data" / "lc"
OUT.mkdir(parents=True, exist_ok=True)

WANT = [  # (sector, author, exptime)
    (18, "TESS-SPOC", 1800), (58, "SPOC", 120), (85, "TESS-SPOC", 200),
    (18, "QLP", 1800), (58, "QLP", 200), (85, "QLP", 200),
]

res = lk.search_lightcurve("TIC 308953703", mission="TESS")
for sec, author, exp in WANT:
    m = res[(res.table["author"] == author) & (res.table["exptime"] == exp)
            & (res.table["mission"] == f"TESS Sector {sec:02d}")]
    if len(m) != 1:
        print(f"S{sec} {author} {exp}: {len(m)} matches, skipped")
        continue
    lc = m.download(download_dir=str(OUT), quality_bitmask="default")
    print(f"S{sec} {author} {exp}s: {len(lc)} cadences, file {Path(lc.meta['FILENAME']).name}")
