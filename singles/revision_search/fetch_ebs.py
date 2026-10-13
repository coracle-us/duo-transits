"""
Fetch TESS light curves (SPOC, TESS-SPOC, QLP) for a few known Prša et al. (2022) TESS eclipsing
binaries, used as positive controls for the unconstrained dip search (search_any.py --ebs).

Reuses scripts/fetch_data._fetch_one, but redirects its output to
singles/revision_search/data/eb_lightcurves/<TIC>/ (gitignored) so data/raw/ is not touched,
and logs to stdout only.

Usage: python singles/revision_search/fetch_ebs.py [TIC ...]
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import fetch_data as fd  # noqa: E402

# Prša+2022 EBs that were themselves uploaded as single-transit CTOIs (dropped from the sample),
# chosen for detached morphology (< 0.3) and periods 8-29 d.
EB_TICS = [92349924, 77320223, 219771660, 362104293]

fd.LC_DIR = HERE / "data" / "eb_lightcurves"
fd.LC_DIR.mkdir(parents=True, exist_ok=True)
fd.log = lambda msg: print(msg, flush=True)
fd.AUTHORS = ("SPOC", "TESS-SPOC", "QLP")

if __name__ == "__main__":
    tics = [int(x) for x in sys.argv[1:]] or EB_TICS
    for i, tic in enumerate(tics, 1):
        print(fd._fetch_one(i, tic), flush=True)
