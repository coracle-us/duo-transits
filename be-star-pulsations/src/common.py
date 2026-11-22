"""Shared paths and constants for the Be star pulsation project.

Set BE_RUN=full for the full run (outputs in results/full/, figures full_*.png);
the default is the pilot (results/pilot/, figures pilot_*.png)."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
RUN = os.environ.get("BE_RUN", "pilot")
PILOT = RESULTS / RUN   # this run's output directory (name kept from the pilot)
FIGURES = ROOT / "figures"
CACHE = Path.home() / "astro-cache" / "be-star-pulsations"
CAT_DIR = CACHE / "catalogue"
WISE_DIR = CACHE / "wise"
TESS_DIR = CACHE / "tess"
SPEC_DIR = CACHE / "spectra"

SEED = 20261104
BTJD_TO_MJD = 2457000.0 - 2400000.5  # MJD = BTJD + 56999.5

for d in (RESULTS, PILOT, FIGURES, CAT_DIR, WISE_DIR, TESS_DIR, SPEC_DIR):
    d.mkdir(parents=True, exist_ok=True)


def figpath(base):
    """figures/{run}_{base}.png; the pilot's spline check keeps its original name."""
    if base == "spline_check" and RUN == "pilot":
        return FIGURES / "spline_check.png"
    return FIGURES / f"{RUN}_{base}.png"
