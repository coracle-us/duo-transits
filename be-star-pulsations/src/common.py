"""Shared paths and constants for the Be star pulsation pilot."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PILOT = RESULTS / "pilot"
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
