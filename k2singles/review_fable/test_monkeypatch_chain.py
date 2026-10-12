"""Does run_full.work's monkeypatch of alias_pipeline.load_lightcurves chain across tasks in a
reused worker process?  Reproduce the exact pattern of run_full.work with mocks, twice in one
process, and inspect what the second target's light curve contains.

Run: .venv/Scripts/python.exe k2singles/review_fable/test_monkeypatch_chain.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "k2singles" / "scripts"))
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
sys.path.insert(0, str(ROOT / "scripts"))

import alias_pipeline as ap  # noqa: E402
import pilot  # noqa: E402


def fake_tess(tic, max_sector=None, products=ap.PRODUCT_RANK):
    t = np.linspace(2500, 2510, 50)
    return pd.DataFrame({"time": t, "flux": 1.0, "err": np.nan, "sector": 44, "product": "SPOC_120",
                         "segment": 44000, "star": f"TESS-{tic}"})


ap.load_lightcurves = fake_tess
pilot.ap.load_lightcurves = fake_tess


def fake_k2(epic):
    t = np.linspace(140, 150, 30)
    return pd.DataFrame({"time": t, "flux": 1.0, "err": np.nan, "sector": 1005, "product": "K2PDC",
                         "segment": 1005000, "star": f"K2-{epic}"})


def work_like_run_full(tic, epic):
    """Verbatim structure of run_full.work lines 241-252 (minus run_target)."""
    k2 = fake_k2(epic)
    orig = ap.load_lightcurves

    def load_with_k2(t, max_sector=None, products=ap.PRODUCT_RANK):
        lc = orig(t, max_sector=max_sector, products=products)
        if k2 is None:
            return lc
        return pd.concat([lc, k2], ignore_index=True).sort_values("time").reset_index(drop=True)

    ap.load_lightcurves = load_with_k2
    pilot.ap.load_lightcurves = load_with_k2
    # what pilot.run_target will see:
    return pilot.ap.load_lightcurves(tic, products=pilot.PRODUCTS)


lc1 = work_like_run_full(111, 201000001)
lc2 = work_like_run_full(222, 201000002)
lc3 = work_like_run_full(333, 201000003)
print("target 1 light curve contents:", lc1.star.value_counts().to_dict())
print("target 2 light curve contents:", lc2.star.value_counts().to_dict())
print("target 3 light curve contents:", lc3.star.value_counts().to_dict())
assert "K2-201000001" in set(lc2.star), "no chaining (unexpected)"
print("\nCONFIRMED: in a reused worker process, every target after the first is searched together with the\n"
      "K2 light curves of all previous targets handled by that worker (orig = previous closure).")
