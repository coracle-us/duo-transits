"""The catalogue pre-registration parses, and the counts it states match the frozen candidate list."""
import re
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
PRE = yaml.safe_load((ROOT / "catalogue/prereg_catalogue.yaml").read_text())
CAND = pd.read_csv(ROOT / "gaia_dr4/candidates_frozen.csv")


def _nums(s):
    return [int(x.replace(",", "")) for x in re.findall(r"\d[\d,]*", s)]


def test_scored_counts_match_frozen_list():
    for cat, text in PRE["scope"]["scored"].items():
        g = CAND[CAND.catalogue == cat]
        assert _nums(text)[:2] == [len(g), g.tic.nunique()], cat


def test_not_scored_counts():
    assert _nums(PRE["scope"]["not_scored"]["k2"])[0] == (CAND.catalogue == "k2").sum()
    own = CAND.catalogue.str.startswith("own").sum()
    assert _nums(PRE["scope"]["not_scored"]["own_ctoi_and_own_salinas"])[0] == own


def test_solved_hawthorn_count():
    st = pd.read_csv(ROOT / "data/catalog/status.csv")
    assert (st.status == "resolved").sum() == 26


def test_prereg_frozen_by_hash():
    """Version 1.0 as approved (catalogue/APPROVALS.md). Any change is a new version and needs Jon."""
    import hashlib
    h = hashlib.sha256((ROOT / "catalogue/prereg_catalogue.yaml").read_bytes()).hexdigest()
    assert h == "1177568b08db3431af02d488cab5bcd7fbee8db7306ebd6d80dab73adfb75d72"
    assert h in (ROOT / "catalogue/APPROVALS.md").read_text()
