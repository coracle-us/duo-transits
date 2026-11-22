"""The pre-registration must parse, and its stated sizes must match the frozen files."""
import re
from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent


def test_prereg_parses_and_has_every_threshold_the_code_reads():
    cfg = yaml.safe_load(open(HERE / "prereg_dr4.yaml"))
    for k in ("ruwe_max", "nsig", "dt_err_d", "max_window", "period_rel_floor", "harmonics"):
        assert k in cfg["thresholds"]
    assert all(isinstance(x, str) for x in cfg["human_review"])
    for key in ("candidates", "hosts"):
        assert (ROOT / cfg["inputs"][key]).exists()


def test_prereg_row_counts_match_the_frozen_files():
    text = (HERE / "prereg_dr4.yaml").read_text()
    n_cand = int(re.search(r"candidates_frozen\.csv\s+# (\d+) rows", text).group(1))
    n_host = int(re.search(r"hosts_frozen\.csv\s+# (\d+) stars", text).group(1))
    assert n_cand == len(pd.read_csv(HERE / "candidates_frozen.csv"))
    assert n_host == len(pd.read_csv(HERE / "hosts_frozen.csv"))
