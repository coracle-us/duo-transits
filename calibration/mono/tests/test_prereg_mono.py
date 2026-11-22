"""The approved single-transit pre-registration must parse and must not change without a new approval."""
import hashlib
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parents[1]


def test_approved_prereg_is_frozen():
    cfg = yaml.safe_load(open(HERE / "prereg_mono_calibration.yaml"))
    assert cfg["status"] == "approved"
    rows = [ln for ln in (HERE / "APPROVALS.md").read_text().splitlines() if ln.startswith("| ")]
    row = [ln for ln in rows if ln.startswith(f"| {cfg['version']} |")]
    assert len(row) == 1
    assert hashlib.sha256((HERE / "prereg_mono_calibration.yaml").read_bytes()).hexdigest() in row[0]


def test_prereg_has_the_fields_the_harness_reads():
    cfg = yaml.safe_load(open(HERE / "prereg_mono_calibration.yaml"))
    b = cfg["arm_B_injections"]
    assert b["draws_per_host"] == 4 and "20261123" in b["hosts"]
    assert "20261122" in cfg["arm_A_confirmed_planets"]["sample"]
    assert cfg["pipeline"]["epoch_rule"].startswith("half3s")


def test_drawn_samples_match_the_prereg():
    import pandas as pd
    a = pd.read_csv(HERE / "arm_a_sample.csv")
    b = pd.read_csv(HERE / "arm_b_hosts.csv")
    assert len(a) == 120 and a.tic.is_unique and (a.depth_bin.value_counts() == 40).all()
    assert a.period_d.between(10, 100).all() and a.disposition.isin(["CP", "KP"]).all()
    assert len(b) == 150 and b.tic.is_unique and (b.tmag_bin.value_counts() == 50).all()
    hosts = set(pd.read_csv(HERE.parent.parent / "gaia_dr4/hosts_frozen.csv").tic)
    assert set(b.tic) <= hosts and not (set(a.tic) & hosts)
