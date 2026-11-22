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
