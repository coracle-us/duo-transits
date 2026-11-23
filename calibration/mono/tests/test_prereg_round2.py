"""The approved round-2 pre-registration must parse and must not change without a new approval."""
import hashlib
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parents[1]


def test_round2_prereg_is_frozen():
    cfg = yaml.safe_load(open(HERE / "prereg_round2.yaml"))
    assert str(cfg["status"]).startswith("approved")
    text = (HERE / "APPROVALS.md").read_text()
    rows = [ln for ln in text[text.index("## Round 2"):].splitlines() if ln.startswith(f"| {cfg['version']} |")]
    assert len(rows) == 1
    assert hashlib.sha256((HERE / "prereg_round2.yaml").read_bytes()).hexdigest() in rows[0]


def test_round2_test_set_is_not_loaded_before_the_freeze():
    """No round-2 test-set result may exist until APPROVALS.md records the step-2 freeze."""
    text = (HERE / "APPROVALS.md").read_text()
    if "Step 2 of 2" not in text:
        assert not (HERE / "round2" / "test").exists()
