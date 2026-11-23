"""Round 3: the approved pre-registration must not change, and the test set is gated on the freeze."""
import hashlib
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parents[1]


def _round3_text():
    text = (HERE / "APPROVALS.md").read_text()
    return text[text.index("## Round 3"):]


def test_round3_prereg_is_frozen():
    cfg = yaml.safe_load(open(HERE / "prereg_round3.yaml"))
    assert str(cfg["status"]).startswith("approved")
    rows = [ln for ln in _round3_text().splitlines() if ln.startswith(f"| {cfg['version']} |")]
    assert len(rows) == 1
    assert hashlib.sha256((HERE / "prereg_round3.yaml").read_bytes()).hexdigest() in rows[0]


def test_round3_test_set_is_not_loaded_before_the_freeze():
    if "Step 2 of 2" not in _round3_text():
        assert not (HERE / "round3" / "test").exists()
