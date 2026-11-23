"""Round 3: the draft must parse, and nothing of round 3 may exist before it is approved."""
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parents[1]


def test_round3_prereg_parses_and_is_gated():
    cfg = yaml.safe_load(open(HERE / "prereg_round3.yaml"))
    assert {"tuning", "test_set", "pass_criteria", "outcomes"} <= set(cfg)
    text = (HERE / "APPROVALS.md").read_text()
    if "## Round 3" not in text:
        assert str(cfg["status"]).startswith("draft")
        assert not (HERE / "round3").exists()
