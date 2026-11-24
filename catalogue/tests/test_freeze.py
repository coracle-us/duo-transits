import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import make_freeze  # noqa: E402


def test_frozen_files_unchanged():
    assert make_freeze.check() == []


def test_step2_recorded():
    assert "Step 2 of 2" in (make_freeze.ROOT / "catalogue/APPROVALS.md").read_text()
