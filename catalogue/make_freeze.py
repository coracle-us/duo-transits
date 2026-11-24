"""
Write or check catalogue/FREEZE.json: SHA-256 of every file frozen at step 2 of the catalogue run.

Inputs  : the files listed in FILES, and calibration/mono/round3/FREEZE3.json
Outputs : catalogue/FREEZE.json (--write); exit status 1 if any hash differs (--check)
Usage   : python catalogue/make_freeze.py --write | --check
"""
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FILES = """prereg_catalogue.yaml defs.csv defs_not_scored.csv tic_info.csv check_periods.csv
validation/v2_truth.csv validation/v3_truth.csv validation/arxiv_2611.03019/published_periods.csv
validation/arxiv_2611.03019/source/anc/table2_full.csv build_defs.py build_check_periods.py
validation/build_truth.py catrun.py run_catalogue.py sheets.py make_sheets.py VETTING.md
merge_verdicts.py report.py""".split()


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def current():
    h = {f"catalogue/{f}": sha(ROOT / "catalogue" / f) for f in FILES}
    f3 = json.loads((ROOT / "calibration/mono/round3/FREEZE3.json").read_text())["sha256"]
    h["calibration/mono/round3/FREEZE3.json"] = sha(ROOT / "calibration/mono/round3/FREEZE3.json")
    return h, f3


def check():
    """List of files whose hash differs from FREEZE.json or from round 3's freeze."""
    fz = json.loads((ROOT / "catalogue/FREEZE.json").read_text())
    h, f3 = current()
    bad = [k for k, v in fz["sha256"].items() if h.get(k) != v]
    bad += [k for k, v in f3.items() if sha(ROOT / k) != v]
    return bad


if __name__ == "__main__":
    if sys.argv[1:] == ["--write"]:
        h, _ = current()
        out = dict(frozen="2026-11-23", approved_by="Jon McBee",
                   note="step 2 of 2 of prereg_catalogue.yaml 1.0; every file of round3/FREEZE3.json is "
                        "frozen too and checked against that file",
                   sha256=h)
        (ROOT / "catalogue/FREEZE.json").write_text(json.dumps(out, indent=2) + "\n")
    else:
        bad = check()
        print("freeze ok" if not bad else "CHANGED: " + ", ".join(bad))
        sys.exit(1 if bad else 0)
