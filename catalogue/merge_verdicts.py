"""
Join the two vetters' verdict files to the sheet key. An event is accepted only if both accept.

Inputs  : <dir>/key.csv (sheet_id, event_id, ...), <dir>/verdicts_1.csv, <dir>/verdicts_2.csv
          (sheet_id, verdict, reason, note - catalogue/VETTING.md)
Outputs : <dir>/verdicts.csv (key columns, vetter_1, reason_1, note_1, vetter_2, reason_2, note_2, accepted)
          For the vetter check (key has at_truth) also <dir>/vetter_check_score.md
Usage   : python catalogue/merge_verdicts.py catalogue/results/v0/vetting
"""
import sys
from pathlib import Path

import pandas as pd

VERDICTS = {"accept", "reject", "cannot_tell"}


def merge(d):
    d = Path(d)
    key = pd.read_csv(d / "key.csv")
    key = key[key.sheet_id.notna() & (key.sheet_id != "")]
    for i in (1, 2):
        v = pd.read_csv(d / f"verdicts_{i}.csv", dtype=str).fillna("")
        assert set(v.verdict) <= VERDICTS, f"vetter {i}: unknown verdict {set(v.verdict) - VERDICTS}"
        assert v.sheet_id.is_unique and set(v.sheet_id) == set(key.sheet_id), f"vetter {i}: sheets do not match the key"
        key = key.merge(v.rename(columns=dict(verdict=f"vetter_{i}", reason=f"reason_{i}", note=f"note_{i}")),
                        on="sheet_id", validate="1:1")
    key["accepted"] = (key.vetter_1 == "accept") & (key.vetter_2 == "accept")
    key.to_csv(d / "verdicts.csv", index=False)
    return key


def score_v1(k):
    """prereg V1: each vetter accepts at least 90% of the true events."""
    t, f = k[k.at_truth.astype(bool)], k[~k.at_truth.astype(bool)]
    lines = ["# Vetter check (V1)", "", f"{len(t)} true events, {len(f)} false events.", "",
             "| | Vetter 1 | Vetter 2 | Both |", "|---|---|---|---|"]
    a1, a2 = (t.vetter_1 == "accept").sum(), (t.vetter_2 == "accept").sum()
    lines.append(f"| True events accepted | {a1} of {len(t)} | {a2} of {len(t)} | {int(t.accepted.sum())} of {len(t)} |")
    r1, r2 = (f.vetter_1 == "reject").sum(), (f.vetter_2 == "reject").sum()
    lines.append(f"| False events rejected | {r1} of {len(f)} | {r2} of {len(f)} | "
                 f"{int(((f.vetter_1 == 'reject') & (f.vetter_2 == 'reject')).sum())} of {len(f)} |")
    lines.append(f"| False events accepted | {(f.vetter_1 == 'accept').sum()} of {len(f)} | "
                 f"{(f.vetter_2 == 'accept').sum()} of {len(f)} | {int(f.accepted.sum())} of {len(f)} |")
    ok = len(t) > 0 and min(a1, a2) >= 0.9 * len(t)
    lines += ["", f"Criterion (each vetter accepts at least 90% of the true events): **{'met' if ok else 'MISSED'}**."]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    k = merge(sys.argv[1])
    print(len(k), "events;", int(k.accepted.sum()), "accepted by both")
    if "at_truth" in k:
        Path(sys.argv[1], "vetter_check_score.md").write_text(score_v1(k))
