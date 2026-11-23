"""
AAS82225 revision, literature search: compare the two independent classifications of the Table 1 stars' hits
(PLAN.md, "Independent check"). A = table1_classes_A.csv, B = table1_classes_B.csv (walled-off agent).
A row disagrees if the class differs, or both are PERIOD/ESTIMATE with periods differing by > 1%.
LISTED vs ANALYSIS disagreements are reported separately: neither carries a period, so they cannot change a verdict.
With --vizier: the same for the fixed-list VizieR rows (table1_vizier_classes_A/B.csv, keyed by catalogue table; a
pair may carry several periods separated by ';', which must match one-to-one within 1%).
Usage:  python k2singles/revision2/compare_classes.py [--vizier]
Output: stdout, data/literature/table1_class_compare.csv (or table1_vizier_class_compare.csv)
"""
import sys
from pathlib import Path

import pandas as pd

LIT = Path(__file__).resolve().parent / "data" / "literature"
NOPER = {"LISTED", "ANALYSIS"}


def periods(x):
    return sorted(float(v) for v in str(x).split(";") if v.strip() not in ("", "nan"))


def same_periods(pa, pb):
    pa, pb = periods(pa), periods(pb)
    return len(pa) == len(pb) and all(abs(x - y) <= 0.01 * x for x, y in zip(pa, pb))


def main():
    viz = "--vizier" in sys.argv
    stem, key = ("table1_vizier_classes", "catalog") if viz else ("table1_classes", "bibcode")
    a = pd.read_csv(LIT / f"{stem}_A.csv")
    b = pd.read_csv(LIT / f"{stem}_B.csv")
    assert len(a) == len(b) == len(a.merge(b, on=["EPIC", key])), "row sets differ"
    m = a.merge(b, on=["EPIC", key], suffixes=("_A", "_B"))

    def kind(r):
        ca, cb = r["class_A"], r["class_B"]
        if ca == cb:
            if ca in ("PERIOD", "ESTIMATE"):
                return "agree" if same_periods(r.period_d_A, r.period_d_B) else "period differs"
            return "agree"
        if {ca, cb} <= NOPER:
            return "LISTED/ANALYSIS only"
        return "class differs"

    m["compare"] = m.apply(kind, axis=1)
    m.to_csv(LIT / ("table1_vizier_class_compare.csv" if viz else "table1_class_compare.csv"), index=False)
    print(m["compare"].value_counts().to_dict())
    print("A:", a["class"].value_counts().to_dict())
    print("B:", b["class"].value_counts().to_dict())
    bad = m[m["compare"].isin(["class differs", "period differs"])]
    note = "note_B" if viz else "notes_B"
    print(bad[["EPIC", key, "class_A", "period_d_A", "class_B", "period_d_B", note]].to_string(index=False))
    if viz:
        la = m[m["compare"] == "LISTED/ANALYSIS only"]
        print(la[["EPIC", key, "class_A", "class_B"]].to_string(index=False))


if __name__ == "__main__":
    main()
