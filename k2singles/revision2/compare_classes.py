"""
AAS82225 revision, literature search: compare the two independent classifications of the Table 1 stars' hits
(PLAN.md, "Independent check"). A = table1_classes_A.csv, B = table1_classes_B.csv (walled-off agent).
A row disagrees if the class differs, or both are PERIOD/ESTIMATE with periods differing by > 1%.
LISTED vs ANALYSIS disagreements are reported separately: neither carries a period, so they cannot change a verdict.
Usage:  python k2singles/revision2/compare_classes.py
Output: stdout, data/literature/table1_class_compare.csv
"""
from pathlib import Path

import pandas as pd

LIT = Path(__file__).resolve().parent / "data" / "literature"
NOPER = {"LISTED", "ANALYSIS"}


def main():
    a = pd.read_csv(LIT / "table1_classes_A.csv")
    b = pd.read_csv(LIT / "table1_classes_B.csv")
    assert len(a) == len(b) == len(a.merge(b, on=["EPIC", "bibcode"])), "row sets differ"
    m = a.merge(b, on=["EPIC", "bibcode"], suffixes=("_A", "_B"))

    def kind(r):
        ca, cb = r["class_A"], r["class_B"]
        if ca == cb:
            if ca in ("PERIOD", "ESTIMATE"):
                pa, pb = float(r.period_d_A), float(r.period_d_B)
                return "agree" if abs(pa - pb) <= 0.01 * pa else "period differs"
            return "agree"
        if {ca, cb} <= NOPER:
            return "LISTED/ANALYSIS only"
        return "class differs"

    m["compare"] = m.apply(kind, axis=1)
    m.to_csv(LIT / "table1_class_compare.csv", index=False)
    print(m["compare"].value_counts().to_dict())
    print("A:", a["class"].value_counts().to_dict())
    print("B:", b["class"].value_counts().to_dict())
    bad = m[m["compare"].isin(["class differs", "period differs"])]
    print(bad[["EPIC", "bibcode", "class_A", "period_d_A", "class_B", "period_d_B", "notes_B"]].to_string(index=False))


if __name__ == "__main__":
    main()
