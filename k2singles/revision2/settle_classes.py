"""
AAS82225 revision, literature search: the settled classes for the Table 1 stars' hits (PLAN.md, "Independent check":
disagreements are settled by reading the paper, and both classes are recorded). Input: table1_class_compare.csv
(compare_classes.py). Each settlement below was decided by reading the source named in it.
Usage:  python k2singles/revision2/settle_classes.py
Output: data/literature/table1_classes_final.csv
"""
from pathlib import Path

import pandas as pd

LIT = Path(__file__).resolve().parent / "data" / "literature"
SETTLED = {
    (202071902, "2021RAA....21..292W"): ("PERIOD", 8.6973,
        "B correct: online Table 3 (VizieR J/other/RAA/21.292/table3, 0.12\" match) lists Type EA, Plc 8.6973 d "
        "from AAVSO VSX; A had searched the PDF text, which prints only a sample of the table"),
}


def main():
    m = pd.read_csv(LIT / "table1_class_compare.csv")
    m["class_final"], m["period_final_d"], m["settled"] = m["class_A"], m["period_d_A"], ""
    for (e, b), (c, p, why) in SETTLED.items():
        i = (m["EPIC"] == e) & (m["bibcode"] == b)
        assert i.sum() == 1
        m.loc[i, "class_final"], m.loc[i, "period_final_d"], m.loc[i, "settled"] = c, p, why
    m.loc[m["compare"] == "LISTED/ANALYSIS only", "settled"] = "no period either way; A's class kept"
    left = m[(m["compare"].isin(["class differs", "period differs"])) & (m["settled"] == "")]
    assert left.empty, f"unsettled disagreements: {left[['EPIC', 'bibcode']].values.tolist()}"
    m[["EPIC", "bibcode", "class_A", "class_B", "class_final", "period_final_d", "compare", "settled",
       "notes_A", "notes_B"]].to_csv(LIT / "table1_classes_final.csv", index=False)
    print(m["class_final"].value_counts().to_dict())


if __name__ == "__main__":
    main()
