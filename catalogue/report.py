"""
Score the validation arms and summarise the catalogue run (prereg_catalogue.yaml, validation, output).

V2: solved Hawthorn periods - not excluded by step 1 in at least 25 of 26; not excluded by the
    duotransit alias table in at least 25 of 26. A solved period that is no alias of the two catalogued
    times ("not an alias") is counted separately; it is not an exclusion by the table.
V3: periods published elsewhere (arXiv:2611.03019) on our stars - at most 5% excluded by step 1,
    none if fewer than 20. Their other resolved periods are reported, not judged.
V4: every period this project stated before, listed with this run's verdict (no criterion).

Inputs  : <out>/catalogue.csv, period_checks_tables.csv, vetting/verdicts.csv, vetter_check/vetter_check_score.md
Outputs : <out>/summary.md, <out>/validation.csv
Usage   : python catalogue/report.py [--out catalogue/results/v0]
"""
import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
COMPLETENESS = ("In round 3 the rule found 54% and 56% of clean repeats, and almost none below a reference SNR "
                "of 6. \"No repeat\" therefore does not mean that no repeat exists.")


def validation(pc):
    """One line per criterion: (arm, criterion, k, n, limit text, verdict)."""
    out = []
    v2 = pc[pc.source == "v2_solved"]
    n = len(v2)
    k1 = int(v2.step1_excluded.sum())
    st = v2.get("alias_status_duo", pd.Series(dtype=str))
    k2, na = int((st == "excluded").sum()), int((st == "not an alias").sum())
    out.append(("V2", "solved period excluded by step 1", k1, n, "at most 1 of 26", "met" if n == 26 and k1 <= 1 else "MISSED"))
    out.append(("V2", "solved period excluded by the duotransit alias table", k2, n, "at most 1 of 26",
                "met" if n == 26 and k2 <= 1 else "MISSED"))
    out.append(("V2", "solved period not an alias of the two catalogued times", na, n, "reported", ""))
    v3 = pc[pc.source == "v3_published"]
    k3, n3 = int(v3.step1_excluded.sum()), len(v3)
    ok = k3 == 0 if n3 < 20 else k3 <= 0.05 * n3
    out.append(("V3", "published period excluded by step 1", k3, n3, "none (fewer than 20)" if n3 < 20 else "at most 5%",
                "met" if ok else "MISSED"))
    vb = pc[pc.source == "v3_best"]
    out.append(("V3", "their resolved period excluded by step 1", int(vb.step1_excluded.sum()), len(vb), "reported", ""))
    for src, label in (("own", "V4"), ("stated", "stated")):
        x = pc[pc.source == src]
        out.append((label, f"{src} period excluded by step 1", int(x.step1_excluded.sum()), len(x), "reported", ""))
    return pd.DataFrame(out, columns=["arm", "criterion", "k", "n", "limit", "verdict"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "results/v0"))
    out = Path(ap.parse_args().out)
    c = pd.read_csv(out / "catalogue.csv")
    pc = pd.read_csv(out / "period_checks_tables.csv")
    val = validation(pc)
    val.to_csv(out / "validation.csv", index=False)
    L = ["# Catalogue run: summary", "", COMPLETENESS, "",
         "Every alias table is conditional on the second event being the same body as the reference.", "",
         "## Rows", "", "| Catalogue | Rows | Scored | No light curve | Reference not in data | Reference measured | "
         "Median ln P excluded (measured) | Passing events | Accepted events | Resolved | Narrowed | All excluded |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for cat, g in c.groupby("catalogue"):
        ok = g[g.status == "ok"]
        m = ok[ok.ref_measured.astype(bool)]
        L.append(f"| {cat} | {len(g)} | {len(ok)} | {(g.status == 'no light curve').sum()} | "
                 f"{(g.status == 'reference not in data').sum()} | {len(m)} | "
                 f"{m.frac_lnP_excluded.median():.2f} | {int(ok.n_passing.sum())} | {int(ok.n_accepted.sum())} | "
                 f"{(ok.call == 'resolved').sum()} | {(ok.call == 'narrowed').sum()} | {(ok.call == 'all excluded').sum()} |")
    other = c[~c.status.isin(["ok", "no light curve", "reference not in data"])]
    L += ["", f"Rows with another status (errors): {len(other)}.", "", "## Validation", "",
          "| Arm | Criterion | Count | Limit | Verdict |", "|---|---|---|---|---|"]
    L += [f"| {r.arm} | {r.criterion} | {r.k} of {r.n} | {r.limit} | {r.verdict} |" for r in val.itertuples()]
    v1 = out / "vetter_check/vetter_check_score.md"
    L += ["", v1.read_text() if v1.exists() else "Vetter check (V1): no score file."]
    (out / "summary.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
