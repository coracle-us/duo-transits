"""
Truth files for validation arms V2 and V3 of prereg_catalogue.yaml. No light curve is read.

V2: Hawthorn candidates with a solved period (data/catalog/status.csv, status "resolved").
V3: Pages Navarrete, arXiv:2611.03019 version 2 - the ancillary table of 187 TOIs
    (arxiv_2611.03019/source/anc/table2_full.csv) and the 12 published periods of its Table
    "Published periods recovered blindly" (main.tex; the 12 rows are kept in
    arxiv_2611.03019/published_periods.csv, the paper's source text is not in the repo), matched to our scored rows: a TOI row by its TOI
    number, any other row by TIC.

Inputs  : catalogue/defs.csv, data/catalog/status.csv, data/catalog/candidates.csv,
          catalogue/validation/arxiv_2611.03019/source/ (e-print v2, SHA-256 of the archive in SOURCE.md)
Outputs : catalogue/validation/v2_truth.csv, v3_truth.csv
Usage   : python catalogue/validation/build_truth.py
"""
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = HERE / "arxiv_2611.03019/source"


def published(tex):
    body = tex[tex.index(r"\label{tab:published}"):]
    body = body[body.index(r"\midrule"):body.index(r"\bottomrule")]
    rows = []
    for line in body.splitlines():
        c = [x.strip() for x in line.rstrip("\\ ").split("&")]
        if len(c) == 7 and re.fullmatch(r"[\d.]+", c[0]):
            rows.append(dict(TOI=float(c[0]), planet=c[1].replace(r"\,", " "), P_published_d=float(c[4]),
                             resolved_by=c[5], reference=re.sub(r"\\citet\{(.*)\}", r"\1", c[6])))
    return pd.DataFrame(rows)


def main():
    defs = pd.read_csv(ROOT / "catalogue/defs.csv")
    st = pd.read_csv(ROOT / "data/catalog/status.csv")
    st = st[st.status == "resolved"]
    h = defs[defs.catalogue == "hawthorn"].merge(
        st[["TIC", "period_claim_d", "tier", "source"]], left_on="tic", right_on="TIC")
    h[["row_id", "tic", "period_claim_d", "tier", "source"]].rename(
        columns=dict(period_claim_d="P_true_d")).to_csv(HERE / "v2_truth.csv", index=False)

    t = pd.read_csv(SRC / "anc/table2_full.csv")
    tex = SRC / "main.tex"                 # not in the repo; published_periods.csv is its extracted table
    pub = published(tex.read_text()) if tex.exists() else pd.read_csv(HERE / "arxiv_2611.03019/published_periods.csv")
    assert len(t) == 187 and len(pub) == 12
    t = t.merge(pub, on="TOI", how="left")
    toi = defs[defs.catalogue == "toi_single"].assign(TOI=lambda d: d.cand_id.str[3:].astype(float))
    a = toi.merge(t, on="TOI")[["row_id", "tic", "TOI", "State", "P_cat_d", "P_best_A_d", "w_best_A",
                                "P_published_d", "resolved_by", "reference"]].assign(match="toi")
    other = defs[defs.catalogue != "toi_single"].merge(t, left_on="tic", right_on="TIC")
    b = other[["row_id", "tic", "TOI", "State", "P_cat_d", "P_best_A_d", "w_best_A", "P_published_d",
               "resolved_by", "reference"]].assign(match="tic")
    out = pd.concat([a, b], ignore_index=True)
    out.to_csv(HERE / "v3_truth.csv", index=False)
    print("V2:", len(h), "rows;", "V3:", len(out), "rows,", out.State.value_counts().to_dict(),
          "; published:", int(out.P_published_d.notna().sum()))


if __name__ == "__main__":
    main()
