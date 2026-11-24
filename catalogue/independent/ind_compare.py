"""Compare the independent search stage (ind_rows.csv, ind_events.csv) with the primary
(catalogue/results/v0/rows_search.csv, events.csv). Prints the numbers COMPARE.md quotes.
Usage: .venv/bin/python -W ignore catalogue/independent/ind_compare.py"""
from pathlib import Path
import numpy as np, pandas as pd
V = Path(__file__).resolve().parents[1] / "results/v0"
I = V / "independent"
P, R = pd.read_csv(V / "rows_search.csv"), pd.read_csv(I / "ind_rows.csv")
PE, RE = pd.read_csv(V / "events.csv"), pd.read_csv(I / "ind_events.csv")
assert list(P.row_id) == list(R.row_id), "row order differs"
m = P.merge(R, on="row_id", suffixes=("_p", "_i"))
dis = []
def cmp(col, kind, rtol=1e-9):
    a, b = m[col + "_p"], m[col + "_i"]
    if kind == "num":
        a, b = a.astype(float), b.astype(float)
        same = (a.isna() & b.isna()) | (np.abs(a - b) <= rtol * np.maximum(np.abs(a), np.abs(b)))
    else:
        same = (a.astype(str) == b.astype(str)) | (a.isna() & b.isna())
    print(f"{col}: {int(same.sum())} of {len(m)} agree; both empty {int((a.isna() & b.isna()).sum())}")
    for r in m[~same].itertuples():
        dis.append(dict(row_id=r.row_id, column=col, primary=getattr(r, col + "_p"), independent=getattr(r, col + "_i")))
print("status primary:", P.status.value_counts().to_dict())
cmp("status", "str")
m["ref_measured_p"] = m.ref_measured_p.map(lambda v: np.nan if pd.isna(v) else float(str(v) == "True"))
m["ref_measured_i"] = m.ref_measured_i.map(lambda v: np.nan if pd.isna(v) else float(str(v) == "True"))
for c in ("ref_measured", "ref_sector", "ref_snr", "depth_ref", "dur_ref_d", "frac_lnP_excluded",
          "p_false_exclusion", "n_passing", "n_events", "n_sectors", "n_points", "baseline_d"):
    cmp(c, "num")
m["sectors_i"] = m.sectors_i.str.replace(";", " ")   # separator only
m.loc[m.status_p != "ok", ["sectors_p", "sectors_i"]] = np.nan   # the primary lists sectors on ok rows only
cmp("sectors", "str")
m["p_false_exclusion_i"] = m.p_false_exclusion_i.round(4)      # the table's rounded p column
m = m.rename(columns={"p_false_exclusion_p": "p_table_p", "p_false_exclusion_i": "p_table_i"})
print("after rounding the independent p to the table's 4 decimals:")
cmp("p_table", "num")
# events: all, and the passing set per row (sector equal, tc within 1e-6 d)
def match(A, B):
    """one-to-one match of event tables on sector and |tc| < 1e-6; returns unmatched of each."""
    ua, ub = [], B.copy()
    ub["used"] = False
    for r in A.itertuples():
        c = ub[(~ub.used) & (ub.row_id == r.row_id) & (ub.sector == r.sector) & ((ub.tc - r.tc).abs() < 1e-6)]
        if len(c):
            ub.loc[c.index[0], "used"] = True
        else:
            ua.append(r)
    return ua, ub[~ub.used]
k = ["row_id", "sector", "tc"]
A, B = PE.sort_values(k).reset_index(drop=True), RE.sort_values(k).reset_index(drop=True)
print("events: primary", len(A), "independent", len(B))
pos = (A.row_id == B.row_id) & (A.sector == B.sector) & ((A.tc - B.tc).abs() < 1e-6)
print("  events matched one to one (row, sector, tc within 1e-6 d):", int(pos.sum()), "of", len(A))
print("  passes agree:", int((A.passes.astype(bool) == B.passes.astype(bool)).sum()), "of", len(A))
for c in ("tc", "snr", "sec_box_snr", "depth", "depth_err", "dur_d", "frac", "edge_dist_d", "anti_snr", "depth_ratio", "dur_ratio"):
    a, b = A[c].astype(float), B[c].astype(float)
    same = (a.isna() & b.isna()) | (np.abs(a - b) <= 1e-9 * np.maximum(np.abs(a), np.abs(b)))
    print(f"  event {c}: {int(same.sum())} of {len(A)} agree; both empty {int((a.isna() & b.isna()).sum())}")
pp, ip = PE[PE.passes.astype(bool)], RE[RE.passes.astype(bool)]
ua, ub = match(pp, ip)
print("passing events: primary", len(pp), "independent", len(ip), "primary-only", len(ua), "independent-only", len(ub))
rows_p = set(pp.row_id) | set(ip.row_id)
bad = set(r.row_id for r in ua) | set(ub.row_id)
print("rows with a passing event on either side:", len(rows_p), "; rows whose passing sets agree:",
      len(m) - len(bad), "of", len(m))
for r in ua:
    dis.append(dict(row_id=r.row_id, column="passing event", primary=f"S{r.sector} tc={r.tc}", independent="absent"))
for r in ub.itertuples():
    dis.append(dict(row_id=r.row_id, column="passing event", primary="absent", independent=f"S{r.sector} tc={r.tc}"))
D = pd.DataFrame(dis, columns=["row_id", "column", "primary", "independent"])
D.to_csv(I / "compare_disagreements.csv", index=False)
print(len(D), "disagreements"); print(D.groupby("column").size().to_string() if len(D) else "")
