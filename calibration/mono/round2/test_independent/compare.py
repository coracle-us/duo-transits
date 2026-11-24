#!/usr/bin/env python3
"""Step 4: per-trial comparison of the independent rows with the primary rows (read only after
ind_*_rows.csv, ind_score.md and AMBIGUITIES.md were on disk; see SHA256_before_compare.txt).
Prints the per-column agreement table and the disagreeing trials as markdown."""
from pathlib import Path
import numpy as np
import pandas as pd

H = Path(__file__).resolve().parent
BOOL = ["truth_excluded", "false_event", "repeat_detected", "truth_in_alias_set", "ref_measured"]
STR = ["call"]
NUM = ["frac_lnP_excluded", "ref_snr", "n_clean_repeats", "n_vetted", "n_alias_ok", "t2", "depth_ref", "dur_ref_d"]
OWN = ["false_event", "repeat_detected", "call", "truth_in_alias_set", "n_vetted", "n_alias_ok", "t2"]
TOL = 1e-9


def norm(s, kind):
    if kind == "bool":
        return s.astype(str).str.strip().str.lower().map({"true": True, "false": False})
    return s


def same(x, y, kind):
    if kind == "num":
        x, y = pd.to_numeric(x), pd.to_numeric(y)
        return ((x - y).abs() <= TOL * np.maximum(1, y.abs())) | (x.isna() & y.isna())
    x, y = norm(x, kind), norm(y, kind)
    return (x == y) | (x.isna() & y.isna())


out = []
for arm, pf, jf, key in (("A", "../test/arm_a_rows.csv", "ind_a_rows.csv", ["tic", "toi", "arm"]),
                         ("B", "../test/arm_b_rows.csv", "ind_b_rows.csv", ["tic", "draw"])):
    p, j = pd.read_csv(H / pf), pd.read_csv(H / jf)
    assert len(p) == len(j) and not p.duplicated(key).any() and not j.duplicated(key).any()
    p = j[key].merge(p, on=key, how="left", validate="1:1")   # primary rows put in the independent (defs) order
    assert p.status.notna().all() and np.allclose(p.t_cat, j.t_cat, rtol=0, atol=1e-12)
    print(f"\n### Arm {arm}: {len(p)} trials, matched one to one on {key}; primary status: {p.status.value_counts().to_dict()}\n")
    print("| column | agree | disagree | max abs difference |\n|---|---|---|---|")
    cols = [(c, "bool") for c in BOOL] + [(c, "str") for c in STR] + [(c, "num") for c in NUM]
    if arm == "B":
        cols += [(c + "_own", "bool" if c in BOOL else "str" if c in STR else "num") for c in OWN]
        cols += [("null_event", "bool"), ("n_null_events", "num")]
    for c, kind in cols:
        ok = same(p[c], j[c], kind)
        md = ""
        if kind == "num":
            d = (pd.to_numeric(p[c]) - pd.to_numeric(j[c])).abs()
            md = f"{d.max():.3g}" if d.notna().any() else ""
        print(f"| {c} | {int(ok.sum())} | {int((~ok).sum())} | {md} |")
        for i in np.flatnonzero(~ok.values):
            out.append(dict(arm=p.arm[i], tic=p.tic[i], id=p[key[1]][i], column=c, primary=p[c][i], independent=j[c][i]))
print("\n### Disagreeing trials\n")
if out:
    print(pd.DataFrame(out).to_markdown(index=False))
else:
    print("None.")
