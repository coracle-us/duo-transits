#!/usr/bin/env python3
"""Step 4: compare the independent rows with the primary rows, column by column.

Input   ind_a_rows.csv, ind_b_rows.csv; ../test/arm_a_rows.csv, ../test/arm_b_rows.csv
Output  compare_columns.csv, compare_disagreements.csv, and a table on stdout
Usage   python compare.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CALL = ["n_vetted", "false_event", "repeat_detected", "t2", "call", "n_alias_ok", "n_alias",
        "truth_in_alias_set", "event_true", "frac_alias_excluded", "frac_alias_excluded_r2"]
BASE = ["ref_snr", "ref_measured", "depth_ref", "dur_ref_d", "n_events", "n_clean_repeats", "truth_excluded"]
# (independent column, primary column)
PAIRS_A = [(c, c) for c in BASE + CALL] + [("n_r4", "r2_n_vetted"), ("false_event_r4", "r2_false_event"),
                                           ("repeat_detected_r4", "r2_repeat_detected"),
                                           ("ref_author", "ref_author"), ("ref_cad_min", "ref_cadence_min")]
PAIRS_B = PAIRS_A + [(c + "_own", c + "_own") for c in CALL] + [
    ("n_null_events", "n_null_events"), ("n_r4_own", "r2_n_vetted_own"),
    ("false_event_r4_own", "r2_false_event_own"), ("repeat_detected_r4_own", "r2_repeat_detected_own")]
RTOL = 1e-6


def same(x, y):
    """Element-wise agreement: both empty, equal strings/booleans, or numbers within RTOL."""
    out = np.zeros(len(x), bool)
    for i, (u, v) in enumerate(zip(x, y)):
        un, vn = pd.isna(u), pd.isna(v)
        if un or vn:
            out[i] = un and vn
            continue
        su, sv = str(u).strip().lower(), str(v).strip().lower()
        if su in ("true", "false") or sv in ("true", "false"):
            out[i] = su == sv
            continue
        try:
            fu, fv = float(u), float(v)
            out[i] = abs(fu - fv) <= RTOL * max(1.0, abs(fu), abs(fv))
        except ValueError:
            out[i] = su == sv
    return out


def main():
    cols, dis = [], []
    for arm, pairs, key in (("a", PAIRS_A, ["tic", "arm"]), ("b", PAIRS_B, ["tic", "draw"])):
        ind = pd.read_csv(HERE / f"ind_{arm}_rows.csv")
        pri = pd.read_csv(HERE.parent / "test" / f"arm_{arm}_rows.csv")
        assert len(ind) == len(pri)
        assert (ind[key].to_numpy() == pri[key].to_numpy()).all(), "row order differs"
        print(f"arm {arm}: {len(ind)} rows; primary status: {pri['status'].value_counts().to_dict()}")
        for ci, cp in pairs:
            ok = same(ind[ci].to_numpy(object), pri[cp].to_numpy(object))
            cols.append(dict(file=arm, ind_col=ci, primary_col=cp, agree=int(ok.sum()), disagree=int((~ok).sum())))
            for i in np.flatnonzero(~ok):
                dis.append(dict(file=arm, **{k: ind.loc[i, k] for k in key}, column=ci,
                                independent=ind.loc[i, ci], primary=pri.loc[i, cp]))
    c = pd.DataFrame(cols)
    c.to_csv(HERE / "compare_columns.csv", index=False)
    d = pd.DataFrame(dis)
    d.to_csv(HERE / "compare_disagreements.csv", index=False)
    print(c.to_string(index=False))
    print(f"\n{len(d)} disagreeing cells")
    if len(d):
        print(d.to_string(index=False))


if __name__ == "__main__":
    main()
