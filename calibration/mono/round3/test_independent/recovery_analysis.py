#!/usr/bin/env python3
"""Step 5 (analysis made AFTER the outcome was known): where repeat recovery is lost.

Uses only the independent rows and the defs. Denominator as the criterion: trials with at least one
true repeat in clean data (A2 one-TOI stars; arm B).

Input   ind_a_rows.csv, ind_b_rows.csv, ../test/arm_a_defs.csv, ../test/arm_b_defs.csv
Output  recovery_analysis.md
Usage   python recovery_analysis.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def tf(s):
    return s.astype(str).str.strip().str.lower() == "true"


def load():
    a = pd.read_csv(HERE / "ind_a_rows.csv")
    da = pd.read_csv(HERE.parent / "test" / "arm_a_defs.csv")
    a["depth_true"], a["dur_true"] = da["depth_cat"], da["dur_cat_d"]
    b = pd.read_csv(HERE / "ind_b_rows.csv")
    db = pd.read_csv(HERE.parent / "test" / "arm_b_defs.csv")
    b["depth_true"], b["dur_true"], b["inj_b"] = db["inj_depth"], db["inj_dur_d"], db["inj_b"]
    b["cat_depth_ratio"] = db["depth_cat"] / db["inj_depth"]
    b["cat_dur_ratio"] = db["dur_cat_d"] / db["inj_dur_d"]
    out = {}
    for name, df in (("A2 one-TOI", a[(a.arm == "A2") & ~tf(a.multi)]), ("B", b)):
        df = df[df.n_clean_repeats >= 1].copy()
        for c in ("repeat_detected", "repeat_detected_r4", "ref_measured", "true_ev_r4"):
            df[c] = tf(df[c])
        df["why"] = np.select(
            [df.repeat_detected, df.repeat_detected_r4, df.n_true_ev > 0],
            ["detected", "lost to T5 only (passes R4, sec_box_snr < 9 or undefined)", "fails R4"],
            "no event at the truth among the top 3 of its sector")
        df["snr_bin"] = np.where(~df.ref_measured, "unmeasured",
                                 pd.cut(df.ref_snr, [-np.inf, 8, 12, 20, np.inf],
                                        labels=["5-8", "8-12", "12-20", "over 20"], right=False).astype(str))
        df["exp_snr"] = df.depth_true / df.sig
        out[name] = df
    return out


def split(df, col, order=None):
    g = df.groupby(col, observed=True)
    t = pd.DataFrame(dict(n=g.size(), detected=g.repeat_detected.sum(), r4_only=g.repeat_detected_r4.sum()))
    t["recovery"] = (100 * t.detected / t.n).round(1)
    t["without_T5"] = (100 * t.r4_only / t.n).round(1)
    t["lost_to_T5"] = t.r4_only - t.detected
    if order:
        t = t.reindex([o for o in order if o in t.index])
    return t.reset_index()


def md_table(t):
    lines = ["| " + " | ".join(map(str, t.columns)) + " |", "|" + "---|" * len(t.columns)]
    for _, r in t.iterrows():
        lines.append("| " + " | ".join(str(x) for x in r) + " |")
    return lines


def main():
    D = load()
    md = ["# Repeat recovery: where it is lost", "",
          "**Analysis made after the outcome was known. Not part of the pre-registered test.**",
          "From the independent rows only (`ind_a_rows.csv`, `ind_b_rows.csv`) and the defs.",
          "Denominator: trials with at least one true repeat in clean data.", ""]
    for name, df in D.items():
        n = len(df)
        k3, k2 = int(df.repeat_detected.sum()), int(df.repeat_detected_r4.sum())
        md += [f"## {name}", "",
               f"Round 3 (R4 and `sec_box_snr >= 9`): {k3} of {n} ({100 * k3 / n:.1f}%). "
               f"Without the round-3 veto (R4 only): {k2} of {n} ({100 * k2 / n:.1f}%).", ""]
        w = df.why.value_counts().rename_axis("fate of the true repeat").reset_index(name="trials")
        w["share"] = (100 * w.trials / n).round(1)
        md += md_table(w) + [""]
        r4 = df[df.why == "fails R4"]
        md += ["R4 conditions failed by the strongest event at the truth (a trial can fail several):", ""]
        md += md_table(r4.true_ev_fail.value_counts().rename_axis("failed").reset_index(name="trials")) + [""]
        t5 = df[df.why.str.startswith("lost to T5")]
        if len(t5):
            q = t5.true_ev_sbs
            md += [f"Lost to T5 only: best `sec_box_snr` of an R4-passing event at the truth - undefined in "
                   f"{int(q.isna().sum())}, otherwise min {q.min():.1f}, median {q.median():.1f}, max {q.max():.1f}; "
                   f"event SNR median {t5.true_ev_snr.median():.1f} (min {t5.true_ev_snr.min():.1f}, max {t5.true_ev_snr.max():.1f}).", ""]
        md += ["By reference SNR:", ""] + md_table(split(df, "snr_bin", ["unmeasured", "5-8", "8-12", "12-20", "over 20"])) + [""]
        df["depth_bin"] = pd.cut(df.depth_true * 1e6, [0, 1000, 2000, 5000, 1e9],
                                 labels=["under 1000 ppm", "1000-2000", "2000-5000", "over 5000"]).astype(str)
        md += ["By true depth (catalogued for arm A, injected for arm B):", ""] + md_table(
            split(df, "depth_bin", ["under 1000 ppm", "1000-2000", "2000-5000", "over 5000"])) + [""]
        df["dur_bin"] = pd.cut(df.dur_true * 24, [0, 3, 6, 12, 1e3], labels=["under 3 h", "3-6 h", "6-12 h", "over 12 h"]).astype(str)
        md += ["By true duration:", ""] + md_table(split(df, "dur_bin", ["under 3 h", "3-6 h", "6-12 h", "over 12 h"])) + [""]
        df["exp_bin"] = pd.cut(df.exp_snr, [0, 5, 8, 12, 20, 1e9], labels=["under 5", "5-8", "8-12", "12-20", "over 20"]).astype(str)
        md += ["By expected SNR of one repeat (true depth / duration-scale noise of the merged data):", ""] + md_table(
            split(df, "exp_bin", ["under 5", "5-8", "8-12", "12-20", "over 20"])) + [""]
        df["rep_bin"] = np.where(df.n_clean_repeats >= 3, "3 or more", df.n_clean_repeats.astype(str))
        md += ["By number of clean repeats:", ""] + md_table(split(df, "rep_bin", ["1", "2", "3 or more"])) + [""]
        md += ["By product and cadence (min) of the sectors holding the clean repeats:", ""] + md_table(
            split(df, "repeat_products").sort_values("n", ascending=False)) + [""]
        df["ref_prod"] = df.ref_author.fillna("none") + ":" + df.ref_cad_min.round(0).astype("Int64").astype(str)
        md += ["By product and cadence of the reference sector:", ""] + md_table(
            split(df, "ref_prod").sort_values("n", ascending=False)) + [""]
        if name == "B":
            df["b_bin"] = pd.cut(df.inj_b, [0, 0.3, 0.6, 0.9], labels=["0-0.3", "0.3-0.6", "0.6-0.9"]).astype(str)
            md += ["By impact parameter:", ""] + md_table(split(df, "b_bin", ["0-0.3", "0.3-0.6", "0.6-0.9"])) + [""]
            df["cd_bin"] = pd.cut(df.cat_depth_ratio, [0, 0.9, 1.1, 9], labels=["under 0.9", "0.9-1.1", "over 1.1"]).astype(str)
            md += ["By catalogued / injected depth (the 10% catalogue error):", ""] + md_table(
                split(df, "cd_bin", ["under 0.9", "0.9-1.1", "over 1.1"])) + [""]
            df["cu_bin"] = pd.cut(df.cat_dur_ratio, [0, 0.9, 1.1, 9], labels=["under 0.9", "0.9-1.1", "over 1.1"]).astype(str)
            md += ["By catalogued / injected duration:", ""] + md_table(
                split(df, "cu_bin", ["under 0.9", "0.9-1.1", "over 1.1"])) + [""]
        df.to_csv(HERE / f"recovery_trials_{name.split()[0]}.csv", index=False)
    # what the veto bought, on the judged false-event samples
    a = pd.read_csv(HERE / "ind_a_rows.csv")
    b = pd.read_csv(HERE / "ind_b_rows.csv")
    a1 = a[(a.arm == "A1") & ~tf(a.multi)]
    md += ["## What the round-3 veto removed on the false-event side (same rows)", "",
           f"- A1 one-TOI stars with a passing event: R4 only {int((a1.n_r4 > 0).sum())} of {len(a1)}, "
           f"round 3 {int((a1.n_vetted > 0).sum())} of {len(a1)}.",
           f"- B trials with a false event, own events removed: R4 only {int(tf(b.false_event_r4_own).sum())} of {len(b)}, "
           f"round 3 {int(tf(b.false_event_own).sum())} of {len(b)}.",
           f"- B trials with a false event, all events: R4 only {int(tf(b.false_event_r4).sum())} of {len(b)}, "
           f"round 3 {int(tf(b.false_event).sum())} of {len(b)}.", ""]
    (HERE / "recovery_analysis.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
