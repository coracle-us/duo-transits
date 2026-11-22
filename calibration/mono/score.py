"""
Score the single-transit calibration against the pass criteria of prereg_mono_calibration.yaml
(version 1.1). Written before any trial outcome was looked at.

Inputs  : calibration/mono/results/arm_a_trials.csv, arm_b_trials.csv, arm_a_stars.csv, arm_b_stars.csv
Outputs : calibration/mono/results/score_rates.csv, score_summary.md
Usage   : python calibration/mono/score.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

RES = Path(__file__).resolve().parent / "results"
SNR_EDGES = [0, 6, 8, 12, 20, np.inf]
SNR_LABELS = ["<6", "6-8", "8-12", "12-20", ">20"]
MAX_FALSE_EXCL, MAX_WRONG = 0.05, 0.03
MIN_A_STARS, MIN_B_TRIALS = 60, 400


def upper95(k, n):
    """One-sided 95% Clopper-Pearson upper bound on a rate."""
    return 1.0 if n == 0 or k == n else float(beta.ppf(0.95, k + 1, n - k))


def rate(df, col, label, **extra):
    n, k = len(df), int(df[col].sum()) if len(df) else 0
    return dict(metric=label, n=n, k=k, rate=k / n if n else np.nan, upper95=upper95(k, n) if n else np.nan, **extra)


def score(trials):
    """Rows of rates for one data frame of trials (one arm or a subset of it)."""
    t = trials[trials.status == "ok"].copy()
    t["snr_bin"] = pd.cut(t.ref_snr, SNR_EDGES, labels=SNR_LABELS, right=False).astype(str)
    rows = [rate(t, "truth_excluded_single", "false_exclusion", snr_bin="all")]
    for b in SNR_LABELS:
        rows.append(rate(t[t.snr_bin == b], "truth_excluded_single", "false_exclusion", snr_bin=b))
    tested = t[t.truth_n_clean_epochs > 0]
    rows.append(rate(tested, "truth_excluded_single", "false_exclusion_given_a_clean_epoch", snr_bin="all"))
    rep = t[t.n_clean_repeats > 0]
    rows.append(rate(rep, "repeat_detected", "repeat_recovery", snr_bin="all"))
    for b in SNR_LABELS:
        rows.append(rate(rep[rep.snr_bin == b], "repeat_detected", "repeat_recovery", snr_bin=b))
    called = t[t.call.isin(["resolved", "narrowed"])].copy()
    called["wrong"] = ~called.truth_in_alias_set.astype(bool)
    rows.append(rate(called, "wrong", "wrong_period", snr_bin="all"))
    rows.append(rate(called[called.call == "resolved"], "wrong", "wrong_period_resolved_only", snr_bin="all"))
    for b in ["all"] + SNR_LABELS:
        g = t if b == "all" else t[t.snr_bin == b]
        rows.append(dict(metric="median_frac_lnP_excluded", n=len(g), k=np.nan,
                         rate=g.frac_lnP_excluded.median() if len(g) else np.nan, upper95=np.nan, snr_bin=b))
    return pd.DataFrame(rows)


def verdict(r, arm, judge_false_excl):
    """Pass or fail lines for one arm from its rate table."""
    out, ok = [], True
    if judge_false_excl:
        fe = r[r.metric == "false_exclusion"]
        for x in fe.itertuples():
            if x.snr_bin == "<6" or x.n == 0:
                continue
            good = x.rate <= MAX_FALSE_EXCL
            ok &= good
            out.append(f"- {arm} false exclusion, SNR {x.snr_bin}: {x.k} of {x.n} ({x.rate:.1%}, "
                       f"95% upper bound {x.upper95:.1%}) - {'pass' if good else 'FAIL'}")
    w = r[r.metric == "wrong_period"].iloc[0]
    if w.n:
        good = w.rate <= MAX_WRONG
        ok &= good
        out.append(f"- {arm} wrong period: {w.k} of {w.n} ({w.rate:.1%}, 95% upper bound {w.upper95:.1%}) "
                   f"- {'pass' if good else 'FAIL'}")
    else:
        out.append(f"- {arm} wrong period: no trial was called resolved or narrowed")
    return out, ok


def main():
    a = pd.read_csv(RES / "arm_a_trials.csv")
    b = pd.read_csv(RES / "arm_b_trials.csv")
    sa, sb = pd.read_csv(RES / "arm_a_stars.csv"), pd.read_csv(RES / "arm_b_stars.csv")
    sets = {"A1": a[a.arm == "A1"], "A2": a[a.arm == "A2"], "B": b,
            "A1 single-TOI stars": a[(a.arm == "A1") & (a.n_toi_on_star == 1)],
            "A2 single-TOI stars": a[(a.arm == "A2") & (a.n_toi_on_star == 1)]}
    tables = {k: score(v).assign(arm=k) for k, v in sets.items()}
    pd.concat(tables.values()).to_csv(RES / "score_rates.csv", index=False)

    n_a = a[a.status == "ok"].tic.nunique()
    n_b = int((b.status == "ok").sum())
    lines = ["# Single-transit calibration: score (pre-registration version 1.1)", "",
             f"Arm A: {n_a} stars scored of {len(sa)} drawn "
             f"(skipped: {sa.skipped.fillna('').replace('', np.nan).dropna().value_counts().to_dict()}; "
             f"errors: {int(sa.get('error', pd.Series(dtype=str)).notna().sum())}).",
             f"Arm B: {n_b} trials scored on {b[b.status == 'ok'].tic.nunique()} of {len(sb)} hosts "
             f"(errors: {int(sb.get('error', pd.Series(dtype=str)).notna().sum())}).", ""]
    complete = n_a >= MIN_A_STARS and n_b >= MIN_B_TRIALS
    oks = []
    for arm, judge in (("A1", False), ("A2", True), ("B", True)):
        ls, ok = verdict(tables[arm], arm, judge)
        lines += ls
        oks.append(ok)
    outcome = "incomplete" if not complete else "pass" if all(oks) else "fail"
    lines += ["", f"**Outcome: {outcome}**", "", "## All rates", "",
              pd.concat(tables.values())[["arm", "metric", "snr_bin", "n", "k", "rate", "upper95"]]
              .to_markdown(index=False, floatfmt=".3f")]
    (RES / "score_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:len(lines) - 4]))


if __name__ == "__main__":
    main()
