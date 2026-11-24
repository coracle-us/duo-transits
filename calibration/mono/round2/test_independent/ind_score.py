#!/usr/bin/env python3
"""Independent scorer for the round-2 test, written from calibration/mono/prereg_round2.yaml (1.1)
(metrics, pass_criteria, outcomes) and the R5 table of FROZEN_RULE.md only.

Input   rows CSVs with the columns of ind_test.py (or any file with the same column names)
Output  a markdown file with the criteria table and the outcome
Usage   python ind_score.py [--a ind_a_rows.csv] [--b ind_b_rows.csv] [--out ind_score.md]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
BINS = [("6-8", 6, 8), ("8-12", 8, 12), ("12-20", 12, 20), (">20", 20, np.inf)]   # pass-criteria bins
MIN_BIN = 60
# FROZEN_RULE.md R5: p = (k + 0.5) / (n + 1) from the tuning set
R5 = [("unmeasured", 0, 507), ("5-8", 0, 211), ("8-12", 0, 202), ("12-20", 0, 192), (">20", 1, 511)]
R5_P = {name: (k + 0.5) / (n + 1) for name, k, n in R5}


def tf(s):
    """Boolean column that may have been read as text or object; missing stays missing."""
    if s.dtype == bool:
        return s
    m = s.astype(str).str.strip().str.lower().map({"true": True, "false": False, "1": True, "0": False,
                                                    "1.0": True, "0.0": False})
    return m


def load(path):
    df = pd.read_csv(path)
    for c in ["multi", "ref_measured", "truth_excluded", "false_event", "repeat_detected", "truth_in_alias_set",
              "false_event_own", "repeat_detected_own", "truth_in_alias_set_own", "null_event"]:
        if c in df:
            df[c] = tf(df[c])
    df["scored"] = df["truth_excluded"].notna()
    return df


def ub95(k, n):
    """One-sided 95% Clopper-Pearson upper bound."""
    if n == 0:
        return np.nan
    return 1.0 if k == n else float(stats.beta.ppf(0.95, k + 1, n - k))


def r5_bin(row):
    if not bool(row["ref_measured"]) or not np.isfinite(row["ref_snr"]) or row["ref_snr"] < 5:
        return "unmeasured"
    s = row["ref_snr"]
    return "5-8" if s < 8 else "8-12" if s < 12 else "12-20" if s < 20 else ">20"


class Table:
    def __init__(self):
        self.rows = []

    def rate(self, group, name, k, n, limit, kind, judged=True, note=""):
        """kind 'max': pass if k/n <= limit; 'min': pass if k/n >= limit."""
        k, n = int(k), int(n)
        r = k / n if n else np.nan
        if not judged:
            v = "not judged"
        elif n == 0:
            v = "no trials"
        else:
            v = "met" if (r <= limit if kind == "max" else r >= limit) else "MISSED"
        self.rows.append(dict(group=group, criterion=name, k=k, n=n, value=r, ub95=ub95(k, n) if kind == "max" else np.nan,
                              limit=("<= " if kind == "max" else ">= ") + f"{limit:g}", verdict=v, note=note))
        return v

    def value(self, group, name, value, n, limit_txt, ok, judged=True, note=""):
        v = "not judged" if not judged else ("no trials" if n == 0 else ("met" if ok else "MISSED"))
        self.rows.append(dict(group=group, criterion=name, k=np.nan, n=int(n), value=value, ub95=np.nan,
                              limit=limit_txt, verdict=v, note=note))
        return v

    def md(self, groups=None):
        out = ["| Criterion | Sample | k | n | value | 95% upper | limit | verdict |", "|---|---|---|---|---|---|---|---|"]
        for r in self.rows:
            if groups and r["group"] not in groups:
                continue
            k = "" if not np.isfinite(r["k"]) else f"{int(r['k'])}"
            val = "" if not np.isfinite(r["value"]) else f"{r['value']:.4f}"
            ub = "" if not np.isfinite(r["ub95"]) else f"{r['ub95']:.4f}"
            out.append(f"| {r['group']} | {r['criterion']}{(' (' + r['note'] + ')') if r['note'] else ''} | {k} | {r['n']} "
                       f"| {val} | {ub} | {r['limit']} | {r['verdict']} |")
        return "\n".join(out)


def score(A, B):
    T = Table()
    a = A[A["scored"]]
    b = B[B["scored"]]
    a1, a2 = a[(a.arm == "A1") & ~a.multi.astype(bool)], a[(a.arm == "A2") & ~a.multi.astype(bool)]
    m1, m2 = a[(a.arm == "A1") & a.multi.astype(bool)], a[(a.arm == "A2") & a.multi.astype(bool)]
    judged = {}

    # ---- false exclusion (A2 one-TOI and B; overall and per SNR bin >= 60 trials)
    fe = []
    for name, d in (("A2 one-TOI", a2), ("B all trials", b)):
        fe.append(T.rate("false_exclusion", f"{name}, overall", d.truth_excluded.sum(), len(d), 0.05, "max"))
        for bn, lo, hi in BINS:
            s = d[(d.ref_snr >= lo) & (d.ref_snr < hi)]
            j = len(s) >= MIN_BIN
            v = T.rate("false_exclusion", f"{name}, ref SNR {bn}", s.truth_excluded.sum(), len(s), 0.05, "max", judged=j,
                       note="" if j else "under 60 trials")
            if j:
                fe.append(v)
        s = d[~(d.ref_snr >= 6)]
        T.rate("false_exclusion", f"{name}, ref SNR under 6 or none", s.truth_excluded.sum(), len(s), 0.05, "max", judged=False)
    judged["false_exclusion_rate"] = all(v == "met" for v in fe)

    # ---- ln P excluded (median at ref SNR >= 8; A1 one-TOI and B)
    lp = []
    for name, d in (("A1 one-TOI", a1), ("B all trials", b)):
        s = d[d.ref_snr >= 8]
        med = float(s.frac_lnP_excluded.median()) if len(s) else np.nan
        lp.append(T.value("log_period_excluded", f"{name}, ref SNR >= 8, median", med, len(s), ">= 0.5", med >= 0.5))
    judged["log_period_excluded"] = all(v == "met" for v in lp)

    # ---- false event
    ev = [T.rate("false_event", "A1 one-TOI stars", a1.false_event.sum(), len(a1), 0.05, "max"),
          T.rate("false_event", "B, own events removed", b.false_event_own.sum(), len(b), 0.05, "max")]
    T.rate("false_event", "B, all events", b.false_event.sum(), len(b), 0.05, "max", judged=False, note="reported beside")
    if "null_event" in b:
        T.rate("false_event", "B null runs with a passing event", b.null_event.astype(bool).sum(), len(b), 0.05, "max",
               judged=False, note="reported, not judged")
    judged["false_event_rate"] = all(v == "met" for v in ev)

    # ---- wrong period
    def wrong(d, sfx=""):
        s = d[d["call" + sfx].isin(["resolved", "narrowed"])]
        return int((~s["truth_in_alias_set" + sfx].astype(bool)).sum()), len(s)
    wp = [T.rate("wrong_period", "A2 one-TOI", *wrong(a2), 0.03, "max"),
          T.rate("wrong_period", "B, own events removed", *wrong(b, "_own"), 0.03, "max")]
    T.rate("wrong_period", "B, all events", *wrong(b), 0.03, "max", judged=False, note="reported beside")
    judged["wrong_period_rate"] = all(v == "met" for v in wp)

    # ---- repeat recovery
    def rec(d, sfx=""):
        s = d[d.n_clean_repeats >= 1]
        return int(s["repeat_detected" + sfx].astype(bool).sum()), len(s)
    rr = [T.rate("repeat_recovery", "A2 one-TOI", *rec(a2), 0.60, "min"),
          T.rate("repeat_recovery", "B, own events removed", *rec(b, "_own"), 0.60, "min")]
    T.rate("repeat_recovery", "B, all events", *rec(b), 0.60, "min", judged=False, note="reported beside")
    judged["repeat_recovery_rate"] = all(v == "met" for v in rr)

    # ---- probability calibration (R5)
    both = pd.concat([a2.assign(_s="A2 one-TOI"), b.assign(_s="B")])
    both["r5"] = both.apply(r5_bin, axis=1)
    cal = []
    for name, _, _ in R5:
        s = both[both.r5 == name]
        cal.append(dict(bin=name, p=R5_P[name], n=len(s), predicted=len(s) * R5_P[name],
                        observed=int(s.truth_excluded.sum()),
                        n_A2=int((s._s == "A2 one-TOI").sum()), n_B=int((s._s == "B").sum())))
    cal = pd.DataFrame(cal)
    mu, obs = float(cal.predicted.sum()), int(cal.observed.sum())
    lo, hi = int(stats.poisson.ppf(0.025, mu)), int(stats.poisson.ppf(0.975, mu))
    v = T.value("probability_calibration", f"A2 one-TOI + B: observed {obs}, predicted {mu:.2f}, 95% Poisson interval {lo}-{hi}",
                obs, len(both), f"in {lo}-{hi}", lo <= obs <= hi)
    judged["probability_calibration"] = v == "met"

    # ---- multi-TOI stars, reported only
    T.rate("multi-TOI (not judged)", "false exclusion, A2", m2.truth_excluded.sum(), len(m2), 0.05, "max", judged=False)
    T.rate("multi-TOI (not judged)", "false event, A1", m1.false_event.sum(), len(m1), 0.05, "max", judged=False)
    T.rate("multi-TOI (not judged)", "wrong period, A2", *wrong(m2), 0.03, "max", judged=False)
    T.rate("multi-TOI (not judged)", "repeat recovery, A2", *rec(m2), 0.60, "min", judged=False)
    s = m1[m1.ref_snr >= 8]
    T.value("multi-TOI (not judged)", "ln P excluded, A1, ref SNR >= 8, median",
            float(s.frac_lnP_excluded.median()) if len(s) else np.nan, len(s), ">= 0.5", True, judged=False)

    # ---- outcome
    n_one = len(set(a1.tic) & set(a2.tic))
    n_b = len(b)
    core_ok = judged["false_exclusion_rate"] and judged["log_period_excluded"]
    events_ok = judged["false_event_rate"] and judged["wrong_period_rate"] and judged["repeat_recovery_rate"]
    if n_one < 80 or n_b < 400:
        outcome = "incomplete"
    elif all(judged.values()):
        outcome = "pass"
    elif not core_ok:
        outcome = "fail"
    elif judged["probability_calibration"] and not events_ok:
        outcome = "partial"
    else:
        outcome = "not defined by the pre-registration (false exclusion and ln P met, probability calibration missed)"
    return dict(table=T, judged=judged, outcome=outcome, n_one=n_one, n_b=n_b, cal=cal,
                n_a_rows=len(A), n_a_scored=len(a), n_b_rows=len(B), n_multi=len(set(m1.tic) | set(m2.tic)),
                failed_a=A[~A["scored"]], failed_b=B[~B["scored"]])


def report(R, title, src):
    L = [f"# {title}", "", f"Rows: {src}", "",
         f"Scored: {R['n_one']} one-TOI arm-A stars (both A1 and A2 scored), {R['n_multi']} multi-TOI stars, "
         f"{R['n_b']} of {R['n_b_rows']} arm-B trials; {R['n_a_scored']} of {R['n_a_rows']} arm-A rows scored.", "",
         f"## Outcome: **{R['outcome']}**", "",
         "| Criterion | met |", "|---|---|"]
    L += [f"| {k} | {'yes' if v else 'NO'} |" for k, v in R["judged"].items()]
    L += ["", "## Criteria table", "", R["table"].md(), "",
          "'95% upper' is the one-sided 95% Clopper-Pearson bound (reported, not part of any criterion).", "",
          "## Probability calibration by R5 bin", "",
          "| R5 bin | p | trials (A2 one-TOI + B) | predicted | observed |", "|---|---|---|---|---|"]
    for r in R["cal"].itertuples():
        L.append(f"| {r.bin} | {r.p:.4f} | {r.n} ({r.n_A2} + {r.n_B}) | {r.predicted:.3f} | {r.observed} |")
    L.append(f"| total | | {int(R['cal'].n.sum())} | {R['cal'].predicted.sum():.3f} | {int(R['cal'].observed.sum())} |")
    for nm, d in (("A", R["failed_a"]), ("B", R["failed_b"])):
        if len(d):
            L += ["", f"## Arm {nm} rows not scored ({len(d)})", ""]
            L += [f"- {r.arm} TIC {r.tic}: {str(r.notes)[:160]}" for r in d.itertuples()]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", default=str(HERE / "ind_a_rows.csv"))
    ap.add_argument("--b", default=str(HERE / "ind_b_rows.csv"))
    ap.add_argument("--out", default=str(HERE / "ind_score.md"))
    ap.add_argument("--title", default="Round 2 test: independent score (ind_test.py rows, ind_score.py)")
    a = ap.parse_args()
    R = score(load(a.a), load(a.b))
    Path(a.out).write_text(report(R, a.title, f"`{Path(a.a).name}`, `{Path(a.b).name}`"))
    print(report(R, a.title, f"{a.a}, {a.b}"))


if __name__ == "__main__":
    main()
