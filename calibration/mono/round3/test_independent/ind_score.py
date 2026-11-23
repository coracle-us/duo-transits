#!/usr/bin/env python3
"""Independent scorer for the round-3 test, written from prereg_round3.yaml only
(prereg_round2.yaml for the definitions it inherits). Does not read the primary scorer or its output.

Input   ind_a_rows.csv, ind_b_rows.csv (ind_test3.py)
Output  ind_score.md, ind_score.json
Usage   python ind_score.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CALLED = ("resolved", "narrowed")


def tf(s):
    return s.astype(str).str.strip().str.lower() == "true"


def rate(k, n):
    return f"{k} of {n} ({100 * k / n:.1f}%)" if n else f"{k} of {n}"


def main():
    a = pd.read_csv(HERE / "ind_a_rows.csv")
    b = pd.read_csv(HERE / "ind_b_rows.csv")
    for df in (a, b):
        df["scored"] = df["truth_excluded"].notna()
        df["multi"] = tf(df["multi"])
    n_fail_a, n_fail_b = int((~a.scored).sum()), int((~b.scored).sum())
    a, b = a[a.scored].copy(), b[b.scored].copy()
    for df, sufs in ((a, [""]), (b, ["", "_own"])):
        df["truth_excluded"] = tf(df["truth_excluded"])
        for s in sufs:
            for c in ("false_event", "repeat_detected", "truth_in_alias_set", "event_true"):
                df[c + s] = tf(df[c + s])
    a1, a2 = a[(a.arm == "A1") & ~a.multi], a[(a.arm == "A2") & ~a.multi]
    a1m, a2m = a[(a.arm == "A1") & a.multi], a[(a.arm == "A2") & a.multi]
    n_a_stars = len(set(a1.tic) & set(a2.tic))

    rows, verdict = [], {}

    def add(crit, arm, k, n, ok, judged=True, note=""):
        v = ("met" if ok else "MISSED") if judged else "not judged"
        rows.append((crit, arm, k, n, (k / n if n else np.nan), v, note))
        if judged:
            verdict[crit] = verdict.get(crit, True) and bool(ok)

    def report(crit, arm, k, n, note=""):
        rows.append((crit, arm, k, n, (k / n if n else np.nan), "reported", note))

    # ---- wrong period
    def wrong(df, s, true_only):
        m = df["call" + s].isin(CALLED)
        if true_only:
            m &= df["event_true" + s]
        return int((~df.loc[m, "truth_in_alias_set" + s]).sum()), int(m.sum())

    for crit, to in (("wrong_period_rate", False), ("wrong_period_rate_true_event", True)):
        k, n = wrong(a2, "", to)
        n_called = wrong(a2, "", False)[1]
        judged = n >= 40
        add(crit, "A2 one-TOI", k, n, n > 0 and k / n <= 0.03, judged,
            "" if judged else f"fewer than 40 trials in the denominator (called trials: {n_called})")
        k, n = wrong(b, "_own", to)
        add(crit, "B own removed", k, n, n > 0 and k / n <= 0.03)
        k, n = wrong(b, "", to)
        report(crit, "B all events", k, n)
        k, n = wrong(a2m, "", to)
        report(crit, "A2 multi-TOI", k, n)

    # ---- false events
    k, n = int((a1.n_vetted > 0).sum()), len(a1)
    add("false_event_rate", "A1 one-TOI", k, n, k / n <= 0.05, note="stars with any passing event")
    report("false_event_rate", "A1 one-TOI, passing event not at a true transit", int(a1.false_event.sum()), n)
    k, n = int(b.false_event_own.sum()), len(b)
    add("false_event_rate", "B own removed", k, n, k / n <= 0.05)
    report("false_event_rate", "B all events", int(b.false_event.sum()), len(b))
    report("false_event_rate", "A1 multi-TOI", int((a1m.n_vetted > 0).sum()), len(a1m))
    report("false_event_rate", "A2 one-TOI (not a criterion)", int(a2.false_event.sum()), len(a2))
    if "n_null_vetted" in b:
        report("null runs with a passing event", "B", int((b.n_null_vetted > 0).sum()), len(b))

    # ---- repeat recovery
    m = a2.n_clean_repeats >= 1
    k, n = int(a2.loc[m, "repeat_detected"].sum()), int(m.sum())
    add("repeat_recovery_rate", "A2 one-TOI", k, n, n > 0 and k / n >= 0.60)
    m = b.n_clean_repeats >= 1
    k, n = int(b.loc[m, "repeat_detected_own"].sum()), int(m.sum())
    add("repeat_recovery_rate", "B own removed", k, n, n > 0 and k / n >= 0.60)
    report("repeat_recovery_rate", "B all events", int(b.loc[m, "repeat_detected"].sum()), n)
    m = a2m.n_clean_repeats >= 1
    report("repeat_recovery_rate", "A2 multi-TOI", int(a2m.loc[m, "repeat_detected"].sum()), int(m.sum()))

    # ---- alias power
    power = []

    def pw(df, s, arm, judged):
        m = df["call" + s].isin(CALLED) & df["event_true" + s]
        r3 = float(df.loc[m, "frac_alias_excluded" + s].median()) if m.any() else np.nan
        r2 = float(df.loc[m, "frac_alias_excluded_r2" + s].median()) if m.any() else np.nan
        ok = bool(m.any() and r3 >= 0.8 * r2)
        power.append((arm, int(m.sum()), r3, r2, (r3 / r2 if r2 else np.nan),
                      ("met" if ok else "MISSED") if judged else "reported"))
        if judged:
            verdict["alias_power"] = verdict.get("alias_power", True) and ok

    pw(a2, "", "A2 one-TOI", True)
    pw(b, "_own", "B own removed", True)
    pw(a[a.arm == "A2"], "", "A2 all stars", False)
    pw(b, "", "B all events", False)

    # ---- false exclusion (replication of round 2's rule)
    k, n = int(a2.truth_excluded.sum()), len(a2)
    add("false_exclusion_rate", "A2 one-TOI", k, n, k / n <= 0.05)
    k, n = int(b.truth_excluded.sum()), len(b)
    add("false_exclusion_rate", "B", k, n, k / n <= 0.05)
    report("false_exclusion_rate", "A2 multi-TOI", int(a2m.truth_excluded.sum()), len(a2m))

    # ---- outcome
    v = verdict
    core4 = all(v[c] for c in ("wrong_period_rate_true_event", "alias_power", "repeat_recovery_rate",
                               "false_exclusion_rate"))
    lines = []
    if n_a_stars < 60 or len(b) < 600:
        outcome = "incomplete"
    elif all(v.values()):
        outcome = "pass"
    elif not v["wrong_period_rate_true_event"] or not v["alias_power"]:
        outcome = "fail"
    elif core4 and (not v["false_event_rate"] or not v["wrong_period_rate"]):
        outcome = "calls_after_vetting"
    else:
        outcome = "NOT DEFINED by the pre-registration"
        missed = [c for c, ok in v.items() if not ok]
        lines.append(
            "The outcomes section of prereg_round3.yaml does not define this combination: missed = "
            + ", ".join(missed) + "; wrong_period_rate_true_event and alias_power are met, so it is not "
            "`fail`; not every criterion is met, so it is not `pass`; `calls_after_vetting` requires "
            "repeat_recovery_rate and false_exclusion_rate to be met, so it is not that either.")
    if not v["false_exclusion_rate"]:
        lines.insert(0, "EXCLUSION REGRESSION: false_exclusion_rate missed (to be reported at the top, "
                        "whatever else happened).")

    md = ["# Round 3 test: independent score", "",
          "From `ind_a_rows.csv` and `ind_b_rows.csv` (`ind_test3.py`), scored by `ind_score.py`, written from",
          "`prereg_round3.yaml` only. Made before the primary rows and score were read.", "",
          f"Scored: {n_a_stars} one-TOI arm-A stars (A1 and A2 both run; {len(a1m)} multi-TOI stars reported, "
          f"not judged), {len(b)} arm-B trials on {b.tic.nunique()} hosts. "
          f"Failed runs: {n_fail_a} arm-A rows, {n_fail_b} arm-B rows.", "",
          f"## Outcome: **{outcome}**", ""]
    md += lines + ([""] if lines else [])
    md += ["| Criterion | Limit | Verdict |", "|---|---|---|"]
    lim = dict(wrong_period_rate="<= 3% (A2 needs >= 40 called)", wrong_period_rate_true_event="<= 3% (same minimum)",
               false_event_rate="<= 5%", repeat_recovery_rate=">= 60%",
               alias_power="round-3 median >= 0.8 x round-2 median", false_exclusion_rate="<= 5%")
    for c in lim:
        md.append(f"| {c} | {lim[c]} | {'met' if v[c] else 'MISSED'} |")
    md += ["", "## Per arm", "", "| Criterion | Arm | k | n | rate | verdict | note |", "|---|---|---|---|---|---|---|"]
    for crit, arm, k, n, r, vv, note in rows:
        md.append(f"| {crit} | {arm} | {k} | {n} | {'' if not n else f'{100 * r:.1f}%'} | {vv} | {note} |")
    md += ["", "## Alias power (called trials whose alias event is a true transit)", "",
           "| Arm | n | median fraction excluded, round 3 | round 2 rule, same trials | ratio | verdict |",
           "|---|---|---|---|---|---|"]
    for arm, n, r3, r2, q, vv in power:
        md.append(f"| {arm} | {n} | {r3:.4f} | {r2:.4f} | {q:.3f} | {vv} |")
    md += ["", "## Calls", "", "| Arm | resolved | narrowed | all excluded | no repeat |", "|---|---|---|---|---|"]
    for name, df, s in (("A1 one-TOI", a1, ""), ("A2 one-TOI", a2, ""), ("B own removed", b, "_own"),
                        ("B all events", b, "")):
        c = df["call" + s].value_counts()
        md.append(f"| {name} | " + " | ".join(str(int(c.get(x, 0))) for x in
                                              ("resolved", "narrowed", "all excluded", "no repeat")) + " |")
    md += ["", "## Trials behind every miss", ""]

    def listing(title, df, cols):
        md.append(f"**{title}** ({len(df)})")
        md.append("")
        if len(df):
            md.append("| " + " | ".join(cols) + " |")
            md.append("|" + "---|" * len(cols))
            for _, r in df.iterrows():
                md.append("| " + " | ".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + " |")
            md.append("")

    base = ["tic", "toi", "ref_snr", "truth_P", "t2", "t2_snr", "call", "n_alias_ok", "n_alias", "event_true"]
    listing("A2 one-TOI wrong calls", a2[a2.call.isin(CALLED) & ~a2.truth_in_alias_set], base)
    bb = b.rename(columns={c + "_own": c + "@" for c in ("t2", "t2_snr", "call", "n_alias_ok", "n_alias", "event_true")})
    bcols = ["tic", "draw", "ref_snr", "truth_P", "t2@", "t2_snr@", "call@", "n_alias_ok@", "n_alias@", "event_true@"]
    listing("B (own removed) wrong calls", bb[b.call_own.isin(CALLED) & ~b.truth_in_alias_set_own], bcols)
    listing("A1 one-TOI stars with a passing event", a1[a1.n_vetted > 0], ["tic", "toi", "ref_snr", "truth_P", "n_vetted", "t2", "t2_snr", "call"])
    listing("B (own removed) trials with a false event", bb[b.false_event_own], bcols)
    listing("False exclusions, A2 one-TOI", a2[a2.truth_excluded], ["tic", "toi", "ref_snr", "truth_P"])
    listing("False exclusions, B", b[b.truth_excluded], ["tic", "draw", "ref_snr", "truth_P"])
    listing("False exclusions, A2 multi-TOI (not judged)", a2m[a2m.truth_excluded], ["tic", "toi", "ref_snr", "truth_P"])
    (HERE / "ind_score.md").write_text("\n".join(md) + "\n")
    (HERE / "ind_score.json").write_text(json.dumps(dict(
        outcome=outcome, verdict=v,
        rows=[dict(criterion=c, arm=ar, k=k, n=n, verdict=vv) for c, ar, k, n, _, vv, _ in rows],
        alias_power=[dict(arm=ar, n=n, r3=r3, r2=r2, verdict=vv) for ar, n, r3, r2, _, vv in power]), indent=1))
    print("\n".join(md))


if __name__ == "__main__":
    main()
