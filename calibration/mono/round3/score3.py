"""
Score the round-3 test against the pass criteria of prereg_round3.yaml (1.0). Written before the
freeze and before any test-set outcome existed.

Inputs  : <dir>/arm_a_rows.csv, arm_b_rows.csv (harness3.py)
Outputs : <dir>/score_rates.csv, score_summary.md
Usage   : python calibration/mono/round3/score3.py [dir, default calibration/mono/round3/test]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
MIN_A_STARS, MIN_B_TRIALS, MIN_A_CALLED = 60, 600, 40
CALLED = ["resolved", "narrowed"]


def upper95(k, n):
    return np.nan if n == 0 else 1.0 if k == n else float(beta.ppf(0.95, k + 1, n - k))


def line(name, k, n, limit, at_most=True, judged=True):
    r = k / n if n else np.nan
    ok = (r <= limit) if at_most else (r >= limit)
    verdict = "reported, not judged" if not judged else "no trials" if not n else "pass" if ok else "FAIL"
    return dict(criterion=name, k=k, n=n, rate=r, upper95=upper95(k, n), limit=limit, verdict=verdict)


def main():
    d = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "test"
    a = pd.read_csv(d / "arm_a_rows.csv")
    b = pd.read_csv(d / "arm_b_rows.csv")
    a, b = a[a.status == "ok"].copy(), b[b.status == "ok"].copy()
    one = ~a.multi.astype(bool)
    a1, a2, a2m = a[(a.arm == "A1") & one], a[(a.arm == "A2") & one], a[(a.arm == "A2") & ~one]
    rows = []
    for name, t, suf, judged in (("A2 one-TOI", a2, "", True), ("B, own events removed", b, "_own", True),
                                 ("A2 multi-TOI", a2m, "", False), ("B, all events", b, "", False)):
        c = t[t["call" + suf].isin(CALLED)]
        enough = judged and (not name.startswith("A2") or len(c) >= MIN_A_CALLED)
        wrong = ~c["truth_in_alias_set" + suf].astype(bool)
        rows.append(line(f"wrong period, {name}", int(wrong.sum()), len(c), 0.03, judged=enough))
        ct = c[c["event_true" + suf].astype(bool)]
        rows.append(line(f"wrong period (true event), {name}",
                         int((~ct["truth_in_alias_set" + suf].astype(bool)).sum()), len(ct), 0.03, judged=enough))
        r = t[t.n_clean_repeats > 0]
        rows.append(line(f"repeat recovery, {name}", int(r["repeat_detected" + suf].astype(bool).sum()), len(r),
                         0.60, at_most=False, judged=judged))
        m3 = float(ct["frac_alias_excluded" + suf].median()) if len(ct) else np.nan
        m2 = float(ct["frac_alias_excluded_r2" + suf].median()) if len(ct) else np.nan
        ratio = m3 / m2 if m2 else np.nan
        rows.append(dict(criterion=f"alias power, {name}: median excluded {m3:.4f}, round-2 rule {m2:.4f}",
                         k=np.nan, n=len(ct), rate=ratio, upper95=np.nan, limit=0.8,
                         verdict="reported, not judged" if not judged else "no trials" if not len(ct)
                         else "pass" if ratio >= 0.8 else "FAIL"))
    rows.append(line("false event, A1 one-TOI", int((a1.n_vetted > 0).sum()), len(a1), 0.05))
    rows.append(line("false event, B, own events removed", int(b.false_event_own.astype(bool).sum()), len(b), 0.05))
    rows.append(line("false event, B, all events", int(b.false_event.astype(bool).sum()), len(b), 0.05, judged=False))
    rows.append(line("false exclusion, A2 one-TOI", int(a2.truth_excluded.sum()), len(a2), 0.05))
    rows.append(line("false exclusion, B", int(b.truth_excluded.sum()), len(b), 0.05))
    r = pd.DataFrame(rows)
    r.to_csv(d / "score_rates.csv", index=False)

    def met(prefix):
        g = r[r.criterion.str.startswith(prefix) & r.verdict.isin(["pass", "FAIL", "no trials"])]
        return bool((g.verdict == "pass").all())
    wrong, wtrue, power = met("wrong period, "), met("wrong period (true event)"), met("alias power")
    rec, fev, fex = met("repeat recovery"), met("false event"), met("false exclusion")
    n_a, n_b = a2.tic.nunique(), len(b)
    if n_a < MIN_A_STARS or n_b < MIN_B_TRIALS:
        outcome = "incomplete"
    elif wrong and wtrue and power and rec and fev and fex:
        outcome = "pass"
    elif not (wtrue and power):
        outcome = "fail"
    elif rec and fex:
        outcome = "calls after vetting"
    else:
        outcome = ("not defined by the pre-registration: wrong period (true event) and alias power met, "
                   "repeat recovery or false exclusion missed")
    flag = "" if fex else " **Exclusion regression: the false-exclusion criterion was missed.**"
    out = ["# Round 3 test: score (pre-registration 1.0)", "",
           f"Arm A: {n_a} one-TOI stars scored ({a[a.arm == 'A2'].tic.nunique()} stars in all). "
           f"Arm B: {n_b} trials on {b.tic.nunique()} hosts.", "", f"**Outcome: {outcome}**{flag}", "",
           r.to_markdown(index=False, floatfmt=".3f")]
    (d / "score_summary.md").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
