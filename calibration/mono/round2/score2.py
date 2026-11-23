"""
Score the round-2 test against the pass criteria of prereg_round2.yaml (1.1). Written before the
freeze and before any test-set star was loaded.

Inputs  : <dir>/arm_a_rows.csv, <dir>/arm_b_rows.csv (harness.py rows)
Outputs : <dir>/score_rates.csv, <dir>/score_summary.md
Usage   : python calibration/mono/round2/score2.py [dir, default calibration/mono/round2/test]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, poisson

HERE = Path(__file__).resolve().parent
SNR_EDGES, SNR_LABELS = [0, 6, 8, 12, 20, np.inf], ["<6", "6-8", "8-12", "12-20", ">20"]
MIN_BIN, MIN_A_STARS, MIN_B_TRIALS = 60, 80, 400
# FROZEN_RULE.md, R5: (upper edge of reference SNR, p); unmeasured references use the first row
P_TABLE = [(5, 0.0010), (8, 0.0024), (12, 0.0025), (20, 0.0026), (np.inf, 0.0029)]


def p_false_exclusion(ref_snr, ref_measured):
    if not ref_measured or not np.isfinite(ref_snr):
        return P_TABLE[0][1]
    return next(p for edge, p in P_TABLE if ref_snr < edge)


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
    for t in (a, b):
        t["snr_bin"] = pd.cut(t.ref_snr, SNR_EDGES, labels=SNR_LABELS, right=False).astype(str)
    one = ~a.multi.astype(bool)
    a1, a2, a2m = a[(a.arm == "A1") & one], a[(a.arm == "A2") & one], a[(a.arm == "A2") & ~one]
    rows = []
    # false exclusion
    for name, t in (("A2 one-TOI", a2), ("B", b)):
        rows.append(line(f"false exclusion, {name}, all", int(t.truth_excluded.sum()), len(t), 0.05))
        for s in SNR_LABELS[1:]:
            g = t[t.snr_bin == s]
            rows.append(line(f"false exclusion, {name}, SNR {s}", int(g.truth_excluded.sum()), len(g), 0.05,
                             judged=len(g) >= MIN_BIN))
    rows.append(line("false exclusion, A2 multi-TOI", int(a2m.truth_excluded.sum()), len(a2m), 0.05, judged=False))
    # ln P excluded
    for name, t in (("A1 one-TOI", a1), ("B", b)):
        g = t[t.ref_snr >= 8]
        med = float(g.frac_lnP_excluded.median()) if len(g) else np.nan
        rows.append(dict(criterion=f"median ln P excluded, {name}, SNR >= 8", k=np.nan, n=len(g), rate=med,
                         upper95=np.nan, limit=0.5, verdict="pass" if med >= 0.5 else "FAIL"))
    # events
    rows.append(line("false event, A1 one-TOI", int((a1.n_vetted > 0).sum()), len(a1), 0.05))
    rows.append(line("false event, B, own events removed", int(b.false_event_own.astype(bool).sum()), len(b), 0.05))
    rows.append(line("false event, B, all events", int(b.false_event.astype(bool).sum()), len(b), 0.05, judged=False))
    rows.append(line("passing event in B null runs", int(b.null_event.astype(bool).sum()), len(b), 0.05, judged=False))
    for name, t, suf, judged in (("A2 one-TOI", a2, "", True), ("B, own events removed", b, "_own", True),
                                 ("B, all events", b, "", False)):
        c = t[t["call" + suf].isin(["resolved", "narrowed"])]
        rows.append(line(f"wrong period, {name}", int((~c["truth_in_alias_set" + suf].astype(bool)).sum()),
                         len(c), 0.03, judged=judged))
        r = t[t.n_clean_repeats > 0]
        rows.append(line(f"repeat recovery, {name}", int(r["repeat_detected" + suf].astype(bool).sum()), len(r),
                         0.60, at_most=False, judged=judged))
    # probability calibration
    both = pd.concat([a2, b])
    mu = float(sum(p_false_exclusion(s, m) for s, m in zip(both.ref_snr, both.ref_measured.astype(bool))))
    k = int(both.truth_excluded.sum())
    lo, hi = int(poisson.ppf(0.025, mu)), int(poisson.ppf(0.975, mu))
    rows.append(dict(criterion=f"probability calibration: observed {k}, predicted {mu:.2f}, 95% interval {lo}-{hi}",
                     k=k, n=len(both), rate=mu, upper95=np.nan, limit=np.nan,
                     verdict="pass" if lo <= k <= hi else "FAIL"))
    r = pd.DataFrame(rows)
    r.to_csv(d / "score_rates.csv", index=False)

    def met(prefixes):
        g = r[r.criterion.str.startswith(tuple(prefixes)) & r.verdict.isin(["pass", "FAIL", "no trials"])]
        return bool((g.verdict == "pass").all())
    excl, lnp, prob = met(["false exclusion"]), met(["median ln P"]), met(["probability"])
    events = met(["false event", "wrong period", "repeat recovery"])
    n_a, n_b = a2.tic.nunique(), len(b)
    if n_a < MIN_A_STARS or n_b < MIN_B_TRIALS:
        outcome = "incomplete"
    elif not (excl and lnp):
        outcome = "fail"
    elif prob and events:
        outcome = "pass"
    elif prob:
        outcome = "partial"
    else:
        outcome = "not defined by the pre-registration: exclusion criteria met, probability calibration missed"
    out = ["# Round 2 test: score (pre-registration 1.1)", "",
           f"Arm A: {n_a} one-TOI stars scored ({a[a.arm == 'A2'].tic.nunique()} stars in all). "
           f"Arm B: {n_b} trials on {b.tic.nunique()} hosts.", "", f"**Outcome: {outcome}**", "",
           r.to_markdown(index=False, floatfmt=".3f")]
    (d / "score_summary.md").write_text("\n".join(out) + "\n")
    print("\n".join(out[:5]))


if __name__ == "__main__":
    main()
