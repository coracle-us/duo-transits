#!/usr/bin/env python3
"""Score the blind independent re-derivation against the pre-registration's metrics.

Purpose
    Collect the per-star files of ind_arm_a.py and ind_arm_b.py, write the trial tables and print
    every metric of prereg_mono_calibration.yaml v1.1 (false_exclusion_rate overall and by
    reference-SNR bin, repeat_recovery_rate, wrong_period_rate, median fraction of ln P excluded),
    then the pass/fail verdict under pass_criteria and outcomes.

Inputs
    work/arm_a/*.json, work/arm_b/*.json

Outputs
    arm_a_trials.csv, arm_b_trials.csv, score_tables.md (the tables quoted in REPORT.md)

Usage
    python ind_score.py
"""
import glob
import json

import numpy as np
import pandas as pd
from scipy.stats import beta

import ind_common as ic

SNR_EDGES = [-np.inf, 6, 8, 12, 20, np.inf]
SNR_LABELS = ["<6", "6-8", "8-12", "12-20", ">20"]
FIRST = ["tic", "arm", "status", "ref_snr", "truth_P", "t_cat", "truth_excluded", "n_passing_events", "call",
         "truth_in_alias_set", "repeat_detected", "n_clean_repeats", "frac_lnP_excluded"]
B_EXTRA = ["draw", "inj_P", "inj_depth_ppm", "inj_b", "inj_duration_h"]
lines = []


def say(s=""):
    print(s)
    lines.append(s)


def load(arm):
    rows = []
    for p in sorted(glob.glob(str(ic.WORK / arm / "*.json"))):
        rows += json.loads(open(p).read())
    df = pd.DataFrame(rows)
    return df.drop(columns=[c for c in ("traceback",) if c in df])


def ub(k, n):
    """One-sided 95% Clopper-Pearson upper bound."""
    return float(beta.ppf(0.95, k + 1, n - k)) if n > 0 and k < n else (1.0 if n else np.nan)


def rate(k, n):
    return f"{k}/{n}" + (f" = {100*k/n:.1f}% (95% UB {100*ub(k, n):.1f}%)" if n else " (n/a)")


def snr_bin(df):
    b = pd.cut(df["ref_snr"], SNR_EDGES, labels=SNR_LABELS, right=False).astype(object)
    return b.where(df["ref_snr"].notna(), "no ref box")


def metrics(df, name, by_extra=()):
    """Returns dict of (k, n) for the pass criteria; prints the tables."""
    df = df.copy()
    df["snr_bin"] = snr_bin(df)
    te = df["truth_excluded"].astype(bool)
    out = {}
    say(f"### {name}: {len(df)} scored trials on {df['tic'].nunique()} stars")
    say()
    say(f"- false_exclusion_rate overall: {rate(int(te.sum()), len(df))}")
    out["fe"] = (int(te.sum()), len(df))
    out["fe_bins"] = {}
    say()
    say("| reference SNR | false exclusions | median frac ln P excluded | median sectors |")
    say("|---|---|---|---|")
    for lab in SNR_LABELS + ["no ref box"]:
        g = df[df["snr_bin"] == lab]
        if lab == "no ref box" and not len(g):
            continue
        k = int(g["truth_excluded"].astype(bool).sum())
        out["fe_bins"][lab] = (k, len(g))
        med = f"{g['frac_lnP_excluded'].median():.3f}" if len(g) else "-"
        ns = f"{g['n_sectors'].median():.0f}" if len(g) else "-"
        say(f"| {lab} | {rate(k, len(g))} | {med} | {ns} |")
    say()
    rr = df[df["n_clean_repeats"] >= 1]
    out["rr"] = (int(rr["repeat_detected"].astype(bool).sum()), len(rr))
    say(f"- repeat_recovery_rate (trials with >= 1 other true transit at a clean epoch): {rate(*out['rr'])}")
    cl = df[df["call"].isin(["resolved", "narrowed"])]
    wrong = ~cl["truth_in_alias_set"].astype(bool)
    out["wp"] = (int(wrong.sum()), len(cl))
    say(f"- wrong_period_rate (called resolved or narrowed, truth not in the set): {rate(*out['wp'])}")
    for c in ("resolved", "narrowed"):
        g = cl[cl["call"] == c]
        say(f"  - {c}: {rate(int((~g['truth_in_alias_set'].astype(bool)).sum()), len(g))}")
    if len(cl):
        ev_true = cl["alias_event_is_true_repeat"].astype(bool)
        say(f"  - of the {int(wrong.sum())} wrong: {int((wrong & ~ev_true).sum())} had an alias event that is not a "
            f"transit of the planet, {int((wrong & ev_true).sum())} had a true repeat whose alias was excluded")
    say(f"- calls: " + ", ".join(f"{k} {v}" for k, v in df["call"].value_counts().items()))
    say(f"- median fraction of ln P (1-1000 d) excluded: {df['frac_lnP_excluded'].median():.3f} "
        f"(quartiles {df['frac_lnP_excluded'].quantile(.25):.3f}-{df['frac_lnP_excluded'].quantile(.75):.3f})")
    df["nsec_bin"] = pd.cut(df["n_sectors"], [0, 2, 5, 10, 20, 1000], labels=["1-2", "3-5", "6-10", "11-20", ">20"])
    say("- median frac ln P excluded by number of sectors: " + ", ".join(
        f"{k}: {g['frac_lnP_excluded'].median():.3f} (n={len(g)})" for k, g in df.groupby("nsec_bin", observed=True)))
    say(f"- reference measured by the pipeline's box (snr >= 5, depth within 0.33-3x): "
        f"{int((df['ref_used'] == 'box').sum())}/{len(df)}; catalogue depth and duration used for the rest")
    for col in by_extra:
        say()
        say(f"| {col} | false exclusions | repeat recovery | wrong period | median frac ln P |")
        say("|---|---|---|---|---|")
        for key, g in df.groupby(col, observed=True):
            r2 = g[g["n_clean_repeats"] >= 1]
            c2 = g[g["call"].isin(["resolved", "narrowed"])]
            say(f"| {key} | {rate(int(g['truth_excluded'].astype(bool).sum()), len(g))} | "
                f"{rate(int(r2['repeat_detected'].astype(bool).sum()), len(r2))} | "
                f"{rate(int((~c2['truth_in_alias_set'].astype(bool)).sum()), len(c2))} | "
                f"{g['frac_lnP_excluded'].median():.3f} |")
    say()
    return out


def write_csv(df, path, extra=()):
    cols = FIRST + list(extra)
    for c in cols:
        if c not in df:
            df[c] = np.nan
    df[cols + [c for c in df.columns if c not in cols]].to_csv(path, index=False)


def main():
    a, b = load("arm_a"), load("arm_b")
    write_csv(a, ic.HERE / "arm_a_trials.csv")
    write_csv(b, ic.HERE / "arm_b_trials.csv", B_EXTRA)
    res = {}

    say("## Arm A")
    say()
    n_star = a["tic"].nunique()
    bad = a[a["status"] != "ok"]
    ok = a[a["status"] == "ok"].copy()
    say(f"Stars in the sample file with an output: {n_star}. Finished: {ok['tic'].nunique()}. "
        f"Skipped or failed: {bad['tic'].nunique()}.")
    for st, g in bad.groupby("status"):
        say(f"- {st}: {len(g)} star(s): TIC {', '.join(str(t) for t in g['tic'])}")
    say()
    ok["multi_toi"] = np.where(ok["n_toi_on_star"] > 1, "n_toi_on_star > 1", "n_toi_on_star == 1")
    for arm in ("A1", "A2"):
        g = ok[ok["arm"] == arm]
        res[arm] = metrics(g, f"{arm} ({'hidden' if arm == 'A1' else 'visible'}), all stars", by_extra=("multi_toi", "depth_bin"))
        for key, gg in g.groupby("multi_toi"):
            res[f"{arm} {key}"] = metrics(gg, f"{arm}, {key}")
    ex = ok[(ok["arm"] == "A2") & ok["truth_excluded"].astype(bool)]
    say(f"### A2 stars where the true period is excluded (for human review): {len(ex)}")
    say()
    if len(ex):
        say("| TIC | TOI | n TOI on star | P true (d) | ref SNR | ref depth (ppm) | ref dur (h) | TOI depth (ppm) | TOI dur (h) "
            "| clean true epochs | excluding | excluding epochs (BTJD: box depth) |")
        say("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in ex.itertuples():
            say(f"| {r.tic} | {r.toi} | {r.n_toi_on_star} | {r.truth_P:.5f} | {r.ref_snr:.1f} | {r.depth_ref_ppm:.0f} | "
                f"{r.dur_ref_h:.2f} | {r.toi_depth_ppm:.0f} | {r.toi_duration_h:.2f} | {r.n_truth_epochs_clean} | "
                f"{r.n_truth_epochs_excluding} | {r.truth_excluding_epochs} |")
        say()

    say("## Arm B")
    say()
    badb = b[b["status"] != "ok"]
    okb = b[b["status"] == "ok"].copy()
    say(f"Hosts with an output: {b['tic'].nunique()}. Trials scored: {len(okb)} on {okb['tic'].nunique()} hosts. "
        f"Trials or hosts not scored: {len(badb)}.")
    for st, g in badb.groupby("status"):
        say(f"- {st}: {len(g)}: TIC {', '.join(str(t) for t in g['tic'])}")
    say(f"- hosts on solar mass and radius: {okb.loc[okb['solar_values'].astype(bool), 'tic'].nunique()}; "
        f"hosts cooler than 3500 K (3500 K coefficients): {okb.loc[okb['cool_star_3500'].astype(bool), 'tic'].nunique()}; "
        f"hosts without a TIC Teff (5772 K used): {okb.loc[okb['teff_missing'].astype(bool), 'tic'].nunique()}")
    say()
    okb["inj_depth_bin"] = pd.cut(okb["inj_depth_ppm"], [500, 2000, 5000, 20000],
                                  labels=["500-2000", "2000-5000", "5000-20000"], include_lowest=True)
    okb["ref_cadence"] = okb["ref_author"].astype(str) + " " + okb["ref_exptime_s"].round(0).astype("Int64").astype(str) + " s"
    okb["inj_P_bin"] = pd.cut(okb["inj_P"], [20, 50, 100, 300, 1000], labels=["20-50", "50-100", "100-300", "300-1000"],
                              include_lowest=True)
    res["B"] = metrics(okb, "B (injections)", by_extra=("inj_depth_bin", "ref_cadence", "tmag_bin", "inj_P_bin"))
    exb = okb[okb["truth_excluded"].astype(bool)]
    say(f"### B trials where the true period is excluded: {len(exb)}")
    say()
    if len(exb):
        say("| TIC | draw | P (d) | depth (ppm) | b | T14 (h) | ref SNR | ref used | ref depth (ppm) | ref dur (h) | excluding epochs |")
        say("|---|---|---|---|---|---|---|---|---|---|---|")
        for r in exb.itertuples():
            say(f"| {r.tic} | {r.draw} | {r.inj_P:.3f} | {r.inj_depth_ppm:.0f} | {r.inj_b:.2f} | {r.inj_duration_h:.2f} | "
                f"{r.ref_snr:.1f} | {r.ref_used} | {r.depth_ref_ppm:.0f} | {r.dur_ref_h:.2f} | {r.truth_excluding_epochs} |")
        say()

    # ---------------- verdict
    say("## Verdict under pass_criteria and outcomes")
    say()
    n_a, n_b = ok["tic"].nunique(), len(okb)
    incomplete = n_a < 60 or n_b < 400
    say(f"- completeness: {n_a} arm-A stars (need >= 60), {n_b} arm-B trials (need >= 400): "
        f"{'INCOMPLETE' if incomplete else 'complete'}")
    fails = []
    for arm in ("A2", "B"):
        k, n = res[arm]["fe"]
        okk = n > 0 and k / n <= 0.05
        say(f"- false exclusion, {arm}, overall: {rate(k, n)} -> {'met' if okk else 'MISSED'}")
        if not okk:
            fails.append(f"false exclusion {arm} overall {k}/{n}")
        for lab in SNR_LABELS[1:]:
            k, n = res[arm]["fe_bins"][lab]
            okk = n == 0 or k / n <= 0.05
            say(f"  - SNR {lab}: {rate(k, n)} -> {'met' if okk else 'MISSED'}" + (" (empty bin)" if n == 0 else ""))
            if not okk:
                fails.append(f"false exclusion {arm} SNR {lab} {k}/{n}")
    for arm in ("A1", "A2", "B"):
        k, n = res[arm]["wp"]
        okk = n == 0 or k / n <= 0.03
        say(f"- wrong period, {arm}: {rate(k, n)} -> {'met' if okk else 'MISSED'}")
        if not okk:
            fails.append(f"wrong period {arm} {k}/{n}")
    say()
    verdict = "INCOMPLETE" if incomplete else ("FAIL" if fails else "PASS")
    say(f"**Outcome: {verdict}**" + (": " + "; ".join(fails) if fails else ""))
    (ic.HERE / "score_tables.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
