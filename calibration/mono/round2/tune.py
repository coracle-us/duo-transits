"""
Round 2 tuning, step 2: score rule variants on the stored tuning-set summaries (dump.py).

Metrics are those of prereg_round2.yaml. Differences from the frozen harness, which only matter for
ranking variants: the fraction of ln P excluded uses every 4th period of the grid, and alias tables
are looked up on the dur/8 time grid instead of at exact epochs. The true period is tested at its
exact epochs, as in pre-registration 1.1.

Inputs  : calibration/mono/work/round2_tuning/arm_a/*.pkl, arm_b/*.pkl
Outputs : calibration/mono/round2/tuning/<name>.csv (one row per variant)
Usage   : python calibration/mono/round2/tune.py exclusion|vetting|final [--workers 12]
"""
import argparse
import itertools
import multiprocessing as mp
import pickle
import sys
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import monocal as mc  # noqa: E402
import rule2  # noqa: E402
from rule2 import Rule  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
WORK = ROOT / "calibration/mono/work/round2_tuning"
OUT = HERE / "tuning"
SNR_EDGES, SNR_LABELS = [0, 6, 8, 12, 20, np.inf], ["<6", "6-8", "8-12", "12-20", ">20"]
GRID_THIN = 4
TRIALS = None


def load():
    global TRIALS
    if TRIALS is None:
        TRIALS = []
        for arm in ("arm_a", "arm_b"):
            for p in sorted((WORK / arm).glob("*.pkl")):
                TRIALS += pickle.loads(p.read_bytes())["trials"]
    return TRIALS


def upper95(k, n):
    return np.nan if n == 0 else 1.0 if k == n else float(beta.ppf(0.95, k + 1, n - k))


def grid_masks(s, rule):
    multi = s["n_toi_on_star"] > 1
    ex = rule2.timeline_excludes(s["g_depth"], s["g_frac"], s["g_clean"], s["g_idx"], s["n_grid"], s["sig"],
                                 s["depth_ref"], s["derr_ref"], rule, s["ref_measured"], multi)
    err = s["sig"] / np.sqrt(np.clip(s["g_frac"], 1e-3, None))
    su = np.zeros(s["n_grid"], bool)
    su[s["g_idx"]] = s["g_clean"] & (np.abs(s["g_depth"] - s["depth_ref"]) < 3 * np.hypot(err, s["derr_ref"]))
    return ex, su


def alias_lookup(s, t1, t2, ex, su):
    """Status of each alias P_n = (t2 - t1)/n >= 1 d from the time grid: (excluded, supported)."""
    dT = t2 - t1
    n = np.arange(1, int(dT // mc.P_FLOOR) + 1)[:3000]
    exc, sup = np.zeros(n.size, bool), np.zeros(n.size, bool)
    for i, nn in enumerate(n):
        P = dT / nn
        k = np.arange(int(np.floor((s["tmin"] - t1) / P)) - 1, int(np.ceil((s["tmax"] - t1) / P)) + 2)
        k = k[(k != 0) & (k != nn)]
        g = np.rint((t1 + k * P - s["c0"]) / s["step"]).astype(np.int64)
        g = g[(g >= 0) & (g < ex.size)]
        exc[i], sup[i] = ex[g].any(), su[g].any()
    return n, dT / n, exc, sup


def score_trial(s, rule, want_grid, want_call):
    multi = s["n_toi_on_star"] > 1
    r = dict(arm=s["arm"], tic=s["tic"], multi=multi, ref_snr=s["ref_snr"], ref_measured=s["ref_measured"],
             ref_cadence_min=s["ref_cadence_min"], n_clean_repeats=s.get("n_clean_repeats", 0))
    te = rule2.epoch_excludes(s["tr_depth"], s["tr_frac"], s["tr_edge"], s["sig"], s["depth_ref"],
                              s["derr_ref"], rule, s["ref_measured"], multi) if s["tr_depth"].size else np.zeros(0, bool)
    r["truth_excluded"] = bool(te.any())
    ev = [e for e in s["events"] if "error" not in e]
    ok = [e for e in ev if rule2.vet(e, rule)]
    r["n_vetted"] = len(ok)
    r["n_false_vetted"] = sum(not e.get("at_truth", False) for e in ok)
    r["repeat_detected"] = any(e.get("at_truth", False) for e in ok)
    if "null_events" in s:
        r["null_event"] = any(rule2.vet(e, rule) for e in s["null_events"])
    best = rule2.best_event(ev, rule)
    r["best_at_truth"] = bool(best and best.get("at_truth", False))
    ex = su = None
    if want_grid and s["ref_snr"] >= 8 and s["arm"] in ("A1", "B", "Bx"):
        ex, su = grid_masks(s, rule)
        per = mc.period_grid(s["tmax"] - s["tmin"], s["dur_ref_d"])[::GRID_THIN]
        r["frac_lnP"] = float(mc.excluded_periods(s["t_cat"], per, s["c0"], s["step"], ex, s["tmin"], s["tmax"]).mean())
    if want_call:
        r["call"] = "no repeat"
        if best:
            if ex is None:
                ex, su = grid_masks(s, rule)
            t1, t2 = sorted([s["t_cat"], best["tc"]])
            n, P, exc, sup = alias_lookup(s, t1, t2, ex, su)
            pmin = s["pmin"]
            live = ~exc & (sup | ~(np.isfinite(pmin) & (P < pmin)))
            r["call"] = "resolved" if live.sum() == 1 else "narrowed" if live.any() else "all excluded"
            dT = t2 - t1
            nt = int(round(dT / s["truth_P"]))
            on = 1 <= nt <= n.size and abs(dT - nt * s["truth_P"]) < s["dur_ref_d"] / 2
            r["truth_in_set"] = bool(on and not exc[nt - 1])
    return r


def _frac(df, col):
    n, k = len(df), int(df[col].sum()) if len(df) else 0
    return k, n


def metrics(rule, want_grid=True, want_call=False):
    t = pd.DataFrame([score_trial(s, rule, want_grid, want_call) for s in load()])
    t["snr_bin"] = pd.cut(t.ref_snr, SNR_EDGES, labels=SNR_LABELS, right=False).astype(str)
    m = dict(asdict(rule))
    sets = {"A2one": t[(t.arm == "A2") & ~t.multi], "A2multi": t[(t.arm == "A2") & t.multi],
            "B": t[t.arm == "B"], "Ball": t[t.arm.isin(["B", "Bx"])], "A1one": t[(t.arm == "A1") & ~t.multi]}
    for name in ("A2one", "A2multi", "B", "Ball"):
        k, n = _frac(sets[name], "truth_excluded")
        m[f"fe_{name}"] = f"{k}/{n}"
        m[f"fe_{name}_rate"] = k / n if n else np.nan
    for b in SNR_LABELS[1:]:
        for name in ("A2one", "Ball"):
            k, n = _frac(sets[name][sets[name].snr_bin == b], "truth_excluded")
            m[f"fe_{name}_{b}"] = f"{k}/{n}"
    if want_grid:
        m["lnP_A1"] = t[t.arm == "A1"].frac_lnP.median()
        m["lnP_B"] = t[t.arm == "B"].frac_lnP.median()
        m["lnP_Ball"] = t[t.arm.isin(["B", "Bx"])].frac_lnP.median()
    a1 = sets["A1one"]
    m["false_event_A1one"] = f"{int((a1.n_vetted > 0).sum())}/{len(a1)}"
    m["false_event_A1one_rate"] = float((a1.n_vetted > 0).mean())
    for name in ("B", "Ball"):
        d = sets[name].dropna(subset=["null_event"])
        m[f"null_event_{name}"] = f"{int(d.null_event.sum())}/{len(d)}"
        m[f"null_event_{name}_rate"] = float(d.null_event.astype(bool).mean())
        m[f"null_event_{name}_hosts"] = f"{int(d.groupby('tic').null_event.any().sum())}/{d.tic.nunique()}"
    for name in ("A2one", "B", "Ball"):
        d = sets[name][sets[name].n_clean_repeats > 0]
        k, n = _frac(d, "repeat_detected")
        m[f"recovery_{name}"] = f"{k}/{n}"
        m[f"recovery_{name}_rate"] = k / n if n else np.nan
    if want_call:
        for name in ("A2one", "B", "Ball"):
            d = sets[name][sets[name].call.isin(["resolved", "narrowed"])]
            k = int((~d.truth_in_set.astype(bool)).sum())
            m[f"wrong_{name}"] = f"{k}/{len(d)}"
            m[f"wrong_{name}_rate"] = k / len(d) if len(d) else np.nan
    return m


def _run(args):
    return metrics(*args)


def run(variants, name, workers, want_grid, want_call):
    OUT.mkdir(exist_ok=True)
    load()
    with mp.get_context("fork").Pool(workers) as pool:
        rows = pool.map(_run, [(v, want_grid, want_call) for v in variants], chunksize=1)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / f"{name}.csv", index=False)
    return df


def exclusion_variants():
    out = []
    for fl, mt, wk, w in itertools.product([np.inf, 3.0, 2.0], [-np.inf, 0.0, 1.0, 2.0], [True, False],
                                           [0.0, 0.5, 1.0]):
        out.append(Rule(floor_sig=fl, min_thr_sig=mt, weak_ref_excludes=wk, w_dur=w))
    return out


def vetting_variants(base):
    out = []
    for sn, dr, du, an, la, cl in itertools.product(
            [7.0, 8.0, 10.0, 12.0, 15.0], [2.0, 1.5, 1.3], [2.0, 1.5, 1.3], [1.5, 2.0, 3.0],
            [np.inf, 5.0, 4.0, 3.0], [False, True]):
        out.append(replace(base, snr_min=sn, depth_ratio_max=dr, dur_ratio_max=du, anti_ratio=an,
                           local_anti_max=la, require_clean=cl))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["exclusion", "vetting", "final"])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--base", default="{}", help="Rule fields as a Python dict, for vetting and final")
    a = ap.parse_args()
    base = Rule(**eval(a.base, {"inf": np.inf}))
    if a.what == "exclusion":
        df = run(exclusion_variants(), "exclusion", a.workers, True, False)
    elif a.what == "vetting":
        df = run(vetting_variants(base), "vetting", a.workers, False, False)
    else:
        df = run([rule2.ROUND1, base], "final", 2, True, True)
    print(len(df), "variants written to", OUT)


if __name__ == "__main__":
    main()
