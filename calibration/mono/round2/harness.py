"""
Round 2 harness: one trial scored exactly under a rule (rule2.Rule), from light curves.

A trial is fully described by a definition (a dict, one CSV row):
  arm (A1, A2 or B), tic, t_cat, dur_cat_d, depth_cat, depth_cat_err, truth_P, truth_t0, mstar, rstar,
  multi (more than one planet candidate on the star), and for arm B the injected planet
  (inj_P, inj_t0, inj_k, inj_ar, inj_b, u1, u2).
make_defs_a / make_defs_b build definitions exactly as round 1 chose its trials; run_def scores one.

Compared with tune.py (which ranks variants from stored summaries) everything here is exact: alias
tables are evaluated at the predicted epochs and the ln P grid is not thinned.

Inputs  : a definitions CSV, or a sample CSV (--arm); MAST light curves; gaia_dr4/candidates_frozen.csv
Outputs : <out> CSV of trial rows
Usage   : python calibration/mono/round2/harness.py --defs defs.csv --out rows.csv [--workers 8]
          python calibration/mono/round2/harness.py --make-defs a|b --sample s.csv --out defs.csv
"""
import argparse
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import monocal as mc  # noqa: E402
import run_arm_a as ra  # noqa: E402
import run_arm_b as rb  # noqa: E402
import dump  # noqa: E402
import rule2  # noqa: E402
from rule2 import Rule  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
core, pl = mc.core, mc.pl
CACHE = ROOT / "calibration/mono/work/lc_cache"

# The rule proposed for the freeze (FROZEN_RULE.md). Not frozen until APPROVALS.md says so.
FROZEN = Rule(floor_sig=3.0, weak_ref_excludes=False, w_dur=0.5, w_dur_multi=1.0,
              snr_min=8.0, depth_ratio_max=1.7, dur_ratio_max=2.0, anti_ratio=1.5, require_clean=True)
STEP_FRAC = 0.125


class Boxes:
    """core.box_stats with the cumulative sums made once (t sorted)."""

    def __init__(self, t, f):
        self.t = t
        self.cs = np.concatenate([[0.0], np.cumsum(f)])
        self.cc = np.concatenate([[0.0], np.cumsum(core.local_cadence(t))])
        self.edges = core.gap_edges(t)

    def __call__(self, centers, dur_d):
        lo = np.searchsorted(self.t, centers - dur_d / 2)
        hi = np.searchsorted(self.t, centers + dur_d / 2)
        with np.errstate(invalid="ignore", divide="ignore"):
            depth = 1.0 - (self.cs[hi] - self.cs[lo]) / (hi - lo)
        return depth, (self.cc[hi] - self.cc[lo]) / dur_d, core.edge_distance(self.edges, centers)


def alias_table(t1, t2, box, tmin, tmax, dur, depth0, derr, sig, rule, ref_measured, multi, pmin,
                p_floor=mc.P_FLOOR, n_max=3000):
    """core.alias_table with the round-2 exclusion rule: status of every alias P_n = (t2 - t1)/n."""
    dT = t2 - t1
    w = rule2.w_steps(rule, STEP_FRAC, multi)
    off = np.arange(-w, w + 1) * dur * STEP_FRAC
    rows = []
    for n in range(1, n_max + 1):
        P = dT / n
        if P < p_floor:
            break
        ks = np.arange(int(np.floor((tmin - t1) / P)) - 1, int(np.ceil((tmax - t1) / P)) + 2)
        ks = ks[(ks != 0) & (ks != n)]
        status = "allowed"
        if ks.size:
            c = (t1 + ks * P)[:, None] + off[None, :]
            d, fr, ed = (x.reshape(c.shape) for x in box(c.ravel(), dur))
            ex = rule2.epoch_excludes(d, fr, ed, sig, depth0, derr, rule, ref_measured, multi, STEP_FRAC)
            d0, f0, e0 = d[:, w], fr[:, w], ed[:, w]
            tot = np.hypot(sig / np.sqrt(np.clip(f0, 1e-3, None)), derr)
            supp = rule2.clean(f0, e0) & (np.abs(d0 - depth0) < 3 * tot) & ~ex
            status = "excluded" if ex.any() else "supported" if supp.any() else "allowed"
        if np.isfinite(pmin) and P < pmin and status == "allowed":
            status = "allowed_below_pmin"
        rows.append((n, P, status))
    return rows


def call_from(events, t_cat, rule, ctx, truth_P):
    """Vetting, the alias event and the call for one list of events."""
    ok = [e for e in events if rule2.vet(e, rule)]
    out = dict(n_vetted=len(ok), false_event=any(not e.get("at_truth", False) for e in ok),
               repeat_detected=any(e.get("at_truth", False) for e in ok), t2=np.nan, call="no repeat",
               n_alias_ok=0, truth_in_alias_set=False)
    best = rule2.best_event(events, rule)
    if best:
        t1, t2 = sorted([t_cat, best["tc"]])
        # as st2_pipeline.step_aliases (round 1, correction 3): the alias table uses the second
        # event's depth error, not the reference's
        rows = alias_table(t1, t2, **dict(ctx, derr=best["depth_err"]))
        n_ok = sum(r[2] in ("allowed", "supported") for r in rows)
        dT = t2 - t1
        n = int(round(dT / truth_P))
        on = 1 <= n <= len(rows) and abs(dT - n * truth_P) < ctx["dur"] / 2
        out.update(t2=best["tc"], n_alias_ok=n_ok,
                   call="resolved" if n_ok == 1 else "narrowed" if n_ok else "all excluded",
                   truth_in_alias_set=bool(on and rows[n - 1][2] != "excluded"))
    return out


def score(prep, t_cat, truth_P, truth_t0, mstar, rstar, multi, rule=FROZEN, null_events=None):
    s = dump.summarize(prep, t_cat, truth_P, truth_t0)
    T, F, dur = prep["T"], prep["F"], prep["dur_ref"]
    d0, derr, meas, sig = prep["depth_ref"], prep["derr_ref"], s["ref_measured"], s["sig"]
    out = {k: s[k] for k in ("status", "ref_snr", "ref_measured", "depth_ref", "dur_ref_d", "ref_sector",
                             "ref_author", "ref_cadence_min", "n_sectors", "n_points", "baseline_d",
                             "n_events", "n_clean_repeats")}
    te = (rule2.epoch_excludes(s["tr_depth"], s["tr_frac"], s["tr_edge"], sig, d0, derr, rule, meas, multi)
          if s["tr_depth"].size else np.zeros(0, bool))
    out.update(truth_excluded=bool(te.any()), truth_n_epochs=int(te.size))
    ex = rule2.timeline_excludes(s["g_depth"], s["g_frac"], s["g_clean"], s["g_idx"], s["n_grid"], sig,
                                 d0, derr, rule, meas, multi)
    per = mc.period_grid(np.ptp(T), dur)
    out["frac_lnP_excluded"] = float(mc.excluded_periods(t_cat, per, s["c0"], s["step"], ex,
                                                         T.min(), T.max()).mean())
    ctx = dict(box=Boxes(T, F), tmin=T.min(), tmax=T.max(), dur=dur, depth0=d0, derr=derr, sig=sig,
               rule=rule, ref_measured=meas, multi=multi,
               pmin=core.pmin_from_duration(dur, mstar, rstar, d0))
    ev = [e for e in s["events"] if "error" not in e]
    out.update(call_from(ev, t_cat, rule, ctx, truth_P))
    if null_events is not None:
        null = np.array([e["tc"] for e in null_events])
        own = [e for e in ev if not e.get("at_truth", False) and null.size
               and np.min(np.abs(null - e["tc"])) < dur]
        keep = [e for e in ev if not any(e is o for o in own)]
        out.update({f"{k}_own": v for k, v in call_from(keep, t_cat, rule, ctx, truth_P).items()})
        out["null_event"] = any(rule2.vet(e, rule) for e in null_events)
        out["n_null_events"] = len(null_events)
    return out


# ------------------------------------------------------------------ trials from definitions

def _events_of(tic):
    c = pd.read_csv(ROOT / "gaia_dr4/candidates_frozen.csv")
    c = c[c.tic == tic]
    return pd.concat([c.t1_btjd, c.t2_btjd]).dropna().to_numpy(float)


def run_star(defs, cache=CACHE, rule=FROZEN):
    """All definitions of one star (they share its light curves)."""
    tic = int(defs[0]["tic"])
    sectors, _ = mc.load_star(tic, Path(cache))
    rows = []
    if not sectors:
        return [dict(d, status="no light curve") for d in defs]
    for d in defs:
        derr = d.get("depth_cat_err", np.nan)
        if d["arm"] in ("A1", "A2"):
            sec = sectors
            if d["arm"] == "A1":
                Traw = np.concatenate([t for _, _, t, _ in sectors])
                _, ep = mc.epochs(d["truth_t0"], d["truth_P"], Traw.min(), Traw.max())
                ref = ep[np.argmin(np.abs(ep - d["t_cat"]))]
                others = ep[np.abs(ep - ref) > d["dur_cat_d"] / 2]
                sec = [(s, a, *mc.mask_windows(t, f, others, 1.5 * d["dur_cat_d"])) for s, a, t, f in sectors]
                sec = [x for x in sec if x[2].size >= 50]
            null_events = None
        else:
            base = [(s, a, *mc.mask_windows(t, f, _events_of(tic), rb.EVENT_MASK_D)) for s, a, t, f in sectors]
            base = [x for x in base if x[2].size >= 50]
            p = dict(t0=d["inj_t0"], P=d["inj_P"], k=d["inj_k"], ar=d["inj_ar"], b=d["inj_b"])
            sec = rb.inject(base, p, (d["u1"], d["u2"]))
            null = mc.prepare(base, d["t_cat"], d["dur_cat_d"], d["depth_cat"])
            null_events = [] if null is None else [
                e for e in null["events"] if "error" not in e and abs(e["tc"] - d["t_cat"]) >= d["dur_cat_d"]]
        prep = mc.prepare(sec, d["t_cat"], d["dur_cat_d"], d["depth_cat"], derr)
        if prep is None:
            rows.append(dict(d, status="reference not in data"))
            continue
        rows.append(dict(d, **score(prep, d["t_cat"], d["truth_P"], d["truth_t0"], d["mstar"], d["rstar"],
                                    bool(d["multi"]), rule, null_events)))
    return rows


def make_defs_a(item, cache=CACHE):
    """Round 1's choice of the reference transit of a confirmed planet (run_arm_a.run_star)."""
    tic = int(item["tic"])
    rng = np.random.default_rng([ra.SEED, tic])
    sectors, _ = mc.load_star(tic, Path(cache))
    info = pl.tic_info(tic)
    if not sectors:
        return [], "no light curve"
    dur, depth, P0 = item["duration_h"] / 24, item["depth_ppm"] * 1e-6, item["period_d"]
    T, F = [], []
    for _, _, t, f in sectors:
        try:
            td, fd = pl.detrend(t, f, dur)
            T.append(td)
            F.append(fd)
        except Exception:
            pass
    T, F = np.concatenate(T), np.concatenate(F)
    o = np.argsort(T)
    perr = item.get("period_err")
    half = max(5 * (perr if perr is not None and np.isfinite(perr) else 0.0), 2e-4 * P0)
    eph = ra._s10().refine_ephemeris(T[o], F[o], P0, item["t0"], dur, half)
    P, t0 = eph["P"], eph["t0"]
    Traw = np.sort(np.concatenate([t for _, _, t, _ in sectors]))
    _, ep = mc.epochs(t0, P, Traw.min(), Traw.max())
    cl = mc.clean_mask(Traw, ep, dur)
    if cl.sum() < ra.MIN_CLEAN:
        return [], "under 3 clean transits"
    t_cat = float(ep[cl][0]) + rng.normal(0, ra.T_CAT_ERR_D)
    derr = item.get("depth_err")
    base = dict(tic=tic, toi=item["toi"], t_cat=t_cat, dur_cat_d=dur, depth_cat=depth,
                depth_cat_err=derr * 1e-6 if derr is not None and np.isfinite(derr) else np.nan,
                truth_P=P, truth_t0=t0, mstar=info.get("mstar", np.nan), rstar=info.get("rstar", np.nan),
                multi=bool(item["n_toi_on_star"] > 1), tmag=item["tmag"], at_grid_edge=bool(eph["at_grid_edge"]))
    return [dict(base, arm="A1"), dict(base, arm="A2")], ""


def make_defs_b(item, cache=CACHE, seed=rb.SEED, draws=rb.DRAWS):
    """Round 1's injections (run_arm_b.run_host) for one host, with the given seed."""
    tic = int(item["tic"])
    rng = np.random.default_rng([seed, tic])
    sectors, _ = mc.load_star(tic, Path(cache))
    info = pl.tic_info(tic)
    if not sectors:
        return [], "no light curve"
    sectors = [(s, a, *mc.mask_windows(t, f, _events_of(tic), rb.EVENT_MASK_D)) for s, a, t, f in sectors]
    sectors = [x for x in sectors if x[2].size >= 50]
    mstar, rstar = info.get("mstar", np.nan), info.get("rstar", np.nan)
    u1, u2, _ = mc.limb_dark(info.get("teff", np.nan), mstar, rstar)
    if not (np.isfinite(mstar) and np.isfinite(rstar) and mstar > 0 and rstar > 0):
        mstar, rstar = 1.0, 1.0
    out = []
    for d in range(draws):
        p = rb.draw(rng, sectors, mstar, rstar)
        if p is None:
            continue
        t_cat = p["t_ref"] + rng.normal(0, rb.T_CAT_ERR_D)
        out.append(dict(arm="B", tic=tic, draw=d, t_cat=t_cat, dur_cat_d=p["dur_d"], depth_cat=p["depth"],
                        depth_cat_err=np.nan, truth_P=p["P"], truth_t0=p["t0"], mstar=mstar, rstar=rstar,
                        multi=False, tmag=info.get("tmag", np.nan), inj_P=p["P"], inj_t0=p["t0"],
                        inj_k=p["k"], inj_ar=p["ar"], inj_b=p["b"], u1=u1, u2=u2))
    return out, ""


def _star(args):
    fn, a = args
    try:
        return fn(*a)
    except Exception as e:
        return dict(error=f"{type(e).__name__}: {str(e)[:300]}", args=str(a[0])[:200])


def _pool(fn, arglist, workers):
    with mp.get_context("spawn").Pool(workers) as pool:
        for i, r in enumerate(pool.imap(_star, [(fn, a) for a in arglist]), 1):
            print(f"[{i}/{len(arglist)}]", flush=True)   # progress only, no outcome
            yield r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--defs")
    ap.add_argument("--make-defs", choices=["a", "b"])
    ap.add_argument("--sample")
    ap.add_argument("--seed", type=int, default=rb.SEED)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache", default=str(CACHE))
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    if a.make_defs:
        items = pd.read_csv(a.sample).to_dict("records")
        fn = make_defs_a if a.make_defs == "a" else make_defs_b
        args = [(it, a.cache) if a.make_defs == "a" else (it, a.cache, a.seed) for it in items]
        defs, notes = [], []
        for it, r in zip(items, _pool(fn, args, a.workers)):
            if isinstance(r, dict):
                notes.append(dict(tic=it["tic"], note=r["error"]))
            else:
                defs += r[0]
                notes.append(dict(tic=it["tic"], note=r[1]))
        pd.DataFrame(defs).to_csv(a.out, index=False)
        pd.DataFrame(notes).to_csv(Path(a.out).with_suffix(".notes.csv"), index=False)
        return
    defs = pd.read_csv(a.defs)
    groups = [(g.to_dict("records"), a.cache) for _, g in defs.groupby("tic", sort=True)]
    rows, errors = [], []
    for r in _pool(run_star, groups, a.workers):
        if isinstance(r, dict):
            errors.append(r)
        else:
            rows += r
    pd.DataFrame(rows).to_csv(a.out, index=False)
    if errors:
        pd.DataFrame(errors).to_csv(Path(a.out).with_suffix(".errors.csv"), index=False)
        print(len(errors), "stars failed; see the .errors.csv file")


if __name__ == "__main__":
    main()
