"""
The catalogue run, one row at a time (prereg_catalogue.yaml, what_is_computed). No rule is defined
here: the single-transit exclusion is round2/harness.FROZEN, event vetting and the alias table are
round3/harness3.FROZEN3, through the frozen functions of calibration/mono.

search_row : steps 1 and 2 - reference, detrending, the excluded periods, every event.
table_rows : step 4 - the alias table through two times.
finish_row : step 4 for one row, given the events both vetters accepted.

No file input or output. Used by run_catalogue.py and the tests.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "calibration/mono/round3"))
sys.path.insert(0, str(ROOT / "calibration/mono/round2"))
import harness3 as h3  # noqa: E402  (frozen)
import harness as h2  # noqa: E402  (frozen)
import rule2  # noqa: E402
import rule3  # noqa: E402
import score2  # noqa: E402

mc, core, dump, dump3 = h2.mc, h2.core, h2.dump, h3.dump3
SIG_T_D = 15.0 / 1440.0                       # T1: every transit time has an error of 15 min
SIG_DT_D = float(np.hypot(SIG_T_D, SIG_T_D))
LIVE = ("allowed", "supported")


def event_id(row_id, ev):
    return f"{row_id}|s{int(ev['sector'])}|{ev['tc']:.4f}"


def _runs(per, flag):
    """[(p_lo, p_hi)] for each run of consecutive grid periods with flag set."""
    if not flag.any():
        return []
    d = np.diff(np.concatenate([[0], flag.astype(int), [0]]))
    return [(float(per[a]), float(per[b - 1])) for a, b in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1))]


def search_row(d, sectors, other_times=()):
    """Steps 1 and 2 for one definition. Returns (row, events, ctx); ctx is None if the reference
    is not in the data. other_times: catalogued times of the star's other scored rows."""
    multi = bool(d["multi"])
    derr_cat = d.get("depth_cat_err", np.nan)
    prep = mc.prepare(sectors, d["t_cat"], d["dur_cat_d"], d["depth_cat"],
                      derr_cat if derr_cat is not None else np.nan)
    if prep is None:
        return dict(row_id=d["row_id"], status="reference not in data"), [], None
    s = dump.summarize(prep, d["t_cat"])
    T, F, dur = prep["T"], prep["F"], prep["dur_ref"]
    d0, derr, meas, sig = prep["depth_ref"], prep["derr_ref"], s["ref_measured"], s["sig"]
    ex = rule2.timeline_excludes(s["g_depth"], s["g_frac"], s["g_clean"], s["g_idx"], s["n_grid"], sig,
                                 d0, derr, h2.FROZEN, meas, multi)
    p_test = float(max(T.max() - d["t_cat"], d["t_cat"] - T.min()))
    per = mc.period_grid(np.ptp(T), dur, hi=max(mc.P_HI, p_test))
    exP = mc.excluded_periods(d["t_cat"], per, s["c0"], s["step"], ex, T.min(), T.max())
    cal = per < mc.P_HI                       # the calibrated grid (1 to 1000 d), as round 2
    tst = per <= p_test
    row = dict(row_id=d["row_id"], status="ok", ref_measured=bool(meas), ref_snr=s["ref_snr"],
               depth_ref=d0, depth_ref_err=derr, dur_ref_d=dur, ref_sector=s["ref_sector"],
               ref_author=s["ref_author"], n_sectors=s["n_sectors"], n_points=s["n_points"],
               baseline_d=s["baseline_d"], tmin=float(T.min()), tmax=float(T.max()),
               sectors=" ".join(f"{sec}:{a}" for sec, a, _, _ in sectors if sec in prep["det"]),
               frac_lnP_excluded=float(exP[cal].mean()),
               p_longest_testable_d=p_test,
               frac_lnP_excluded_testable=float(exP[tst].mean()) if tst.any() else np.nan,
               p_false_exclusion=score2.p_false_exclusion(s["ref_snr"], meas),
               pmin_d=core.pmin_from_duration(dur, d["mstar"], d["rstar"], d0), n_events=0, n_passing=0)
    events = []
    t2c = d.get("t2_cat", np.nan)
    for e in s["events"]:
        if "error" in e:
            continue
        e = dict(e)
        e["sec_box_snr"] = dump3.sector_box_snr(prep["det"], e, dur)
        e["passes"] = bool(rule3.vet3(e, h3.FROZEN3, s["ref_snr"]))
        e["at_t2_cat"] = bool(np.isfinite(t2c) and abs(e["tc"] - t2c) < dur)
        e["at_other_row"] = bool(any(abs(e["tc"] - x) < dur for x in other_times))
        e["event_id"] = event_id(d["row_id"], e)
        e["row_id"] = d["row_id"]
        events.append(e)
    row.update(n_events=len(events), n_passing=sum(e["passes"] for e in events))
    ctx = dict(T=T, F=F, det=prep["det"], dur=dur, depth0=d0, derr=derr, sig=sig, ref_measured=bool(meas),
               ref_snr=s["ref_snr"], tmin=float(T.min()), tmax=float(T.max()), pmin=row["pmin_d"], multi=multi,
               excluded_runs=_runs(per[cal | tst], exP[cal | tst]), box=h2.Boxes(T, F))
    return row, events, ctx


def table_rows(ctx, ta, tb, derr_event):
    """T1 alias table through two transit times: [(n, P, status)]."""
    t1, t2 = sorted([float(ta), float(tb)])
    return rule3.alias_table3(t1, t2, SIG_DT_D, ctx["box"], ctx["tmin"], ctx["tmax"], ctx["dur"], ctx["depth0"],
                              derr_event, ctx["sig"], ctx["ref_measured"], ctx["multi"], ctx["pmin"], h3.FROZEN3)


def summarize_table(rows):
    n_ok = sum(r[2] in LIVE for r in rows)
    return dict(n_alias=len(rows), n_alias_ok=n_ok,
                call="resolved" if n_ok == 1 else "narrowed" if n_ok else "all excluded",
                frac_alias_excluded=float(np.mean([r[2] == "excluded" for r in rows])) if rows else np.nan)


def fits_event(ctx, t1, t2, n, t_ev):
    """Does alias n of the table through t1 < t2 put an epoch on t_ev, within that epoch's T1 window?
    'yes', 'no', or 'unconstrained' (the window at the nearest epoch is over the cap)."""
    P = (t2 - t1) / n
    k = int(round((t_ev - t1) / P))
    step = ctx["dur"] * rule3.STEP_FRAC
    w0 = rule2.w_steps(rule3.BASE, rule3.STEP_FRAC, ctx["multi"])
    h = w0 + int(np.ceil(h3.FROZEN3.c_drift * min(abs(k), abs(k - n)) * (SIG_DT_D / n) / step - 1e-9))
    if h > h3.FROZEN3.cap_dur / rule3.STEP_FRAC:
        return "unconstrained"
    return "yes" if abs(t_ev - (t1 + k * P)) <= h * step else "no"


def _depth_err_at(ctx, events, t, dur):
    """The pipeline's depth error of a transit at catalogued time t: that of its event within half a
    duration, else core.refine_event at t in the sector that holds t; None if neither exists."""
    near = [e for e in events if abs(e["tc"] - t) < dur / 2]
    if near:
        return max(near, key=lambda e: e["snr"])["depth_err"], "event"
    for td, fd in ctx["det"].values():
        if td.min() - 0.5 <= t <= td.max() + 0.5:
            r = core.refine_event(td, fd, t, dur)
            if r:
                return r["depth_err"], "refit"
    return None, "none"


def finish_row(d, ctx, events, accepted_ids):
    """Step 4. events: the row's events from search_row; accepted_ids: event_ids both vetters accepted.
    Returns (columns, tables) with tables = [(kind, t1, t2, [(n, P, status, fits_third)])]."""
    acc = [e for e in events if e["passes"] and e["event_id"] in accepted_ids]
    out, tables = dict(n_accepted=len(acc)), []
    t_cat, t2c, dur = d["t_cat"], d.get("t2_cat", np.nan), ctx["dur"]
    if np.isfinite(t2c):                                    # Hawthorn duotransit
        derr2, how = _depth_err_at(ctx, events, t2c, dur)
        out["second_depth_err_from"] = how
        out["second_recovered"] = any(e["passes"] and abs(e["tc"] - t2c) < dur / 2 for e in events)
        third = [e for e in acc if abs(e["tc"] - t2c) >= dur]
        best, _ = rule3.best_event3(third, h3.FROZEN3, ctx["ref_snr"])
        rows = table_rows(ctx, t_cat, t2c, ctx["derr"] if derr2 is None else derr2)
        t1, t2 = sorted([t_cat, t2c])
        fit = [fits_event(ctx, t1, t2, n, best["tc"]) if best else "" for n, _, _ in rows]
        tables.append(("duo", t1, t2, [(n, P, st, f) for (n, P, st), f in zip(rows, fit)]))
        out.update(summarize_table(rows))
        if best:
            r3 = table_rows(ctx, t_cat, best["tc"], best["depth_err"])
            a, b = sorted([t_cat, best["tc"]])
            tables.append(("ref_third", a, b, [(n, P, st, "") for n, P, st in r3]))
            live_fit = sum(st in LIVE and f == "yes" for (_, _, st), f in zip(rows, fit))
            out.update(third_event_id=best["event_id"], t3=best["tc"], n_alias_ok_fitting_third=live_fit,
                       **{k + "_ref_third": v for k, v in summarize_table(r3).items()})
        return out, tables
    best, _ = rule3.best_event3(acc, h3.FROZEN3, ctx["ref_snr"])
    if best is None:
        out.update(call="no repeat", n_alias=0, n_alias_ok=0, frac_alias_excluded=np.nan)
        return out, tables
    rows = table_rows(ctx, t_cat, best["tc"], best["depth_err"])
    t1, t2 = sorted([t_cat, best["tc"]])
    tables.append(("single", t1, t2, [(n, P, st, "") for n, P, st in rows]))
    out.update(alias_event_id=best["event_id"], t2=best["tc"], **summarize_table(rows))
    return out, tables
