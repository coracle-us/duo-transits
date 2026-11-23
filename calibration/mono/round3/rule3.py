"""
Round 3: the second-transit call with its tunable parts (prereg_round3.yaml, T1-T5).

The single-transit exclusion rule is round 2's frozen rule (round2/harness.FROZEN) and is not
changed here: a box excludes, is clean, is "deep", exactly as there. Round 3 changes only
  T1  how wide a window around a predicted epoch of an alias must be free of a deep box,
  T3  how the minimum period from the duration enters the call,
  T5  which second events are accepted.
ROUND2 (all round-3 switches off) reproduces round2/harness.alias_table and call_from.

No file input or output. Usage: imported by tune3.py, harness3.py and the tests.
"""
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "round2"))
import rule2  # noqa: E402

STEP_FRAC = 0.125
P_FLOOR = 1.0
BASE = rule2.Rule(floor_sig=3.0, weak_ref_excludes=False, w_dur=0.5, w_dur_multi=1.0, snr_min=8.0,
                  depth_ratio_max=1.7, dur_ratio_max=2.0, anti_ratio=1.5, require_clean=True)


@dataclass(frozen=True)
class Rule3:
    # T1 period tolerance. Half-width of the window at epoch k of alias n, in days:
    #   W dur + c_drift * |k - k_near| * sigma_P,   sigma_P = hypot(sigma_t1, sigma_t2) / n,
    # k_near the nearer of 0 and n, sigma_t = max(sig_floor_min minutes, sig_dur_snr * dur_event / snr).
    c_drift: float = 0.0
    sig_floor_min: float = 15.0
    sig_dur_snr: float = 0.0
    cap_dur: float = np.inf        # a window wider than this many durations: the epoch is unconstrained
    need_window: bool = False      # an uncovered box inside the window: the epoch is unconstrained
    use_ref_time: bool = False     # first transit time = the pipeline's fit, if the reference is measured
    # T3 minimum period from the duration: "hard" (round 2), "off", or "cat" (catalogued duration)
    pmin_mode: str = "hard"
    pmin_factor: float = 1.0       # the limit is pmin * factor
    # T5 event vetting beyond round 2's R4
    lowref_snr: float = 0.0        # if the reference SNR is below this ...
    lowref_event_snr: float = 0.0  # ... the event needs at least this SNR
    sec_box_snr_min: float = 0.0   # event depth over the scatter of same-duration boxes in its sector
    # two or more passing events that fit no common period with the reference: none of them is used.
    # Periods tried: (t_nearest_event - t_ref)/n down to max(1 d, cp_pmin_frac * pmin); an event fits
    # if it lies within cp_tol_dur durations of an epoch.
    common_period: bool = False
    cp_pmin_frac: float = 0.0
    cp_tol_dur: float = 1.0
    local_anti_max: float = np.inf


ROUND2 = Rule3()


def vet3(ev, r, ref_snr):
    if not rule2.vet(ev, BASE):
        return False
    if ref_snr < r.lowref_snr and ev["snr"] < r.lowref_event_snr:
        return False
    if r.sec_box_snr_min > 0 and not (ev.get("sec_box_snr", np.nan) >= r.sec_box_snr_min):
        return False
    return ev.get("local_anti_snr", 0.0) <= r.local_anti_max


def best_event3(events, r, ref_snr):
    ok = [e for e in events if vet3(e, r, ref_snr)]
    if not ok:
        return None, ok
    return max(ok, key=lambda e: (e["edge_dist_d"] >= rule2.CLEAN_EDGE_D and e["frac"] >= rule2.CLEAN_FRAC,
                                  e["snr"])), ok


def fits_common_period(times, t_ref, dur, pmin, r):
    """True if some period puts every time within cp_tol_dur durations of an epoch counted from t_ref."""
    dt = np.asarray(times, float) - t_ref
    if dt.size < 2:
        return True
    base = np.abs(dt).min()
    floor = max(P_FLOOR, r.cp_pmin_frac * pmin if np.isfinite(pmin) else P_FLOOR)
    n = np.arange(1, max(1, int(base / floor)) + 1)
    P = base / n
    x = dt[:, None] / P[None, :]
    return bool((np.abs(x - np.round(x)) * P[None, :] < r.cp_tol_dur * dur).all(axis=0).any())


def sigma_t(dur_event, snr, r):
    return max(r.sig_floor_min / 1440.0, r.sig_dur_snr * dur_event / max(snr, 1e-9))


def alias_table3(t1, t2, sig_dt, box, tmin, tmax, dur, depth0, derr, sig, ref_measured, multi, pmin, r,
                 n_max=3000, chunk=64):
    """Status of every alias P_n = (t2 - t1)/n above P_FLOOR: (n, P, status).

    An epoch excludes if its window is not wider than the cap, its central box is clean and excludes
    (round 2's box rule), and no covered box in the window is as deep as its threshold. sig_dt is
    the error of t2 - t1 in days.
    """
    dT = t2 - t1
    step = dur * STEP_FRAC
    w0 = rule2.w_steps(BASE, STEP_FRAC, multi)
    rows = []
    for n in range(1, n_max + 1):
        P = dT / n
        if P < P_FLOOR:
            break
        ks = np.arange(int(np.floor((tmin - t1) / P)) - 1, int(np.ceil((tmax - t1) / P)) + 2)
        ks = ks[(ks != 0) & (ks != n)]
        status = "allowed"
        if ks.size:
            far = np.minimum(np.abs(ks), np.abs(ks - n))
            hw = w0 + np.ceil(r.c_drift * far * (sig_dt / n) / step - 1e-9).astype(int)
            ok = hw <= r.cap_dur / STEP_FRAC
            o = np.argsort(hw[ok], kind="stable")
            ks, hw = ks[ok][o], hw[ok][o]
            excluded, supported = False, False
            for i in range(0, ks.size, chunk):
                k, h = ks[i:i + chunk], hw[i:i + chunk]
                m = int(h.max())
                off = np.arange(-m, m + 1)
                c = (t1 + k * P)[:, None] + off[None, :] * step
                d, fr, ed = (x.reshape(c.shape) for x in box(c.ravel(), dur))
                err = sig / np.sqrt(np.clip(fr, 1e-3, None))
                thr = rule2.threshold(err, depth0, derr)
                inwin = np.abs(off)[None, :] <= h[:, None]
                covered = fr >= rule2.MIN_FRAC
                ex0 = rule2.box_excludes(d[:, m], fr[:, m], sig, depth0, derr, BASE, ref_measured) \
                    & rule2.clean(fr[:, m], ed[:, m])
                deep = (covered & (d >= thr) & inwin).any(axis=1)
                ex = ex0 & ~deep
                if r.need_window:
                    ex &= ~(inwin & ~covered).any(axis=1)
                if ex.any():
                    excluded = True
                    break
                tot = np.hypot(err[:, m], derr)
                supported |= bool((rule2.clean(fr[:, m], ed[:, m]) & (np.abs(d[:, m] - depth0) < 3 * tot)).any())
            status = "excluded" if excluded else "supported" if supported else "allowed"
        if r.pmin_mode != "off" and np.isfinite(pmin) and P < pmin * r.pmin_factor and status == "allowed":
            status = "allowed_below_pmin"
        rows.append((n, P, status))
    return rows


def call3(events, tr, r, box, truth_P):
    """Vetting, the alias event and the call for one list of events of one stored trial `tr`."""
    best, ok = best_event3(events, r, tr["ref_snr"])
    out = dict(n_vetted=len(ok), false_event=any(not e.get("at_truth", False) for e in ok),
               repeat_detected=any(e.get("at_truth", False) for e in ok), t2=np.nan, call="no repeat",
               n_alias_ok=0, truth_in_alias_set=False, event_true=False, frac_alias_excluded=np.nan,
               n_alias=0)
    d = tr["d"]
    if best is not None and r.common_period and len(ok) >= 2 and not fits_common_period(
            [e["tc"] for e in ok], d["t_cat"], tr["dur"], tr["pmin"], r):
        out.update(n_vetted=0, false_event=False, repeat_detected=False, call="events fit no common period")
        return out
    if best is None:
        return out
    use_ref = r.use_ref_time and tr["ref_measured"] and np.isfinite(tr["ref_tc"])
    tref = tr["ref_tc"] if use_ref else d["t_cat"]
    s_ref = sigma_t(tr["ref_dur"] if np.isfinite(tr["ref_dur"]) else tr["dur"], tr["ref_snr"], r) if use_ref \
        else max(r.sig_floor_min / 1440.0, 15 / 1440.0)
    s_ev = sigma_t(best.get("dur_d", tr["dur"]), best["snr"], r)
    t1, t2 = sorted([tref, best["tc"]])
    pmin = tr["pmin_cat"] if r.pmin_mode == "cat" else tr["pmin"]
    rows = alias_table3(t1, t2, float(np.hypot(s_ref, s_ev)), box, tr["tmin"], tr["tmax"], tr["dur"],
                        tr["depth0"], best["depth_err"], tr["sig"], tr["ref_measured"], bool(d["multi"]), pmin, r)
    n_ok = sum(x[2] in ("allowed", "supported") for x in rows)
    dT = t2 - t1
    n = int(round(dT / truth_P))
    on = 1 <= n <= len(rows) and abs(dT - n * truth_P) < tr["dur"] / 2
    out.update(t2=best["tc"], n_alias_ok=n_ok, n_alias=len(rows),
               call="resolved" if n_ok == 1 else "narrowed" if n_ok else "all excluded",
               truth_in_alias_set=bool(on and rows[n - 1][2] != "excluded"),
               event_true=bool(best.get("at_truth", False)),
               frac_alias_excluded=float(np.mean([x[2] == "excluded" for x in rows])) if rows else np.nan)
    return out
