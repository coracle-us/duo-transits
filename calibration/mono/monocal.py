"""
Single-transit calibration harness (calibration/mono/prereg_mono_calibration.yaml).

One trial = one star's raw sectors, with one transit taken as the catalogued single transit:
  1. the caller changes the raw flux (inject a planet, or remove a known planet's other transits);
  2. every sector is detrended and searched exactly as st2_pipeline does (pl.detrend, pl.sector_events,
     pl.passes with the pipeline's default cuts, core.alias_table);
  3. the trial records whether the true period was excluded, whether a repeat was found, and how
     much of ln P between 1 and 1000 d is excluded.

Step 3's single-transit excluded set is not a function st2core has; it is defined here from
core.epoch_verdicts (see exclusion_timeline / excluded_periods) and listed in PROPOSED_v1.1.md.

No network in this module except through load_star(). Usage: see run_arm_b.py.
"""
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "salinas-second-transits"))
import st2_pipeline as pl  # noqa: E402
import st2core as core  # noqa: E402

# st2_pipeline.main() defaults; the pipeline has them only as argparse defaults
CUTS = SimpleNamespace(snr_min=7.0, depth_lo=0.5, depth_hi=2.0, dur_lo=0.5, dur_hi=2.0, anti_ratio=1.5)
P_FLOOR = 1.0
P_LO, P_HI = 1.0, 1000.0
LIVE = ("allowed", "supported", "allowed_below_pmin")


# ------------------------------------------------------------------ transit model

def transit_duration_d(P_d, mstar, rstar, k, b):
    """First-to-fourth-contact duration of a circular orbit (days)."""
    ar = core.a_over_rstar(P_d, mstar, rstar)
    x = np.sqrt(max((1 + k) ** 2 - b ** 2, 0.0)) / np.sqrt(max(ar ** 2 - b ** 2, 1e-12))
    return P_d / np.pi * np.arcsin(min(x, 1.0))


def transit_flux(t, t0, P_d, k, ar, b, u, cad_d):
    """batman quadratic-limb-darkening model, integrated over the cadence."""
    import batman
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp, pm.a = t0, P_d, k, ar
    pm.inc = float(np.degrees(np.arccos(min(b / ar, 1.0))))
    pm.ecc, pm.w, pm.u, pm.limb_dark = 0.0, 90.0, list(u), "quadratic"
    ss = max(1, int(round(cad_d / (2.0 / 1440))))
    return batman.TransitModel(pm, np.ascontiguousarray(t, float), supersample_factor=ss,
                               exp_time=cad_d).light_curve(pm)


_LD = None


def limb_dark(teff, mstar, rstar):
    """Quadratic TESS coefficients (Claret 2017, ATLAS, solar metallicity, LSM) at the nearest grid
    point in Teff and log g. log g comes from the TIC mass and radius. Returns (u1, u2, how)."""
    global _LD
    if _LD is None:
        import pandas as pd
        _LD = pd.read_csv(Path(__file__).resolve().parent / "data/claret2017_table25_Z0_xi2.csv")
    how = "claret"
    if not (np.isfinite(teff) and teff > 0):
        teff, how = 5750.0, "claret, solar Teff assumed"
    if np.isfinite(mstar) and np.isfinite(rstar) and mstar > 0 and rstar > 0:
        logg = 4.438 + np.log10(mstar / rstar ** 2)
    else:
        logg, how = 4.5, how + ", log g 4.5 assumed"
    if teff < _LD.Teff.min():
        how += ", below the grid (3500 K)"
    d = ((_LD.Teff - teff) / 250.0) ** 2 + ((_LD.logg - logg) / 0.5) ** 2
    r = _LD.loc[d.idxmin()]
    return float(r.aLSM), float(r.bLSM), how


def mask_windows(t, f, centers, half_d):
    keep = np.ones(t.size, bool)
    for c in np.atleast_1d(centers):
        keep &= np.abs(t - c) > half_d
    return t[keep], f[keep]


def epochs(t0, P_d, tmin, tmax):
    k = np.arange(int(np.floor((tmin - t0) / P_d)) - 1, int(np.ceil((tmax - t0) / P_d)) + 2)
    return k, t0 + k * P_d


def clean_mask(t, ep, dur_d):
    """Epochs with the pipeline's clean coverage: CLEAN_FRAC covered, CLEAN_EDGE_D from a gap."""
    _, frac, _ = core.box_stats(t, np.ones(t.size), ep, dur_d)
    return (frac >= core.CLEAN_FRAC) & (core.edge_distance(core.gap_edges(t), ep) >= core.CLEAN_EDGE_D)


# ------------------------------------------------------------------ single-transit excluded set

def exclusion_timeline(t, f, dur_d, depth0, depth0_err):
    """core.epoch_verdicts on a time grid of step dur/8: (first centre, step, 'excludes' mask)."""
    step = dur_d / 8
    centers = np.arange(t.min(), t.max() + step, step)
    v = core.epoch_verdicts(centers, t, f, dur_d, depth0, depth0_err)
    return float(centers[0]), step, v["verdict"] == "excludes"


def period_grid(baseline_d, dur_d, lo=P_LO, hi=P_HI):
    """Log-spaced periods; a step moves the farthest epoch by at most a quarter duration."""
    dln = dur_d / (4 * max(baseline_d, dur_d))
    return np.exp(np.arange(np.log(lo), np.log(hi), dln))


def excluded_periods(t_ref, periods, c0, step, excl, tmin, tmax, max_elems=4_000_000):
    """True for each period with an epoch t_ref + k P (k != 0) on an 'excludes' grid point."""
    periods = np.asarray(periods, float)
    kmin = np.floor((tmin - t_ref) / periods).astype(np.int64)
    kmax = np.ceil((tmax - t_ref) / periods).astype(np.int64)
    cnt = kmax - kmin + 1
    out = np.zeros(periods.size, bool)
    cum = np.concatenate([[0], np.cumsum(cnt)])
    i = 0
    while i < periods.size:
        j = max(i + 1, int(np.searchsorted(cum, cum[i] + max_elems, side="right")) - 1)
        c = cnt[i:j]
        idx = np.repeat(np.arange(j - i), c)
        k = np.arange(c.sum()) - np.repeat(cum[i:j] - cum[i], c) + np.repeat(kmin[i:j], c)
        g = np.rint((t_ref + k * periods[i:j][idx] - c0) / step).astype(np.int64)
        ok = (k != 0) & (g >= 0) & (g < excl.size)
        out[i:j] = np.bincount(idx[ok], excl[g[ok]], minlength=j - i) > 0
        i = j
    return out


# ------------------------------------------------------------------ one trial

def run_trial(sectors, t_cat, dur_cat_d, depth_cat, mstar=np.nan, rstar=np.nan, truth_P=None,
              truth_t0=None, with_grid=True):
    """sectors: [(sector, author, t, f)] raw, normalised (as pl.load_sector_lc returns them).

    t_cat, dur_cat_d, depth_cat play the catalogue's time, duration and depth of the single transit.
    truth_P / truth_t0: the real ephemeris, used only to score the trial.
    """
    out = dict(status="ok")
    ref_i = [i for i, (_, _, t, _) in enumerate(sectors) if t.min() - 0.5 <= t_cat <= t.max() + 0.5]
    if not ref_i:
        return dict(status="reference not in data")
    ri = ref_i[0]
    s, a, t, f = sectors[ri]
    tr, fr = pl.detrend(t, f, dur_cat_d)
    ref = core.refine_event(tr, fr, t_cat, dur_cat_d)
    # as st2_pipeline.search_one: our own depth if it is solid, else the catalogue's
    if ref and ref["snr"] >= 5 and 0.33 < ref["depth"] / depth_cat < 3:
        depth_ref, dur_ref, derr_ref = ref["depth"], ref["dur_d"], ref["depth_err"]
    else:
        depth_ref, dur_ref, derr_ref = depth_cat, dur_cat_d, 0.2 * depth_cat
    out.update(ref_snr=ref["snr"] if ref else np.nan, ref_measured=bool(ref and depth_ref == ref["depth"]),
               depth_ref=depth_ref, dur_ref_d=dur_ref, ref_sector=s, ref_author=a)

    T, F, events = [tr], [fr], []
    for i, (s, a, t, f) in enumerate(sectors):
        if i == ri:
            td, fd = tr, fr
        else:
            try:
                td, fd = pl.detrend(t, f, dur_ref)
            except Exception as e:  # as search_sector: a sector that fails is skipped
                events.append(dict(sector=s, error=str(e)))
                continue
            T.append(td)
            F.append(fd)
        for ev in pl.sector_events(td, fd, s, a, dur_ref, depth_ref):
            if i == ri and abs(ev["tc"] - t_cat) < dur_ref:
                continue  # the reference itself
            events.append(ev)
    T, F = np.concatenate(T), np.concatenate(F)
    o = np.argsort(T)
    T, F = T[o], F[o]
    passing = [e for e in events if pl.passes(e, CUTS)]
    out.update(n_sectors=len(sectors), n_points=int(T.size), baseline_d=float(np.ptp(T)),
               n_events=len(events), n_passing=len(passing))

    pmin = core.pmin_from_duration(dur_ref, mstar, rstar, depth_ref)
    if passing:  # as st2_pipeline.alias_events: best clean event, else best event
        key = lambda e: (e["edge_dist_d"] >= core.CLEAN_EDGE_D and e["frac"] >= core.CLEAN_FRAC, e["snr"])
        best = max(passing, key=key)
        rows = core.alias_table(t_cat, best["tc"], dur_ref, depth_ref, derr_ref, T, F,
                                p_floor=P_FLOOR, pmin_density=pmin)
        ok = [r for r in rows if r["status"] in ("allowed", "supported")]
        out.update(t2=best["tc"], n_aliases=len(rows), n_alias_ok=len(ok),
                   call="resolved" if len(ok) == 1 else "narrowed" if ok else "all excluded")
    else:
        rows = []
        out.update(t2=np.nan, n_aliases=0, n_alias_ok=0, call="no repeat")

    if truth_P is not None:
        t0 = t_cat if truth_t0 is None else truth_t0
        k, ep = epochs(t0, truth_P, T.min(), T.max())
        other = np.abs(ep - t_cat) > dur_ref / 2
        cl = clean_mask(T, ep, dur_ref) & other
        out["n_clean_repeats"] = int(cl.sum())
        out["repeat_detected"] = bool(any(np.min(np.abs(ep[other] - e["tc"])) < dur_ref / 2
                                          for e in passing)) if other.any() else False
        # the single-transit test: epochs counted from the catalogued time, as a user would
        k1, ep1 = epochs(t_cat, truth_P, T.min(), T.max())
        v = core.epoch_verdicts(ep1[k1 != 0], T, F, dur_ref, depth_ref, derr_ref)
        out["truth_excluded_single"] = bool((v["verdict"] == "excludes").any())
        out["truth_n_clean_epochs"] = int(np.isin(v["verdict"], ["excludes", "supports", "shallow", "deeper"]).sum())
        if rows:
            dT = out["t2"] - t_cat
            n = int(round(dT / truth_P))
            on_grid = n >= 1 and n <= len(rows) and abs(dT - n * truth_P) < dur_ref / 2
            out["truth_alias_status"] = rows[n - 1]["status"] if on_grid else "not an alias"
            out["truth_in_alias_set"] = bool(on_grid and rows[n - 1]["status"] in LIVE)

    if with_grid:
        c0, step, excl = exclusion_timeline(T, F, dur_ref, depth_ref, derr_ref)
        per = period_grid(np.ptp(T), dur_ref)
        ex = excluded_periods(t_cat, per, c0, step, excl, T.min(), T.max())
        out.update(frac_lnP_excluded=float(ex.mean()), n_grid_periods=int(per.size))
        if truth_P is not None:
            out["truth_excluded_grid"] = bool(ex[min(np.searchsorted(per, truth_P), per.size - 1)])
    return out


# ------------------------------------------------------------------ data access (network)

def load_star(tic, cache_dir):
    """Raw sectors of one star, one product per sector, as the pipeline chooses them."""
    sr, prods = pl.lc_products(int(tic))
    best = pl.best_per_sector(prods)
    sectors, errors = [], []
    for s, row in sorted(best.items()):
        if s in pl.SPECIAL_SECTORS:
            continue
        try:
            t, f = pl.load_sector_lc(sr, row, cache_dir)
            if t.size >= 50:
                sectors.append((s, row["author"], t, f))
        except Exception as e:
            errors.append(f"S{s}: {str(e)[:80]}")
    return sectors, errors
