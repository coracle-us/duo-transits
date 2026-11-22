#!/usr/bin/env python3
"""Round 3, section 11: ephemeris, odd/even, secondary eclipse and shape checks on four periodic stars.

For each star (STARS below, periods from round 2):

1. Ephemeris. Every predicted transit that is clean (box coverage >= core.CLEAN_FRAC and at least
   core.CLEAN_EDGE_D from a gap) gets a trapezoid fit (mid-time, depth, total duration, ingress, linear
   baseline) on a light curve detrended with all transits masked. Light curves are SPOC or TESS-SPOC where
   MAST has them, QLP otherwise before S94, and our own TESScut photometry on the target's barycentric
   clock (round2/s6_qlp_times.py) instead of QLP from S94 on. Mid-times whose error falls outside
   st2_pipeline.EVENT_TIME_ERR_MIN are degenerate fits and dropped, as in round 2. P and T0 come from a
   weighted least-squares line through the mid-times (errors scaled by sqrt(reduced chi2) when above 1).
   T0 is the paper event, so epoch numbers match round 2's k.
2. Odd/even. One trapezoid fitted to all clean transits together (shared duration and ingress), with
   separate odd and even depths. Fit errors are inflated for correlated noise (binned residuals) and by
   the excess scatter of single-transit depths about their odd and even means. The same data at 2P.
   Depths come from SPOC, TESS-SPOC and QLP flux at every sector, S94+ QLP included: our TESScut
   aperture photometry is diluted (its depths ran 10-40% under QLP's for the same events in
   results/round2/qlp_s94_timing.csv). QLP's S94+ clock offset (up to ~6 min) doesn't matter here,
   because each transit window is centred on a mid-time fitted in the same light curve.
3. Secondary eclipse. Transits masked (+/- one duration), folded at P, a box of the transit duration slid
   over phase 0.1-0.9 in steps of a tenth of the duration. Noise is core.duration_scale_noise of the
   masked light curve, divided by sqrt(total covered durations) in the box.
4. Shape and density. Ingress / total duration and Rp from the single-depth joint fit and TIC R*; stellar
   density from P, depth and the durations (Seager & Mallen-Ornelas 2003, circular orbit), Monte Carlo
   over the fit covariance, against the TIC density.
5. Catalogs, read-only: TESS EB catalog (Prsa et al. 2022, VizieR J/ApJS/258/16) and VSX (B/vsx/vsx)
   within 30 arcsec.

  python round3/s11_oddeven.py            # all four stars
  python round3/s11_oddeven.py 275647095  # one star
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "round2"))
import st2_pipeline as pl  # noqa: E402
import st2core as core  # noqa: E402
import tesscut_phot as tp  # noqa: E402
from s6_qlp_times import roemer_s  # noqa: E402

warnings.filterwarnings("ignore")
WORK, OUT = HERE / "work", HERE / "results" / "round3"
RAW_CACHE = WORK / "lc_cache" / "round3"
RJ_PER_RSUN = 6.957e8 / 7.1492e7
RHO_SUN = 1.41  # g cm^-3
QLP_TESSCUT_FROM = 94
STARS = [  # tic, cand_id, period (d), box depth (ppm), duration (h) from round 2
    (275647095, "275647095_2315.517", 17.0643, 9272, 7.57),
    (241550976, "241550976_2339.445", 21.1971, 16550, 4.61),
    (107113048, "107113048_2312.944", 49.7018, 6688, 6.69),
    (230133616, "230133616_2396.131", 23.7535, 5296, 8.19),
]
NSUB = 7  # supersamples per exposure
VIZIER = "vizier.cfa.harvard.edu"  # the CDS server timed out from here


# ---------------------------------------------------------------- light curves

def raw_segments(tic, ra, dec):
    """Normalised, quality-masked flux per sector before detrending: list of (sector, source, t, f)."""
    RAW_CACHE.mkdir(parents=True, exist_ok=True)
    path = RAW_CACHE / f"{tic}.npz"
    if path.exists():
        d = np.load(path, allow_pickle=True)
        return list(d["segs"])
    sr, prods = pl.lc_products(tic)
    best = pl.best_per_sector(prods)
    segs = []
    for s, row in sorted(best.items()):
        try:
            t, f = pl.load_sector_lc(sr, row, WORK / "lc_cache")
            segs.append((int(s), row["author"], t, f))
            print(f"  TIC {tic} S{s} {row['author']}: {t.size} points")
            if row["author"] == "QLP" and s >= QLP_TESSCUT_FROM:
                t, f = tesscut_lightcurve(tic, ra, dec, s)
                segs.append((int(s), "TESScut", t, f))
                print(f"  TIC {tic} S{s} TESScut: {t.size} points")
        except Exception as e:
            print(f"  TIC {tic} S{s}: {e}")
    np.savez_compressed(path, segs=np.array(segs, dtype=object))
    return segs


def tesscut_lightcurve(tic, ra, dec, sector):
    """TESScut aperture photometry on the target's own barycentric clock (as round2/s6_qlp_times.py)."""
    tpf = tp.cutout(ra, dec, sector, WORK / "tesscut")
    d = tpf.hdu[1].data
    good = np.isfinite(d["TIME"]) & np.isfinite(d["TIMECORR"]) & (d["TIME"] > 0)
    t_sc = d["TIME"] - d["TIMECORR"]
    t_true = t_sc[good] + roemer_s(t_sc[good] + 2457000, ra, dec) / 86400
    t, f, _, _ = tp.sector_lightcurve(tic, ra, dec, sector, WORK / "tesscut")
    t = np.interp(t, d["TIME"][good], t_true)
    cad = core.cadence_days(t)
    if cad < 9 / 1440:  # same 10-min binning as pl.load_sector_lc
        edges = np.arange(t.min(), t.max() + 10 / 1440, 10 / 1440)
        idx = np.digitize(t, edges)
        n = np.bincount(idx)
        ok = n > 0
        t, f = (np.bincount(idx, t)[ok] / n[ok]), (np.bincount(idx, f)[ok] / n[ok])
    return core.drop_common_mode(t, f, sector, "TESScut")


def transit_mask(t, P, T0, half_d):
    ph = (t - T0 + P / 2) % P - P / 2
    return np.abs(ph) < half_d


def timing_set(segs):
    return [g for g in segs if not (g[1] == "QLP" and g[0] >= QLP_TESSCUT_FROM)]


def phot_set(segs):
    return [g for g in segs if g[1] != "TESScut"]


def detrend_masked(segs, dur_d, P, T0):
    """Biweight detrend (window as the pipeline) with transits masked; trend interpolated across them.

    Returns concatenated t, f, exposure (d), source, sector arrays, sorted by time.
    """
    from wotan import flatten
    w = max(4.0 * dur_d, 1.0)
    T, F, X, S, N = [], [], [], [], []
    for sector, src, t, f in segs:
        m = transit_mask(t, P, T0, 0.75 * dur_d)
        if (~m).sum() < 50:
            continue
        _, trend = flatten(t[~m], f[~m], method="biweight", window_length=w, edge_cutoff=0.0,
                           break_tolerance=0.5, return_trend=True)
        ok = np.isfinite(trend)
        tr = np.interp(t, t[~m][ok], trend[ok])
        flat = f / tr
        # drop masked points that sit beyond the ends of the trend or across a gap
        ts = t[~m][ok]
        j = np.clip(np.searchsorted(ts, t), 1, ts.size - 1)
        bridge = np.minimum(np.abs(t - ts[j - 1]), np.abs(ts[j] - t)) > 0.75 * dur_d + 0.1
        keep = np.isfinite(flat) & ~bridge
        r = flat[keep & ~m] - 1
        s = core.robust_std(r)
        keep &= flat - 1 < 5 * s  # upward outliers only
        T.append(t[keep]); F.append(flat[keep])
        X.append(np.full(keep.sum(), core.cadence_days(t)))
        S.append(np.full(keep.sum(), src, dtype=object)); N.append(np.full(keep.sum(), sector))
    t, f, x, s, n = (np.concatenate(a) for a in (T, F, X, S, N))
    o = np.argsort(t)
    return t[o], f[o], x[o], s[o], n[o]


# ---------------------------------------------------------------- models

def trap(x, depth, t14, t12):
    t12 = np.clip(t12, 1e-4, t14 / 2)
    a = np.abs(x)
    m = np.ones_like(x)
    full = a <= t14 / 2 - t12
    ing = (a > t14 / 2 - t12) & (a < t14 / 2)
    m[full] -= depth
    m[ing] -= depth * (t14 / 2 - a[ing]) / t12
    return m


def trap_ss(x, exp, depth, t14, t12):
    """Trapezoid averaged over each point's exposure (depth may be an array)."""
    off = (np.arange(NSUB) + 0.5) / NSUB - 0.5
    xs = x[:, None] + exp[:, None] * off[None, :]
    d = np.broadcast_to(depth, x.shape)[:, None]
    return (1 - d * (1 - trap(xs.ravel(), 1.0, t14, t12).reshape(xs.shape))).mean(1)


def beta_factor(res, t, scale_d):
    """Red-noise inflation: rms of residuals binned at scale_d vs the white expectation."""
    o = np.argsort(t)
    res = res[o]
    nb = max(int(round(scale_d / core.cadence_days(t))), 1)
    if nb < 2 or res.size < 4 * nb:
        return 1.0
    k = res.size // nb
    b = res[: k * nb].reshape(k, nb).mean(1)
    return max(1.0, float(np.std(b) / (np.std(res) / np.sqrt(nb))))


def cov_from(r, npar):
    s2 = np.sum(r.fun ** 2) / max(r.fun.size - npar, 1)
    return np.linalg.pinv(r.jac.T @ r.jac) * s2


def fit_one(t, f, exp, tc0, dur0, depth0):
    """Trapezoid + linear baseline on one transit window. Returns dict with tc, tc_err (d), depth, t14, t12."""
    w = np.abs(t - tc0) < max(2.0 * dur0, 0.5)
    t, f, exp = t[w], f[w], exp[w]
    out = dict(tc=np.nan, tc_err=np.nan, depth=np.nan, depth_err=np.nan, t14=np.nan, t12=np.nan, n=int(t.size))
    if t.size < 20:
        return out

    def resid(p):
        tc, d, t14, t12, a, b = p
        return trap_ss(t - tc, exp, d, t14, t12) * (a + b * (t - tc0)) - f

    p0 = [tc0, depth0, dur0, dur0 / 8, 1.0, 0.0]
    lo = [tc0 - dur0 / 2, 0, dur0 / 3, 1e-3, 0.9, -1]
    hi = [tc0 + dur0 / 2, 5 * depth0, 3 * dur0, dur0, 1.1, 1]
    r = least_squares(resid, p0, bounds=(lo, hi), x_scale=[dur0 / 10, depth0, dur0, dur0 / 8, 1e-3, 1e-3])
    cov = cov_from(r, 6)
    beta = beta_factor(r.fun, t, max(r.x[3], 15 / 1440))
    out.update(tc=r.x[0], tc_err=np.sqrt(cov[0, 0]) * beta, depth=r.x[1], depth_err=np.sqrt(cov[1, 1]) * beta,
               t14=r.x[2], t12=min(r.x[3], r.x[2] / 2), beta=beta)
    return out


def joint_fit(x, f, exp, groups, dur0, depth0, tabs):
    """Shared offset, t14, t12 and baseline; one depth per group (0..G-1). tabs (absolute times) orders
    the residuals for the red-noise factor. Returns params, cov (red-noise inflated), result."""
    G = int(groups.max()) + 1

    def resid(p):
        dx, t14, t12, c = p[:4]
        d = p[4:][groups]
        return trap_ss(x - dx, exp, d, t14, t12) * (1 + c) - f

    p0 = [0.0, dur0, dur0 / 8, 0.0] + [depth0] * G
    lo = [-dur0 / 4, dur0 / 3, 1e-3, -0.01] + [0] * G
    hi = [dur0 / 4, 3 * dur0, 1.5 * dur0, 0.01] + [5 * depth0] * G
    r = least_squares(resid, p0, bounds=(lo, hi), x_scale=[dur0 / 10, dur0, dur0 / 8, 1e-4] + [depth0] * G)
    cov = cov_from(r, len(p0)) * beta_factor(r.fun, tabs, max(r.x[2], 15 / 1440)) ** 2
    return r.x, cov, r


# ---------------------------------------------------------------- per-star analysis

def ephemeris(t, f, exp, src, sec, P, T0, dur_d, depth):
    edges = core.gap_edges(t)
    k = np.arange(int(np.floor((t.min() - T0) / P)) - 1, int(np.ceil((t.max() - T0) / P)) + 2)
    ep = T0 + k * P
    _, frac, _ = core.box_stats(t, f, ep, dur_d)
    clean = (frac >= core.CLEAN_FRAC) & (core.edge_distance(edges, ep) >= core.CLEAN_EDGE_D)
    rows = []
    for kk, e, fr, cl in zip(k, ep, frac, clean):
        if fr <= 0:
            continue
        near = np.abs(t - e) < dur_d
        row = dict(epoch=int(kk), t_pred=e, coverage=fr, edge_d=float(core.edge_distance(edges, e)[0]), clean=bool(cl),
                   sector=int(np.median(sec[near])) if near.any() else -1,
                   source=pd.Series(src[near]).mode().iloc[0] if near.any() else "")
        if cl:
            row.update(fit_one(t, f, exp, e, dur_d, depth))
            lo, hi = pl.EVENT_TIME_ERR_MIN
            row["degenerate"] = not (lo <= row["tc_err"] * 1440 <= hi)
        rows.append(row)
    tab = pd.DataFrame(rows)
    use = tab["clean"] & ~tab.get("degenerate", pd.Series(False, index=tab.index)).fillna(True).astype(bool)
    use &= np.isfinite(tab.get("tc", np.nan))
    n, tc, err = tab.loc[use, "epoch"].to_numpy(float), tab.loc[use, "tc"].to_numpy(), tab.loc[use, "tc_err"].to_numpy()
    A = np.column_stack([np.ones_like(n), n]) / err[:, None]
    sol, *_ = np.linalg.lstsq(A, tc / err, rcond=None)
    cov = np.linalg.inv(A.T @ A)
    resid = tc - (sol[0] + sol[1] * n)
    chi2 = float(np.sum((resid / err) ** 2))
    dof = max(n.size - 2, 1)
    scale = max(1.0, np.sqrt(chi2 / dof))
    T0f, Pf = sol
    tab["used"] = use
    tab["oc_min"] = (tab["tc"] - (T0f + tab["epoch"] * Pf)) * 1440
    tab["oc_err_min"] = tab["tc_err"] * 1440
    tab["oc_sigma"] = tab["oc_min"] / tab["oc_err_min"]
    tab["off_3sigma"] = use & (tab["oc_sigma"].abs() > 3)
    eph = dict(P=Pf, P_err=np.sqrt(cov[1, 1]) * scale, T0=T0f, T0_err=np.sqrt(cov[0, 0]) * scale,
               chi2_red=chi2 / dof, n_times=int(n.size), n_off_3sigma=int(tab["off_3sigma"].sum()))
    return eph, tab


def windows(t, f, exp, src, centres, epochs, dur_d):
    """Points within 2 durations of each centre: x (time from that centre), f, exp, src, epoch, t."""
    X, Fs, E, S, K, T = [], [], [], [], [], []
    for c, k in zip(centres, epochs):
        m = np.abs(t - c) < max(2 * dur_d, 0.5)
        X.append(t[m] - c); Fs.append(f[m]); E.append(exp[m]); S.append(src[m]); K.append(np.full(m.sum(), k))
        T.append(t[m])
    return [np.concatenate(a) for a in (X, Fs, E, S, K, T)]


def phot_transits(t, f, exp, src, P, T0, dur_d, depth):
    """Clean transits in the photometry set, each with its own fitted mid-time (same light curve)."""
    edges = core.gap_edges(t)
    k = np.arange(int(np.floor((t.min() - T0) / P)) - 1, int(np.ceil((t.max() - T0) / P)) + 2)
    ep = T0 + k * P
    _, frac, _ = core.box_stats(t, f, ep, dur_d)
    clean = (frac >= core.CLEAN_FRAC) & (core.edge_distance(edges, ep) >= core.CLEAN_EDGE_D)
    rows = []
    for kk, e in zip(k[clean], ep[clean]):
        r = fit_one(t, f, exp, e, dur_d, depth)
        lo, hi = pl.EVENT_TIME_ERR_MIN
        near = np.abs(t - e) < dur_d
        rows.append(dict(epoch=int(kk), phot_source=pd.Series(src[near]).mode().iloc[0], phot_tc=r["tc"],
                         phot_tc_err_min=r["tc_err"] * 1440, phot_ok=bool(lo <= r["tc_err"] * 1440 <= hi)))
    return pd.DataFrame(rows)


def single_depths(x, f, exp, kw, tabs, shape, dur0, depth0):
    """Depth of each transit with the joint shape (offset, t14, t12) fixed; error inflated for red noise."""
    dx, t14, t12 = shape
    out = {}
    for k in np.unique(kw):
        m = kw == k

        def resid(p):
            return trap_ss(x[m] - dx, exp[m], p[0], t14, t12) * (1 + p[1]) - f[m]
        r = least_squares(resid, [depth0, 0.0], x_scale=[depth0, 1e-4])
        cov = cov_from(r, 2)
        out[int(k)] = (r.x[0], np.sqrt(cov[0, 0]) * beta_factor(r.fun, tabs[m], max(t12, 15 / 1440)))
    return out


def density_mc(P, pars, cov, idx, rng, n=4000):
    """Stellar density (g/cc) from P and (depth, t14, t12) by Seager & Mallen-Ornelas, circular orbit."""
    sub = np.asarray(idx)
    draws = rng.multivariate_normal(pars[sub], cov[np.ix_(sub, sub)], size=n)
    out = []
    for d, t14, t12 in np.column_stack([draws[:, 0], draws[:, 1], draws[:, 2]]):
        if d <= 0 or t14 <= 0 or t12 <= 0:
            continue
        tF = t14 - 2 * t12
        k = np.sqrt(d)
        sT, sF = np.sin(np.pi * t14 / P), np.sin(np.pi * max(tF, 0) / P)
        q = (sF / sT) ** 2
        b2 = ((1 - k) ** 2 - q * (1 + k) ** 2) / (1 - q) if q < 1 else np.nan
        b2 = max(b2, 0.0) if np.isfinite(b2) else np.nan
        if not np.isfinite(b2):
            continue
        aR = np.sqrt(((1 + k) ** 2 - b2 * (1 - sT ** 2)) / sT ** 2)
        rho = 3 * np.pi / (6.674e-8 * (P * 86400) ** 2) * aR ** 3
        out.append(rho)
    out = np.array(out)
    if out.size < 100:
        return np.nan, np.nan, np.nan
    return tuple(np.percentile(out, [50, 16, 84]))


def secondary_scan(t, f, exp, P, T0, dur_d):
    m = ~transit_mask(t, P, T0, dur_d)
    t, f = t[m], f[m]
    lc = core.local_cadence(t)
    sig = core.duration_scale_noise(t, f, dur_d)
    ph = ((t - T0) / P) % 1.0
    o = np.argsort(ph)
    ph, f, lc = ph[o], f[o], lc[o]
    cs, cc = np.concatenate([[0], np.cumsum(f)]), np.concatenate([[0], np.cumsum(lc)])
    half = dur_d / P / 2
    cen = np.arange(0.1, 0.9 + 1e-9, dur_d / P / 10)
    lo, hi = np.searchsorted(ph, cen - half), np.searchsorted(ph, cen + half)
    nn = hi - lo
    with np.errstate(invalid="ignore", divide="ignore"):
        depth = 1 - (cs[hi] - cs[lo]) / nn
        cover = (cc[hi] - cc[lo]) / dur_d  # covered durations in this box, summed over orbits
        err = sig / np.sqrt(cover)
    ok = cover >= 0.5
    scan = pd.DataFrame(dict(phase=cen, depth=depth, err=err, snr=depth / err, cover=cover))[ok]
    return scan, (t, f, ph, sig)


def drop_top_orbit(t, f, sec, P, T0, dur_d, phase, sig):
    """The box at this phase with its largest single-orbit contribution removed: is one orbit the dip?"""
    m = ~transit_mask(t, P, T0, dur_d)
    t, f, sec = t[m], f[m], sec[m]
    lc = core.local_cadence(t)
    n = np.round((t - T0) / P - phase)
    w = np.abs(t - (T0 + (n + phase) * P)) < dur_d / 2
    df = pd.DataFrame(dict(orbit=n[w], sec=sec[w], r=1 - f[w], c=lc[w]))
    g = df.groupby("orbit").agg(sum_r=("r", lambda v: v.sum()), n=("r", "size"), c=("c", "sum"), sec=("sec", "first"))
    g["depth"] = g["sum_r"] / g["n"]
    top = g["sum_r"].idxmax()
    rest = g.drop(top)
    d = rest["sum_r"].sum() / rest["n"].sum()
    err = sig / np.sqrt(rest["c"].sum() / dur_d)
    return dict(top_orbit_sector=int(g.loc[top, "sec"]), top_orbit_depth_ppm=g.loc[top, "depth"] * 1e6,
                without_top_depth_ppm=d * 1e6, without_top_snr=d / err)


def catalog_matches(ra, dec, tic):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier
    v = Vizier(columns=["**", "+_r"], row_limit=50, vizier_server=VIZIER)
    c = SkyCoord(ra * u.deg, dec * u.deg)
    out = dict(tess_eb="none", vsx="none")
    try:
        r = v.query_region(c, radius=30 * u.arcsec, catalog="J/ApJS/258/16")
        hits = []
        for tab in r:
            cols = tab.colnames
            for row in tab:
                tid = row["TIC"] if "TIC" in cols else "?"
                per = row["Per"] if "Per" in cols else (row["period"] if "period" in cols else "?")
                hits.append(f"TIC {tid} P={per} d r={float(row['_r']) * 60:.1f}\"")
        out["tess_eb"] = "; ".join(sorted(set(hits))) if hits else "none within 30\""
        out["tess_eb_target"] = any(f"TIC {tic} " in h for h in hits)
    except Exception as e:
        out["tess_eb"] = f"query failed: {e}"[:200]
    try:
        r = v.query_region(c, radius=30 * u.arcsec, catalog="B/vsx/vsx")
        hits = []
        for tab in r:
            for row in tab:
                hits.append(f"{row['Name']} {row['Type']} P={row['Period']} d r={float(row['_r']) * 60:.1f}\"")
        out["vsx"] = "; ".join(hits) if hits else "none within 30\""
    except Exception as e:
        out["vsx"] = f"query failed: {e}"[:200]
    return out


def analyse(tic, cand_id, P0, depth_ppm, dur_h, cands, cov_tab):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(11)
    info = cov_tab.loc[tic]
    T0 = float(cands.loc[cands["cand_id"] == cand_id, "t0"].iloc[0])
    dur_d, depth = dur_h / 24, depth_ppm * 1e-6
    segs = raw_segments(tic, info.ra, info.dec)
    P = P0
    for _ in range(2):  # detrend with the latest ephemeris, refit, repeat once
        t, f, exp, src, sec = detrend_masked(timing_set(segs), dur_d, P, T0)
        eph, tab = ephemeris(t, f, exp, src, sec, P, T0, dur_d, depth)
        P, T0 = eph["P"], eph["T0"]

    # photometry set: every transit centred on its own mid-time in the same light curve
    tp_, fp, ep_, sp, secp = detrend_masked(phot_set(segs), dur_d, P, T0)
    ph = phot_transits(tp_, fp, ep_, sp, P, T0, dur_d, depth)
    ph = ph[ph["phot_ok"]]
    x, fw, ew, sw, kw, tw = windows(tp_, fp, ep_, sp, ph["phot_tc"].to_numpy(), ph["epoch"].to_numpy(), dur_d)
    epochs = ph["epoch"].to_numpy(int)
    p1, c1, _ = joint_fit(x, fw, ew, np.zeros(x.size, int), dur_d, depth, tw)  # single depth
    par = (kw % 2).astype(int)  # 0 even, 1 odd
    p2, c2, _ = joint_fit(x, fw, ew, par, dur_d, depth, tw)
    d_even, d_odd = p2[4], p2[5]
    diff = d_odd - d_even
    diff_err_fit = float(np.sqrt(c2[4, 4] + c2[5, 5] - 2 * c2[4, 5]))
    # single-transit depths with the joint shape: their scatter about the odd and even means
    sd = single_depths(x, fw, ew, kw, tw, p1[:3], dur_d, depth)
    ph["depth_single_ppm"] = [sd[k][0] * 1e6 for k in ph["epoch"]]
    ph["depth_single_err_ppm"] = [sd[k][1] * 1e6 for k in ph["epoch"]]
    ph["parity"] = np.where(ph["epoch"] % 2 == 1, "odd", "even")
    chi2 = 0.0
    for g, gg in ph.groupby("parity"):
        d, e = gg["depth_single_ppm"].to_numpy(), gg["depth_single_err_ppm"].to_numpy()
        mu = np.sum(d / e ** 2) / np.sum(1 / e ** 2)
        chi2 += float(np.sum(((d - mu) / e) ** 2))
    dof = len(ph) - ph["parity"].nunique()
    scat = max(1.0, np.sqrt(chi2 / dof)) if dof > 0 else np.nan
    diff_err = diff_err_fit * scat if np.isfinite(scat) else diff_err_fit
    sigma = diff / diff_err
    ev_e, ev_o = epochs[epochs % 2 == 0], epochs[epochs % 2 == 1]
    # odd / even within the most common reduction only
    main_src = pd.Series(sw).mode().iloc[0]
    m1 = sw == main_src
    oe1 = dict(diff=np.nan, err=np.nan)
    n1e, n1o = np.unique(kw[m1 & (par == 0)]).size, np.unique(kw[m1 & (par == 1)]).size
    if n1e and n1o:
        p3, c3, _ = joint_fit(x[m1], fw[m1], ew[m1], par[m1], dur_d, depth, tw[m1])
        oe1 = dict(diff=p3[5] - p3[4], err=float(np.sqrt(c3[4, 4] + c3[5, 5] - 2 * c3[4, 5])) * (scat if np.isfinite(scat) else 1))
    verdict_2p = (f"two equal eclipses at 2P (odd - even {sigma:+.1f} sigma)" if abs(sigma) < 3 else
                  f"one deeper, one shallower eclipse at 2P (odd - even {sigma:+.1f} sigma)")

    # shape and density
    dpt, t14, t12 = p1[4], p1[1], p1[2]
    dpt_err, t14_err, t12_err = np.sqrt(np.diag(c1))[[4, 1, 2]]
    rp = np.sqrt(dpt) * info.rstar * RJ_PER_RSUN
    rho_tr, rho_lo, rho_hi = density_mc(P, p1, c1, [4, 1, 2], rng)
    rho_tic = info.mstar / info.rstar ** 3 * RHO_SUN

    # secondary, on the photometry set
    scan, (ts, fs, phs, sig) = secondary_scan(tp_, fp, ep_, P, T0, dur_d)
    best = scan.loc[scan["snr"].idxmax()]
    deepest = scan.loc[scan["depth"].idxmax()]
    half = scan.iloc[(scan["phase"] - 0.5).abs().argmin()]
    jk = drop_top_orbit(tp_, fp, secp, P, T0, dur_d, best.phase, sig)
    cat = catalog_matches(info.ra, info.dec, tic)

    tab = tab.merge(ph, on="epoch", how="left")
    tab.insert(0, "tic", tic)
    tab.to_csv(OUT / f"tic{tic}_transits.csv", index=False)
    used = tab[tab["used"]]
    srcs = used["source"].value_counts().to_dict()
    psrcs = ph["phot_source"].value_counts().to_dict()
    row = dict(tic=tic, cand_id=cand_id, P_d=P, P_err_d=eph["P_err"], T0_btjd=T0, T0_err_d=eph["T0_err"],
               n_transits_timed=eph["n_times"], timing_chi2_red=eph["chi2_red"], n_oc_over_3sigma=eph["n_off_3sigma"],
               sources_timed="; ".join(f"{k} {v}" for k, v in srcs.items()),
               sources_depth="; ".join(f"{k} {v}" for k, v in psrcs.items()),
               n_even=int(ev_e.size), n_odd=int(ev_o.size),
               depth_even_ppm=d_even * 1e6, depth_odd_ppm=d_odd * 1e6, depth_odd_minus_even_ppm=diff * 1e6,
               depth_diff_err_ppm=diff_err * 1e6, depth_diff_err_fit_ppm=diff_err_fit * 1e6,
               single_depth_scatter_factor=scat, odd_even_sigma=sigma,
               odd_even_one_reduction=f"{main_src}: {n1e} even, {n1o} odd",
               odd_even_one_reduction_ppm=oe1["diff"] * 1e6, odd_even_one_reduction_err_ppm=oe1["err"] * 1e6,
               verdict_2p=verdict_2p,
               sec_best_phase=best.phase, sec_best_depth_ppm=best.depth * 1e6, sec_best_snr=best.snr,
               sec_best_top_orbit_sector=jk["top_orbit_sector"],
               sec_best_top_orbit_depth_ppm=jk["top_orbit_depth_ppm"],
               sec_best_without_top_orbit_depth_ppm=jk["without_top_depth_ppm"],
               sec_best_without_top_orbit_snr=jk["without_top_snr"],
               sec_deepest_phase=deepest.phase, sec_deepest_depth_ppm=deepest.depth * 1e6, sec_deepest_snr=deepest.snr,
               sec_phase05_depth_ppm=half.depth * 1e6, sec_phase05_err_ppm=half.err * 1e6,
               sec_phase05_upper3s_ppm=(half.depth + 3 * half.err) * 1e6,
               depth_ppm=dpt * 1e6, depth_err_ppm=dpt_err * 1e6, t14_h=t14 * 24, t14_err_h=t14_err * 24,
               t12_h=t12 * 24, t12_err_h=t12_err * 24, ingress_ratio=t12 / t14,
               rstar_tic=info.rstar, mstar_tic=info.mstar, contratio_tic=info.contratio, rp_rj=rp,
               rho_transit_gcc=rho_tr, rho_transit_lo=rho_lo, rho_transit_hi=rho_hi, rho_tic_gcc=rho_tic,
               density_ratio=rho_tr / rho_tic, tess_eb=cat["tess_eb"], vsx=cat["vsx"])

    # plots
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    bw = dur_d / 12
    xm = np.linspace(x.min(), x.max(), 600)
    ylim = (1 - 1.8 * depth, 1 + 0.8 * depth)
    for g, lab, col in ((0, f"even ({ev_e.size})", "C0"), (1, f"odd ({ev_o.size})", "C3")):
        m = par == g
        xb, fb = binned(x[m], fw[m], bw)
        ax[0, 0].plot(x[m] * 24, fw[m], ",", color=col, alpha=0.3)
        ax[0, 0].plot(xb * 24, fb, "o", ms=3, color=col, label=lab)
        ax[0, 0].plot(xm * 24, trap(xm - p2[0], p2[4 + g], p2[1], p2[2]) * (1 + p2[3]), color=col, lw=1)
    ax[0, 0].set(xlabel="hours from mid-transit", ylabel="flux", ylim=ylim,
                 title=f"odd - even = {diff * 1e6:+.0f} +/- {diff_err * 1e6:.0f} ppm ({sigma:+.1f} sigma)")
    ax[0, 0].legend(fontsize=8)
    ph2 = x / (2 * P) + 0.5 * par  # odd transits sit at phase 0.5 of 2P
    for g, col in ((0, "C0"), (1, "C3")):
        m = par == g
        xb, fb = binned(ph2[m], fw[m], bw / (2 * P))
        ax[0, 1].plot(ph2[m], fw[m], ",", color="0.6")
        ax[0, 1].plot(xb, fb, "o", ms=2, color=col)
    ax[0, 1].set(xlabel=f"phase at 2P = {2 * P:.4f} d (even at 0, odd at 0.5)", ylabel="flux", ylim=ylim,
                 title=verdict_2p)
    ax[0, 1].set_xlim(-0.1, 0.6)
    phall = ((tp_ - T0) / P + 0.25) % 1 - 0.25
    xb, fb = binned(phall, fp, dur_d / P / 3)
    ax[1, 0].plot(xb, (fb - 1) * 1e6, ".", ms=3, color="k", label="all data, binned")
    ax[1, 0].plot(scan.phase, -scan.depth * 1e6, color="C1", lw=1, label="box scan (transits masked)")
    ax[1, 0].axvspan(-dur_d / P, dur_d / P, color="C3", alpha=0.1)
    lim = max(8 * scan["err"].median() * 1e6, 300)
    ax[1, 0].set(xlabel="phase", ylabel="flux - 1 (ppm)", ylim=(-lim, lim / 2),
                 title=f"best dip phase {best.phase:.3f}: {best.depth * 1e6:.0f} ppm, SNR {best.snr:.1f}; "
                       f"phase 0.5 < {(half.depth + 3 * half.err) * 1e6:.0f} ppm (3 sigma)")
    ax[1, 0].legend(fontsize=8, loc="lower right")
    u = tab[tab["used"]]
    for s_, g_ in u.groupby("source"):
        ax[1, 1].errorbar(g_["epoch"], g_["oc_min"], g_["oc_err_min"], fmt="o", ms=4, label=s_)
    ax[1, 1].axhline(0, color="k", lw=0.8)
    ax[1, 1].set(xlabel="epoch (0 = paper event)", ylabel="O - C (min)",
                 title=f"P = {P:.5f} +/- {eph['P_err']:.5f} d, chi2_red {eph['chi2_red']:.2f}")
    ax[1, 1].legend(fontsize=8)
    fig.suptitle(f"TIC {tic} ({cand_id})")
    fig.tight_layout()
    fig.savefig(OUT / f"tic{tic}_oddeven.png", dpi=110)
    plt.close(fig)
    scan.assign(tic=tic).to_csv(OUT / f"tic{tic}_secondary_scan.csv", index=False)
    return row


def binned(x, y, w):
    if x.size == 0:
        return x, y
    e = np.arange(x.min(), x.max() + w, w)
    i = np.digitize(x, e)
    n = np.bincount(i)
    ok = n >= 3
    return (np.bincount(i, x)[ok] / n[ok]), (np.bincount(i, y)[ok] / n[ok])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pick = {int(a) for a in sys.argv[1:]}
    cands = pd.read_csv(HERE / "candidates.csv")
    cands["cand_id"] = [f"{t}_{x:.3f}" for t, x in zip(cands["tic"], cands["t0"])]
    cov_tab = pd.read_csv(WORK / "coverage.csv").set_index("tic")
    rows = []
    for s in STARS:
        if pick and s[0] not in pick:
            continue
        print(f"TIC {s[0]}", flush=True)
        r = analyse(*s, cands, cov_tab)
        print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
        rows.append(r)
    path = OUT / "oddeven_summary.csv"
    new = pd.DataFrame(rows)
    if pick and path.exists():
        old = pd.read_csv(path)
        new = pd.concat([old[~old["tic"].isin(new["tic"])], new]).sort_values("tic", key=lambda c: c.map(
            {s[0]: i for i, s in enumerate(STARS)}))
    new.to_csv(path, index=False)


if __name__ == "__main__":
    main()
