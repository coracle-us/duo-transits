"""
Eclipse-timing library for the EB timing side quest.

Light curves: scripts/alias_pipeline.load_lightcurves (imported, pointed at our cache dir),
one product per sector (SPOC 120 s > TESS-SPOC > QLP), quality == 0, sector-median normalised.

Template: the whole phase-folded light curve (not just the eclipses), binned finely in phase
(median per bin) from the highest-cadence data, so the ellipsoidal / reflection shape of
close binaries is part of the model. For a long-cadence (FFI) point the template is averaged
over the exposure (sub-sampled box), so 200/600/1800-s data are fitted with a matching
template instead of a 2-min one.

Per eclipse: data within +-h of the predicted minimum; model
    f(t) = a + b*(t-T) + c*(tmpl(t - dt) - 1)
(c rescales the eclipse+variability amplitude, absorbing dilution differences between
products/sectors). For each trial dt on a fine grid (a, b, c) are solved linearly, giving a
profile chi^2(dt); the minimum is refined by a parabola and sigma_dt comes from the profile
curvature (Delta chi^2 = 1) with chi^2 scaled to reduced chi^2 = 1 (white-noise estimate).
Requirements per eclipse: >= MIN_FRAC of the expected cadences on EACH side of minimum,
and the minimum not at the edge of the search grid.

BTJD throughout.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import alias_pipeline as ap  # noqa: E402

ap.LC = HERE / "cache" / "lc"          # read our cache, same file format as data/raw
PRODUCTS = ["SPOC_120", "TESS-SPOC", "QLP"]
MIN_FRAC = 0.6


def load(tic):
    lc = ap.load_lightcurves(tic, products=PRODUCTS)
    if not len(lc):
        return lc
    # exposure per sector = median cadence (FFI exposure ~ cadence)
    cad = lc.groupby("sector").time.apply(lambda t: np.median(np.diff(t.values)))
    lc["exp"] = lc.sector.map(cad).values
    # per-segment robust re-normalisation (sector median already applied)
    return lc


# ------------------------------------------------------------------ template
class Template:
    """Periodic phase template; phase 0 = primary minimum."""

    def __init__(self, phase, flux, nbins):
        edges = np.linspace(0, 1, nbins + 1)
        idx = np.clip(np.digitize(phase % 1, edges) - 1, 0, nbins - 1)
        med = np.full(nbins, np.nan)
        cnt = np.bincount(idx, minlength=nbins)
        order = np.argsort(idx)
        splits = np.split(flux[order], np.cumsum(cnt)[:-1])
        for i, s in enumerate(splits):
            if len(s) >= 3:
                med[i] = np.median(s)
        x = 0.5 * (edges[1:] + edges[:-1])
        ok = np.isfinite(med)
        med = np.interp(x, x[ok], med[ok], period=1)
        # light 3-bin running mean (periodic) to suppress bin noise
        med = (np.roll(med, 1) + med + np.roll(med, -1)) / 3
        self.x, self.y, self.nbins, self.count = x, med, nbins, cnt

    def __call__(self, phase):
        return np.interp(phase % 1, self.x, self.y, period=1)

    def smeared(self, phase, exp_phase):
        """Average over an exposure of exp_phase (phase units) centred on phase."""
        exp_phase = np.atleast_1d(exp_phase)
        nsub = int(np.clip(np.ceil(np.max(exp_phase) * self.nbins * 2), 1, 41))
        if nsub == 1:
            return self(phase)
        off = (np.arange(nsub) + 0.5) / nsub - 0.5
        return np.mean([self(phase + o * exp_phase) for o in off], axis=0)


def build_template(t, f, P, T0, nbins=None):
    if nbins is None:
        nbins = int(np.clip(P * 86400 / 90, 400, 6000))   # ~90 s bins, 400..6000
    return Template((t - T0) / P, f, nbins)


def find_minima(tmpl):
    """Phase of primary (deepest) and secondary minimum (deepest >= 0.15 phase away)."""
    y, x = tmpl.y, tmpl.x
    i1 = int(np.argmin(y))
    d = np.abs(((x - x[i1]) + 0.5) % 1 - 0.5)
    y2 = np.where(d >= 0.15, y, np.inf)
    i2 = int(np.argmin(y2))
    # local-minimum sanity: secondary must be a local minimum of the template
    is_local = y[i2] <= y[(i2 - 3) % len(y)] and y[i2] <= y[(i2 + 3) % len(y)]
    return x[i1], (x[i2] if is_local else np.nan), y[i1], (y[i2] if is_local else np.nan)


def eclipse_halfwidth(tmpl, ph0):
    """Half-width (phase) of the eclipse at ph0: distance to where the template recovers
    90% of the way from the minimum to the local maximum within +-0.25 phase."""
    x = ((tmpl.x - ph0) + 0.5) % 1 - 0.5
    o = np.argsort(x)
    x, y = x[o], tmpl.y[o]
    ymin = np.interp(0, x, y)
    hw = []
    for side in (1, -1):
        m = (side * x > 0) & (side * x < 0.25)
        xs, ys = np.abs(x[m]), y[m]
        oo = np.argsort(xs)
        xs, ys = xs[oo], ys[oo]
        ymax = np.max(ys)
        lev = ymin + 0.9 * (ymax - ymin)
        k = np.argmax(ys >= lev)
        hw.append(xs[k])
    return float(np.mean(hw))


# ------------------------------------------------------------------ one eclipse
def fit_eclipse(t, f, exp, tmpl, P, T0ph, Tpred, h, search, step):
    """Profile-chi^2 timing of one eclipse. Phase of model = (t - dt - T0ph)/P, where T0ph is
    the epoch of phase 0 of the template; Tpred the predicted minimum. Returns dict or None."""
    x = t - Tpred
    n = len(t)

    def profile(grid):
        chi = np.empty(len(grid))
        sols = []
        for j, dt in enumerate(grid):
            ph = (t - dt - T0ph) / P
            m = tmpl.smeared(ph, exp / P) - 1.0
            A = np.column_stack([np.ones(n), x, m])
            coef, *_ = np.linalg.lstsq(A, f, rcond=None)
            r = f - A @ coef
            chi[j] = r @ r
            sols.append(coef)
        return chi, sols

    coarse = search / 25
    g0 = np.arange(-search, search + coarse / 2, coarse)
    c0, _ = profile(g0)
    j0 = int(np.argmin(c0))
    if j0 == 0 or j0 == len(g0) - 1:
        return None
    step = min(step, coarse / 10)
    grid = np.arange(g0[j0] - 1.5 * coarse, g0[j0] + 1.5 * coarse + step / 2, step)
    chi, sols = profile(grid)
    j = int(np.argmin(chi))
    if j == 0 or j == len(grid) - 1:
        return None
    # parabola through 3 points around minimum
    y0, y1, y2 = chi[j - 1], chi[j], chi[j + 1]
    den = y0 - 2 * y1 + y2
    if den <= 0:
        return None
    off = 0.5 * (y0 - y2) / den
    dt_best = grid[j] + off * step
    curv = den / step ** 2                       # d2chi/ddt2 (unscaled chi = sum r^2)
    dof = max(n - 4, 1)
    s2 = y1 / dof                                # white-noise variance estimate
    sig = np.sqrt(2 * s2 / curv)                 # Delta chi^2(scaled) = 1
    c = sols[j][2]
    return dict(dt=dt_best, sig=sig, amp=c, rms=np.sqrt(s2), n=n)


def time_eclipses(lc, tmpl, P, T0, ph_ecl, hw, search, kind, min_frac=MIN_FRAC, shift_fn=None):
    """Time every eclipse of one type. ph_ecl: template phase of this eclipse; hw: its
    half-width (phase). Window h = 1.5 hw (capped at 0.3 P). shift_fn(t) (days), if given, is
    a smooth prior O-C added to the linear-ephemeris prediction (and to the template
    phase), so the search can stay narrow when the ETV is large. Returns DataFrame; T is
    the measured mid-time, independent of the shift used."""
    t_all, f_all = lc.time.values, lc.flux.values
    exp_all, sec_all, prod_all = lc.exp.values, lc.sector.values, lc["product"].values
    h = min(1.5 * hw * P, 0.3 * P)
    Tc0 = T0 + ph_ecl * P
    E = np.arange(np.floor((t_all.min() - Tc0) / P) - 1, np.ceil((t_all.max() - Tc0) / P) + 2)
    rows = []
    Tlin = Tc0 + E * P
    sh = shift_fn(Tlin) if shift_fn is not None else np.zeros(len(E))
    lo = np.searchsorted(t_all, Tlin + sh - h)
    hi = np.searchsorted(t_all, Tlin + sh + h)
    for e, i0, i1, s_e in zip(E, lo, hi, sh):
        if i1 - i0 < 10:
            continue
        Tp = Tc0 + e * P + s_e
        t, f, ex = t_all[i0:i1], f_all[i0:i1], exp_all[i0:i1]
        cad = np.median(ex)
        expect = h / cad
        if (np.sum(t < Tp) < min_frac * expect) or (np.sum(t > Tp) < min_frac * expect):
            continue
        # core coverage: points within +-0.5 hw on both sides
        core = 0.5 * hw * P
        if np.sum((t > Tp - core) & (t < Tp)) < 2 or np.sum((t < Tp + core) & (t > Tp)) < 2:
            continue
        step = max(search / 250, 1e-5)
        r = fit_eclipse(t, f, ex, tmpl, P, T0 + s_e, Tp, h, search, step)
        if r is None:
            continue
        rows.append(dict(kind=kind, epoch=int(e), T=Tp + r["dt"], sig=r["sig"], amp=r["amp"],
                         rms=r["rms"], n=r["n"], sector=int(sec_all[i0]), product=prod_all[i0],
                         cad_s=round(cad * 86400)))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ O-C models
def ltte(t, P3, T3, e, w, asini_d):
    """Light-travel-time delay (days) of the EB in its orbit around the triple's centre of
    mass (Irwin 1952 form, zero-point removed by the linear ephemeris terms):
    dt = asini/c * [ (1-e^2)/(1+e cos v) * sin(v+w) + e sin w ]   (asini_d in light-days)."""
    M = 2 * np.pi * (t - T3) / P3
    Ecc = M.copy()
    for _ in range(30):
        Ecc = Ecc - (Ecc - e * np.sin(Ecc) - M) / (1 - e * np.cos(Ecc))
    v = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(Ecc / 2), np.sqrt(1 - e) * np.cos(Ecc / 2))
    return asini_d * ((1 - e ** 2) / (1 + e * np.cos(v)) * np.sin(v + w) + e * np.sin(w))


def amp_ltte(asini_d, e, w):
    """Borkovits-convention LTTE amplitude A = asini/c * sqrt(1 - e^2 cos^2 w) (days)."""
    return asini_d * np.sqrt(1 - e ** 2 * np.cos(w) ** 2)


def min_m3(f_m, m12):
    """Minimum tertiary mass (sin i = 1) from f(m3) = m3^3/(m12+m3)^2, solar masses."""
    lo, hi = 0.0, 100.0
    for _ in range(100):
        m = 0.5 * (lo + hi)
        if m ** 3 / (m12 + m) ** 2 < f_m:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def mass_function(asini_d, P3_d):
    """f(m3) in Msun from a12 sin i (light-days) and P3 (days)."""
    a_au = asini_d * 86400 * 299792.458 / 1.495978707e8
    return a_au ** 3 / (P3_d / 365.25) ** 2
