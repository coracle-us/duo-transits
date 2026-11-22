"""WISE light curves and disk phase at TESS sector times (stage 2 helpers)."""
import numpy as np
import pandas as pd
from scipy.interpolate import splrep, splev

BUILD_SLOPE = -0.02   # mag / 100 d; below -> building (W1 brightening)
DISS_SLOPE = 0.02     # above -> dissipating
QUIET_EXCESS = 0.05   # mag
MAX_GAP = 250.0       # days between bracketing WISE visits


# ---------------------------------------------------------------- IRSA single exposures
def irsa_cuts(t, w1="w1mpro", w2="w2mpro"):
    """Jian et al. 2024 cuts on a single-exposure table (pandas)."""
    ok = (t["qual_frame"] > 0) & (t["nb"] == 1) & (t["saa_sep"] > 0)
    ok &= t["cc_flags"].astype(str).str.startswith("0")
    ok &= t["moon_masked"].astype(str).str.startswith("0")
    ok &= np.isfinite(t[w1]) & np.isfinite(t[w2])
    return t[ok]


def bin_visits(mjd, w1, w2, gap=10.0):
    """Split where gaps exceed `gap` days; median mag and standard error per visit."""
    o = np.argsort(mjd)
    mjd, w1, w2 = np.asarray(mjd)[o], np.asarray(w1)[o], np.asarray(w2)[o]
    edges = np.where(np.diff(mjd) > gap)[0] + 1
    rows = []
    for idx in np.split(np.arange(len(mjd)), edges):
        if len(idx) == 0:
            continue
        n = len(idx)
        def se(x):
            return np.std(x, ddof=1) / np.sqrt(n) * 1.2533 if n > 1 else np.nan
        a, b = w1[idx], w2[idx]
        rows.append(dict(mjd=np.mean(mjd[idx]), w1mpro_ep=np.median(a), w1sigmpro_ep=se(a),
                         w2mpro_ep=np.median(b), w2sigmpro_ep=se(b), n_exp=n))
    v = pd.DataFrame(rows)
    if len(v):
        v["w1_w2"] = v.w1mpro_ep - v.w2mpro_ep
    return v


# ---------------------------------------------------------------- spline and disk phase
def fit_spline(mjd, w1, s=0.04):
    """Cubic FITPACK smoothing spline, unit weights (as in Jia et al.)."""
    o = np.argsort(mjd)
    x, y = np.asarray(mjd, float)[o], np.asarray(w1, float)[o]
    keep = np.r_[True, np.diff(x) > 0]
    return splrep(x[keep], y[keep], k=3, s=s)


class DiskCurve:
    """A star's W1 spline plus its visit times, for evaluating disk phase at any time."""

    def __init__(self, mjd, w1, w1w2=None, s=0.04):
        o = np.argsort(mjd)
        self.mjd = np.asarray(mjd, float)[o]
        self.w1 = np.asarray(w1, float)[o]
        self.w1w2 = None if w1w2 is None else np.asarray(w1w2, float)[o]
        self.tck = fit_spline(self.mjd, self.w1, s=s)
        grid = np.linspace(self.mjd[0], self.mjd[-1], 4000)
        # faintest spline value, only within 100 d of a real visit: the spline swings
        # freely across the 2011-2013 WISE hibernation gap (figures/spline_check.png)
        near = np.abs(grid[:, None] - self.mjd[None, :]).min(axis=1) <= 100.0
        self.w1_quiet = float(np.max(splev(grid[near], self.tck)))
        self.t0, self.t1 = self.mjd[0], self.mjd[-1]

    def at(self, t_mjd):
        t = np.atleast_1d(np.asarray(t_mjd, float))
        w1_fit = splev(t, self.tck)
        d = splev(t, self.tck, der=1) * 100.0
        inside = (t >= self.t0) & (t <= self.t1)
        j = np.clip(np.searchsorted(self.mjd, t), 1, len(self.mjd) - 1)
        gap = self.mjd[j] - self.mjd[j - 1]
        gap = np.where(inside, gap, np.inf)
        excess = self.w1_quiet - w1_fit
        lab = np.where(d < BUILD_SLOPE, "building",
              np.where(d > DISS_SLOPE, "dissipating",
              np.where(excess < QUIET_EXCESS, "quiescent", "plateau")))
        out = dict(w1_fit=w1_fit, disk_excess=excess, dw1dt=d, wise_gap=gap,
                   phase_label=lab, wise_usable=inside & (gap <= MAX_GAP))
        if self.w1w2 is not None:
            near = np.abs(self.mjd[None, :] - t[:, None]).argmin(axis=1)
            out["w1w2"] = self.w1w2[near]
        return out
