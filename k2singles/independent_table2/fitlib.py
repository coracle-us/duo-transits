"""Trapezoid event fitting with local polynomial baseline; errors from scaled covariance and
residual-permutation (prayer-bead) bootstrap; the larger of the two is adopted."""
import numpy as np
from scipy.optimize import least_squares
import lib2


def fit_event(t, f, e, tg, T14g, depthg, exptime, win, npoly=2, fix_shape=None):
    m = np.abs(t - tg) < win
    t, f, e = t[m], f[m], e[m]
    if len(t) < 10:
        return None
    x = t - tg

    def model(p, tt=t):
        t0, dep, T14, frac = p[:4]
        c = p[4:]
        T12 = frac * T14 / 2
        base = np.polyval(c, tt - tg)
        return base * lib2.trapz_binned(tt, t0, dep, T14, T12, exptime)

    p0 = [tg, depthg, T14g, 0.5] + [0] * npoly + [1.0]
    lo = [tg - T14g, 0, 0.2 * T14g, 0.02] + [-np.inf] * (npoly + 1)
    hi = [tg + T14g, 1, 3 * T14g, 1.0] + [np.inf] * (npoly + 1)
    best = None
    for fr in (0.2, 0.5, 0.9):
        for dt in (-0.1, 0, 0.1):
            p0[0] = tg + dt * T14g
            p0[3] = fr
            r = least_squares(lambda p: (f - model(p)) / e, p0, bounds=(lo, hi), x_scale='jac')
            if best is None or r.cost < best.cost:
                best = r
    p = best.x
    res = f - model(p)
    chi2r = np.sum((res / e) ** 2) / (len(t) - len(p))
    try:
        J = best.jac
        cov = np.linalg.inv(J.T @ J) * max(chi2r, 1)
        perr = np.sqrt(np.diag(cov))
    except Exception:
        perr = np.full(len(p), np.nan)
    # prayer-bead: cyclically shift residuals
    mod = model(p)
    t0s = []
    nshift = min(60, len(t) - 1)
    for k in np.linspace(1, len(t) - 1, nshift).astype(int):
        fs = mod + np.roll(res, k)
        r = least_squares(lambda q: (fs - model(q)) / e, p, bounds=(lo, hi), x_scale='jac')
        t0s.append(r.x[0])
    pb = np.std(t0s)
    return dict(t0=p[0], t0_err_cov=perr[0], t0_err_pb=pb, t0_err=max(perr[0], pb),
                depth=p[1], depth_err=perr[1], T14=p[2], T14_err=perr[2], T12=p[3] * p[2] / 2,
                chi2r=chi2r, n=len(t), tt=t, ff=f, mm=mod, base=np.polyval(p[4:], t - tg))
