"""Transit model + likelihood helpers (batman), written for this analysis."""
import numpy as np
import batman

PREF = 100.0  # reference period (d) used only to parametrise the transit shape


def model(t, tc, k, b, aR, u, exptime_d=None, P=PREF, ss=None):
    p = batman.TransitParams()
    p.t0 = tc; p.per = P; p.rp = k; p.a = aR
    inc_cos = b / aR
    if inc_cos >= 1:
        return np.ones_like(t)
    p.inc = np.degrees(np.arccos(inc_cos)); p.ecc = 0.; p.w = 90.
    p.u = list(u); p.limb_dark = "quadratic"
    if exptime_d is not None and exptime_d > 5 / 1440.:
        n = ss or max(3, int(round(exptime_d / (2 / 1440.))))
        m = batman.TransitModel(p, t, supersample_factor=n, exp_time=exptime_d)
    else:
        m = batman.TransitModel(p, t)
    return m.light_curve(p)


def profile_baseline(t, f, mod, tc, deg=2):
    """Best-fit polynomial baseline multiplying the model (linear least squares). Returns baseline."""
    x = t - tc
    A = np.vstack([mod * x ** i for i in range(deg + 1)]).T
    c, *_ = np.linalg.lstsq(A, f, rcond=None)
    return np.vstack([x ** i for i in range(deg + 1)]).T @ c


def t14(P, aR, k, b):
    arg = np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / aR / np.sqrt(1 - (b / aR) ** 2)
    return P / np.pi * np.arcsin(np.clip(arg, 0, 1))


def t23(P, aR, k, b):
    arg = np.sqrt(np.clip((1 - k) ** 2 - b ** 2, 0, None)) / aR / np.sqrt(1 - (b / aR) ** 2)
    return P / np.pi * np.arcsin(np.clip(arg, 0, 1))


G = 6.674e-8
def rho_star(P, aR):
    """Stellar density (g/cc) for circular orbit."""
    Ps = P * 86400.
    return 3 * np.pi / (G * Ps ** 2) * aR ** 3
