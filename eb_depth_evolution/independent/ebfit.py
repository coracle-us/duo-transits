"""Global eclipsing-binary fit machinery shared by the TESS and K2 fits.

Parameters live in a dict; a list `free` selects those varied.  Each eclipse can carry its own
time shift 'dt_<ecl>' (model evaluated at t - dt), which absorbs eclipse-timing variations.
"""
import numpy as np
from scipy.optimize import least_squares
import ebmodel as M

BOUNDS = dict(P=(21.70, 21.85), T0=(2488.0, 2489.2), esw=(-0.7, 0.7), r1=(0.003, 0.30), r2=(0.003, 0.30),
              J=(0.01, 5.0), i2014=(60, 90), i2021=(60.0, 90.0), i2023=(60.0, 90.0), i2026=(60.0, 90.0),
              ecw2014=(-0.6, 0.2), ecw2021=(-0.6, 0.2), ecw2023=(-0.6, 0.2), ecw2026=(-0.6, 0.2))


def bounds(name):
    if name.startswith("esw"):
        return (-0.7, 0.7)
    if name.startswith("dt_"):
        return (-0.2, 0.2)
    return BOUNDS[name]


class Fitter:
    def __init__(self, df, pars, free, u1=(0.20, 0.30), u2=(0.30, 0.25), l3=0.035, nsub=None):
        self.df = df.reset_index(drop=True)
        self.pars = dict(pars)
        self.free = list(free)
        self.u1, self.u2, self.l3 = u1, u2, l3
        self.groups = []
        for ecl, d in self.df.groupby("ecl", sort=False):
            texp = d.texp.iloc[0]
            ns = nsub if nsub else (5 if texp > 400 else 3)
            if texp > 1000:
                ns = 15
            self.groups.append(dict(ecl=ecl, t=d.time.values, y=d.fn.values, e=d.err.values,
                                    season=int(d.season.iloc[0]), texp=texp, nsub=ns,
                                    tc=d.tc_rough.values[0], typ=ecl[0]))
            self.pars.setdefault("dt_" + ecl, 0.0)
        self.lo = np.array([bounds(n)[0] for n in self.free])
        self.hi = np.array([bounds(n)[1] for n in self.free])

    def full(self, th):
        p = dict(self.pars)
        p.update(zip(self.free, th))
        return p

    def season_pars(self, p, s):
        ecw = p["ecw%d" % s]
        esw = p.get("esw%d" % s, p["esw"])      # optional per-season e sin w
        e = np.hypot(ecw, esw)
        w = np.arctan2(esw, ecw)
        return dict(P=p["P"], T0=p["T0"], e=e, w=w, inc=p["i%d" % s], r1=p["r1"], r2=p["r2"], J=p["J"])

    def model_group(self, p, g, t=None, inc=None, l3=None, u1=None, u2=None):
        sp = self.season_pars(p, g["season"])
        if inc is not None:
            sp["inc"] = inc
        tt = (g["t"] if t is None else t) - p["dt_" + g["ecl"]]
        u1 = self.u1 if u1 is None else u1
        u2 = self.u2 if u2 is None else u2
        return M.lc(tt, sp["P"], sp["T0"], sp["e"], sp["w"], sp["inc"], sp["r1"], sp["r2"], sp["J"],
                    u1[0], u1[1], u2[0], u2[1], l3=self.l3 if l3 is None else l3,
                    texp=g["texp"] / 86400.0, nsub=g["nsub"])

    @staticmethod
    def baseline(m, g):
        dt = g["t"] - g["tc"]
        A = np.vstack([m, m * dt]).T / g["e"][:, None]
        c, *_ = np.linalg.lstsq(A, g["y"] / g["e"], rcond=None)
        return c

    def resid_p(self, p):
        out = []
        for g in self.groups:
            m = self.model_group(p, g)
            c = self.baseline(m, g)
            out.append((g["y"] - m * (c[0] + c[1] * (g["t"] - g["tc"]))) / g["e"])
        return np.concatenate(out)

    def resid(self, th):
        return self.resid_p(self.full(th))

    def fit(self, **kw):
        th0 = np.array([self.pars[n] for n in self.free], float)
        th0 = np.clip(th0, self.lo + 1e-9, self.hi - 1e-9)
        r = least_squares(self.resid, th0, bounds=(self.lo, self.hi), x_scale="jac", **kw)
        self.pars = self.full(r.x)
        return r

    def lnprob(self, th):
        th = np.asarray(th)
        if np.any(th < self.lo) or np.any(th > self.hi):
            return -np.inf
        r = self.resid(th)
        return -0.5 * np.sum(r ** 2)


def ecw_from_interval(P, esw, dsec, ecw0=-0.26):
    """Solve for e cos w given the secondary-minus-primary conjunction interval dsec."""
    from scipy.optimize import brentq
    def f(ecw):
        e = np.hypot(ecw, esw); w = np.arctan2(esw, ecw)
        return M.secondary_offset(P, e, w) - dsec
    return brentq(f, -0.6, 0.3)
