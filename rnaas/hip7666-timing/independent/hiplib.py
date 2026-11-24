"""HIP 7666 blind eclipse-timing library (written from scratch for this re-analysis).

Conventions
- Times in BTJD (BJD_TDB - 2457000).
- Phase phi in [-0.5, 0.5) about the primary eclipse of Lee et al. Eq. 1.
- Flux model: flux = trend(t) * (B(phi) + Puls(t)), B = binary template (~1 out of eclipse).
"""
import warnings
import numpy as np
from scipy.interpolate import LSQUnivariateSpline
from scipy.optimize import least_squares, minimize_scalar

warnings.filterwarnings("ignore")

T0_LEE = 2451769.56299 - 2457000.0   # BTJD
P_LEE = 2.37222188                   # d
ECL_HALF = 0.044                     # eclipse half-width in phase (Lee)
FIT_HALF = 0.06                      # fit window half-width in phase

# Lee et al. Table 5 frequencies (d^-1), f2 = f_orb removed -> 14 frequencies
FREQS_ALL = {
    "f1": 25.4656, "f2": 0.4203, "f3": 21.1939, "f4": 22.4278, "f5": 21.5836,
    "f6": 20.3533, "f7": 24.6340, "f8": 28.0866, "f9": 21.8091, "f10": 27.1883,
    "f11": 22.6767, "f12": 28.0289, "f13": 21.0911, "f14": 20.9630, "f15": 37.8159,
}
FREQS = np.array([v for k, v in FREQS_ALL.items() if k != "f2"])

LC_DIR = "lc"
MAIN = {
    18: dict(file="lc/hlsp_tess-spoc_tess_phot_0000000308953703-s0018_tess_v1_lc.fits", texp=1800.0),
    58: dict(file="lc/tess2022302161335-s0058-0000000308953703-0247-s_lc.fits", texp=120.0),
    85: dict(file="lc/hlsp_tess-spoc_tess_phot_0000000308953703-s0085_tess_v1_lc.fits", texp=200.0),
}


def phase(t, T0=T0_LEE, P=P_LEE):
    return np.mod((t - T0) / P + 0.5, 1.0) - 0.5


def wrap(x):
    return np.mod(x + 0.5, 1.0) - 0.5


def in_eclipse(t, half=ECL_HALF):
    ph = phase(t)
    return (np.abs(ph) < half) | (np.abs(wrap(ph - 0.5)) < half)


def load(sector, flux_column="pdcsap_flux"):
    import lightkurve as lk
    info = MAIN[sector]
    lc = lk.read(info["file"], flux_column=flux_column, quality_bitmask="default")
    lc = lc.remove_nans()
    t = np.asarray(lc.time.value, float)
    f = np.asarray(lc.flux.value, float)
    e = np.asarray(lc.flux_err.value, float)
    good = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e > 0)
    t, f, e = t[good], f[good], e[good]
    med = np.median(f)
    return t, f / med, e / med


def chunks(t, gap=0.5):
    """Index ranges of continuous stretches separated by gaps > gap days."""
    br = np.where(np.diff(t) > gap)[0]
    starts = np.r_[0, br + 1]
    ends = np.r_[br + 1, len(t)]
    return list(zip(starts, ends))


def spline_trend(t, y, use, knot_spacing):
    """Cubic least-squares spline per continuous chunk, fitted to y[use]."""
    trend = np.full_like(y, np.nan)
    for a, b in chunks(t):
        tt, yy, uu = t[a:b], y[a:b], use[a:b]
        if uu.sum() < 10:
            trend[a:b] = np.median(yy[uu]) if uu.any() else 1.0
            continue
        tu, yu = tt[uu], yy[uu]
        L = tu[-1] - tu[0]
        nk = int(round(L / knot_spacing))
        if nk >= 2:
            knots = np.linspace(tu[0], tu[-1], nk + 1)[1:-1]
            # keep only knots with data on both sides (Schoenberg-Whitney safety)
            ok = []
            for k in knots:
                if np.sum((tu > k - knot_spacing / 2) & (tu < k + knot_spacing / 2)) > 5:
                    ok.append(k)
            try:
                s = LSQUnivariateSpline(tu, yu, ok, k=3)
                trend[a:b] = s(tt)
                continue
            except Exception:
                pass
        c = np.polyfit(tu - tu.mean(), yu, 3 if L > 0.5 else 1)
        trend[a:b] = np.polyval(c, tt - tu.mean())
    return trend


def puls_design(t, freqs=FREQS, texp=None, nsub=None):
    """Design matrix [1, sin, cos ...] for fixed frequencies, optionally integrated over texp (s)."""
    t = np.asarray(t)
    cols = [np.ones_like(t)]
    if texp is None or texp <= 130:
        for f in freqs:
            arg = 2 * np.pi * f * t
            cols += [np.sin(arg), np.cos(arg)]
    else:
        # boxcar integration: sinusoid of frequency f averaged over width w -> * sinc(f w)
        w = texp / 86400.0
        for f in freqs:
            arg = 2 * np.pi * f * t
            fac = np.sinc(f * w)
            cols += [fac * np.sin(arg), fac * np.cos(arg)]
    return np.vstack(cols).T


def fit_pulsations(t, r, use, texp=None):
    X = puls_design(t[use], texp=texp)
    coef, *_ = np.linalg.lstsq(X, r[use], rcond=None)
    model = puls_design(t, texp=texp) @ coef
    model -= coef[0]  # remove constant, keep only oscillating part
    return coef, model


def fold_template(ph, y, binw=0.001):
    """Median in phase bins, then a cubic smoothing spline, periodic by wrapping.

    Smoothing: scipy UnivariateSpline with weights 1/sigma_bin and s = number of bins
    (chi^2 of the spline about the bin medians ~ its expectation). sigma_bin is the
    error of a median, 1.2533 * 1.4826 * MAD / sqrt(n).
    (make_smoothing_spline with GCV was tried first and returned a near-constant curve.)
    """
    from scipy.interpolate import UnivariateSpline
    edges = np.arange(-0.5, 0.5 + binw / 2, binw)
    idx = np.digitize(ph, edges) - 1
    xc, yc, sc = [], [], []
    for i in range(len(edges) - 1):
        m = idx == i
        if m.sum() >= 3:
            xc.append(0.5 * (edges[i] + edges[i + 1]))
            yy = y[m]
            yc.append(np.median(yy))
            mad = np.median(np.abs(yy - np.median(yy)))
            sc.append(1.2533 * 1.4826 * mad / np.sqrt(m.sum()))
    xc, yc, sc = np.array(xc), np.array(yc), np.array(sc)
    sc = np.maximum(sc, 0.5 * np.median(sc))
    sel_hi, sel_lo = xc > 0.4, xc < -0.4
    xw = np.r_[xc[sel_hi] - 1.0, xc, xc[sel_lo] + 1.0]
    yw = np.r_[yc[sel_hi], yc, yc[sel_lo]]
    sw = np.r_[sc[sel_hi], sc, sc[sel_lo]]
    spl = UnivariateSpline(xw, yw, w=1.0 / sw, k=3, s=len(xw))
    return spl, xc, yc


class Template:
    """Binary template from a phase fold; full-phase spline plus symmetrised eclipse shapes."""

    def __init__(self, spl):
        self.spl = spl
        self.centers = {}
        for kind, c0 in (("P", 0.0), ("S", 0.5)):
            self.centers[kind] = self._sym_center(c0)

    def full(self, ph):
        return self.spl(wrap(ph))

    def _sym_center(self, c0):
        x = np.linspace(0.0, ECL_HALF, 400)

        def cost(c):
            return np.sum((self.spl(wrap(c + x)) - self.spl(wrap(c - x))) ** 2)
        r = minimize_scalar(cost, bounds=(c0 - 0.01, c0 + 0.01), method="bounded",
                            options=dict(xatol=1e-9))
        return float(r.x)

    def eclipse(self, kind, x):
        """Symmetrised eclipse shape, x = phase offset from the eclipse's symmetry centre."""
        c = self.centers[kind]
        return 0.5 * (self.spl(wrap(c + x)) + self.spl(wrap(c - x)))

    def full_integrated(self, ph, texp):
        if texp <= 130:
            return self.full(ph)
        nsub = int(np.ceil(texp / 60.0)) | 1
        offs = (np.arange(nsub) - (nsub - 1) / 2) / nsub * (texp / 86400.0) / P_LEE
        return np.mean([self.full(ph + o) for o in offs], axis=0)


def eclipse_model_factory(tmpl, kind, texp):
    if texp <= 130:
        offs = np.array([0.0])
    else:
        nsub = int(np.ceil(texp / 60.0)) | 1
        offs = (np.arange(nsub) - (nsub - 1) / 2) / nsub * (texp / 86400.0)

    def model(dt, shift, scale, c0, c1):
        # dt: time minus predicted (symmetry-centre) eclipse time, days
        vals = np.mean([tmpl.eclipse(kind, (dt + o - shift) / P_LEE) for o in offs], axis=0)
        return c0 + c1 * dt + scale * (vals - 1.0)
    return model
