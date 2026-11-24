"""Shared code for the HIP 7666 timing analysis (see RULES.md).

Times in BTJD. Flux is relative (median 1). Section numbers refer to RULES.md.
"""
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.interpolate import make_lsq_spline
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parents[1]
LCDIR = HERE / "data" / "lc"
OUTDIR = HERE / "results"
OUTDIR.mkdir(exist_ok=True)

# Lee et al. 2611.17590 Eq. 1 (BJD_TDB), in BTJD.
T0_LEE = 2451769.56299 - 2457000.0
P_LEE = 2.37222188
ECL_HALF = 0.044        # eclipse half-width in phase (Lee et al. out-of-eclipse boundary)
FIT_HALF = 0.06         # eclipse fit window half-width in phase
# Lee et al. Table 5 frequencies (d^-1); f2 = f_orb is excluded (belongs to the binary model).
FREQS = {"f1": 25.4656, "f3": 21.1939, "f4": 22.4278, "f5": 21.5836, "f6": 20.3533,
         "f7": 24.6340, "f8": 28.0866, "f9": 21.8091, "f10": 27.1883, "f11": 22.6767,
         "f12": 28.0289, "f13": 21.0911, "f14": 20.9630, "f15": 37.8159}
TIMING_MODES = ["f1", "f4", "f3"]
# Binary (Lee et al. Tables 3-4)
K_A, K_B = 93.16, 136.17          # km/s
SINI = np.sin(np.radians(81.366))
M12 = 1.822 + 1.246               # Msun
C_KMS = 299792.458

PRODUCTS = {  # RULES.md section 1: main product per sector, QLP as check
    ("main", 18): ("hlsp_tess-spoc_tess_phot_0000000308953703-s0018_tess_v1_lc.fits", 1800.0),
    ("main", 58): ("tess2022302161335-s0058-0000000308953703-0247-s_lc.fits", 120.0),
    ("main", 85): ("hlsp_tess-spoc_tess_phot_0000000308953703-s0085_tess_v1_lc.fits", 200.0),
    ("qlp", 18): ("hlsp_qlp_tess_ffi_s0018-0000000308953703_tess_v01_llc.fits", 1800.0),
    ("qlp", 58): ("hlsp_qlp_tess_ffi_s0058-0000000308953703_tess_v01_llc.fits", 200.0),
    ("qlp", 85): ("hlsp_qlp_tess_ffi_s0085-0000000308953703_tess_v01_llc.fits", 200.0),
}
DEFAULT_BITMASK = 175        # lightkurve "default" quality bitmask
QLP_BAD = 4096 | 2048 | 1 | 2 | 4 | 8 | 16 | 32 | 128


def load(kind, sector, flux_col=None):
    """Return t, flux, err (median-normalised), exptime in s."""
    if kind == "sap":                       # RULES 1: SAP check of the main product
        kind, flux_col = "main", "SAP_FLUX"
    name, exptime = PRODUCTS[(kind, sector)]
    path = next(LCDIR.rglob(name))
    with fits.open(path) as h:
        d = h[1].data
        t = np.array(d["TIME"], float)
        if kind == "qlp":
            col = flux_col or "SAP_FLUX"
            f = np.array(d[col], float)
            ecol = "KSPSAP_FLUX_ERR" if "KSPSAP_FLUX_ERR" in d.columns.names else "DET_FLUX_ERR"
            e = np.array(d[ecol], float)
            q = np.array(d["QUALITY"], int)
            good = (q & QLP_BAD) == 0
        else:
            col = flux_col or "PDCSAP_FLUX"
            f = np.array(d[col], float)
            e = np.array(d[col + "_ERR"], float)
            q = np.array(d["QUALITY"], int)
            good = (q & DEFAULT_BITMASK) == 0
    ok = good & np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e > 0)
    t, f, e = t[ok], f[ok], e[ok]
    med = np.median(f)
    return t, f / med, e / med, exptime


def phase(t, t0=T0_LEE, p=P_LEE):
    """Phase in [-0.5, 0.5), primary eclipse at 0."""
    return ((t - t0) / p + 0.5) % 1.0 - 0.5


def in_eclipse(t, half=ECL_HALF):
    ph = phase(t)
    return (np.abs(ph) < half) | (np.abs(np.abs(ph) - 0.5) < half)


# ---------------------------------------------------------------- binary template
class Template:
    """Periodic LSQ cubic spline in phase: dense knots in eclipses, sparse outside."""

    def __init__(self, ph, y, w=None):
        kn = np.concatenate([np.arange(-0.5, 0.5, 0.01)])
        dense = np.concatenate([np.arange(-0.07, 0.0701, 0.0025),
                                np.arange(0.43, 0.5, 0.0025), np.arange(-0.5, -0.43, 0.0025)])
        kn = np.unique(np.round(np.concatenate([kn[(np.abs(kn) > 0.07) & (np.abs(np.abs(kn) - 0.5) > 0.07)],
                                                dense]), 6))
        # wrap the data by one period on each side so the spline is periodic in effect
        x = np.concatenate([ph - 1, ph, ph + 1])
        yy = np.concatenate([y, y, y])
        ww = None if w is None else np.concatenate([w, w, w])
        o = np.argsort(x)
        x, yy = x[o], yy[o]
        ww = None if ww is None else ww[o]
        inner = np.concatenate([kn - 1, kn, kn + 1])
        inner = inner[(inner > -1.2) & (inner < 1.2)]
        sel = (x > -1.25) & (x < 1.25)
        x, yy = x[sel], yy[sel]
        ww = None if ww is None else ww[sel]
        k = 3
        tk = np.r_[[x[0]] * (k + 1), inner[(inner > x[0]) & (inner < x[-1])], [x[-1]] * (k + 1)]
        self.spl = make_lsq_spline(x, yy, tk, k=k, w=ww)
        self.shift = {"p": 0.0, "s": 0.0}     # recentring offsets in phase

    def __call__(self, ph):
        ph = (np.asarray(ph) + 0.5) % 1.0 - 0.5
        return self.spl(ph)

    def eclipse_profile(self, which):
        """Function of phase offset from that eclipse's own symmetry centre."""
        c = (0.0 if which == "p" else 0.5) + self.shift[which]
        return lambda dph: self(c + dph)

    def recentre(self):
        """Find each eclipse's symmetry axis (RULES.md 2), store as phase offsets."""
        x = np.linspace(-ECL_HALF, ECL_HALF, 401)
        for which, c in (("p", 0.0), ("s", 0.5)):
            def asym(s):
                return np.sum((self(c + s + x) - self(c + s - x)) ** 2)
            grid = np.linspace(-0.005, 0.005, 2001)
            v = [asym(s) for s in grid]
            self.shift[which] = float(grid[int(np.argmin(v))])
        return self.shift


def smear(fun, t, exptime, nsub=None):
    """Average fun over a boxcar exposure centred on each t (exptime in s)."""
    if exptime <= 130:
        return fun(t)
    nsub = nsub or max(3, int(np.ceil(exptime / 60.0)) | 1)
    off = (np.arange(nsub) - (nsub - 1) / 2) / nsub * exptime / 86400.0
    return np.mean([fun(t + o) for o in off], axis=0)


# ---------------------------------------------------------------- trends
def chunks(t, gap=0.5):
    br = np.where(np.diff(t) > gap)[0]
    return np.split(np.arange(len(t)), br + 1)


def trend_spline(t, r, mask, spacing=1.5):
    """Per-chunk LSQ cubic spline with knots every `spacing` d, fit to r[mask]; returns trend at t."""
    out = np.ones_like(t)
    for idx in chunks(t):
        tt, rr, mm = t[idx], r[idx], mask[idx]
        if mm.sum() < 20:
            out[idx] = np.median(rr[mm]) if mm.any() else 1.0
            continue
        x, y = tt[mm], rr[mm]
        nk = int((x[-1] - x[0]) / spacing)
        inner = np.linspace(x[0], x[-1], nk + 1)[1:-1] if nk >= 2 else np.array([])
        tk = np.r_[[x[0]] * 4, inner, [x[-1]] * 4]
        try:
            s = make_lsq_spline(x, y, tk, k=3)
            out[idx] = s(np.clip(tt, x[0], x[-1]))
        except Exception:
            c = np.polyfit(x - x.mean(), y, 1)
            out[idx] = np.polyval(c, tt - x.mean())
    return out


# ---------------------------------------------------------------- pulsations
def design(t, freqs):
    cols = []
    for f in freqs:
        cols += [np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t)]
    return np.vstack(cols).T


def fit_sines(t, y, freqs, w=None):
    """Linear LSQ of sum of sines at fixed frequencies + constant. Returns (coef, cov, rms)."""
    X = np.column_stack([design(t, freqs), np.ones_like(t)])
    if w is None:
        w = np.ones_like(t)
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    res = y - X @ coef
    rms = np.sqrt(np.sum(w * res ** 2) / np.sum(w)) if len(t) else np.nan
    dof = max(len(t) - X.shape[1], 1)
    s2 = np.sum(w * res ** 2) / dof
    cov = np.linalg.pinv((X * w[:, None]).T @ X) * s2
    return coef, cov, rms


def amp_phase(a, b):
    """a sin + b cos = A sin(x + phi)."""
    return np.hypot(a, b), np.arctan2(b, a)


def sines_model(t, freqs, coef):
    return design(t, freqs) @ coef[: 2 * len(freqs)]


# ---------------------------------------------------------------- eclipse fit
def eclipse_model(par, t, tc, prof, exptime):
    s, d, c0, c1 = par
    shape = smear(lambda tt: prof((tt - tc - s) / P_LEE), t, exptime)
    return (c0 + c1 * (t - tc)) * (1.0 + d * (shape - 1.0))


def fit_eclipse(t, y, e, tc, prof, exptime):
    p0 = [0.0, 1.0, 1.0, 0.0]
    fun = lambda p: (eclipse_model(p, t, tc, prof, exptime) - y) / e
    r = least_squares(fun, p0, x_scale=[1e-3, 0.1, 1e-3, 1e-3], method="lm")
    return r.x, r


def block_bootstrap(t, y, e, tc, prof, exptime, best, nboot=500, block=0.05, seed=1):
    rng = np.random.default_rng(seed)
    model = eclipse_model(best, t, tc, prof, exptime)
    res = y - model
    n = len(t)
    dt = np.median(np.diff(t))
    bl = max(1, int(round(block / dt)))
    starts = np.arange(0, max(n - bl, 1))
    out = np.empty(nboot)
    for i in range(nboot):
        idx = np.concatenate([np.arange(s, s + bl) for s in rng.choice(starts, n // bl + 1)])[:n]
        idx = np.clip(idx, 0, n - 1)
        yb = model + res[idx]
        out[i] = fit_eclipse(t, yb, e, tc, prof, exptime)[0][0]
    return np.std(out, ddof=1)
