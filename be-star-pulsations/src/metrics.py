"""Per-sector pulsation and flicker metrics (stage 3)."""
import numpy as np
from astropy.timeseries import LombScargle
from astropy.stats import sigma_clip
from scipy.interpolate import LSQUnivariateSpline

BIN_D = 30.0 / 1440.0
FMIN_SPEC, FMAX = 0.05, 24.0
FMIN = 0.25                 # all metrics start here (5-d-knot spline kills < ~0.2 c/d)
OVERSAMPLE = 10
BANDS = {"B1": (0.25, 0.8), "B2": (0.8, 2.0), "B3": (2.0, 4.0), "B4": (4.0, 12.0)}
LOG_GRID = np.geomspace(0.25, 24.0, 600)
SNR_MIN, NPEAK_MAX = 4.0, 10


# ---------------------------------------------------------------- light-curve prep
def prepare(t, f, q=None):
    """Drop flagged/NaN, normalise, 5-sigma clip, bin to 30 min, 5-d-knot spline detrend.
    Returns (t_binned, relative flux in ppt)."""
    t, f = np.asarray(t, float), np.asarray(f, float)
    ok = np.isfinite(t) & np.isfinite(f)
    if q is not None:
        ok &= np.asarray(q) == 0
    t, f = t[ok], f[ok]
    if len(t) < 50:
        return t, f
    f = f / np.nanmedian(f)
    m = sigma_clip(f, sigma=5, maxiters=None).mask
    t, f = t[~m], f[~m]
    tb, fb = bin30(t, f)
    return tb, detrend(tb, fb) * 1e3


def bin30(t, f):
    k = np.floor((t - t[0]) / BIN_D + 0.5).astype(int)  # bins centred on the first point
    uk, inv, cnt = np.unique(k, return_inverse=True, return_counts=True)
    fs = np.bincount(inv, weights=f)
    ts = np.bincount(inv, weights=t)
    # at native cadences shorter than 30 min, require at least half the expected points
    cad = np.median(np.diff(t)) if len(t) > 1 else BIN_D
    need = max(1, int(0.5 * round(BIN_D / cad)))
    keep = cnt >= need
    return (ts / cnt)[keep], (fs / cnt)[keep]


def running_median(t, f, window):
    out = np.empty_like(f)
    lo = np.searchsorted(t, t - window / 2)
    hi = np.searchsorted(t, t + window / 2, side="right")
    for i in range(len(t)):
        out[i] = np.median(f[lo[i]:hi[i]])
    return out


def segments(t, gap=0.5):
    edges = np.where(np.diff(t) > gap)[0] + 1
    return np.split(np.arange(len(t)), edges)


def detrend(t, f, knot=5.0, gap=0.5):
    """Divide by a least-squares cubic spline with 5-d knots, per segment (split at gaps
    > 0.5 d). Chosen over a 5-d running median, whose response ripples by +/-30% at
    0.25-0.5 c/d and inflates a 1.1 c/d sinusoid by 10% (tests/test_synthetic.py).
    Spline transfer: 0.88 at 0.25 c/d, 0.97 at 0.5, 0.99 at 1.1, 0.05 at 0.1 c/d."""
    r = np.empty_like(f)
    for idx in segments(t, gap):
        tt, ff = t[idx], f[idx]
        if len(tt) < 10:
            r[idx] = ff / np.median(ff) - 1.0
            continue
        nk = int(np.floor((tt[-1] - tt[0]) / knot))
        kn = np.linspace(tt[0], tt[-1], nk + 1)[1:-1] if nk >= 1 else []
        r[idx] = ff / LSQUnivariateSpline(tt, ff, kn, k=3)(tt) - 1.0
    return r


# ---------------------------------------------------------------- spectra
def freq_grid(t):
    T = t.max() - t.min()
    return np.arange(FMIN_SPEC, FMAX, 1.0 / (OVERSAMPLE * T))


def amp_spectrum(t, y, freq, method="fast"):
    """Amplitude spectrum in the units of y: A = sqrt(4 P / N), P = classical periodogram."""
    p = LombScargle(t, y, fit_mean=False, center_data=True, normalization="psd").power(
        freq, method=method)
    return np.sqrt(4.0 * p / len(t))


def sinefit(t, y, freqs):
    """Joint linear least-squares sinusoids at fixed frequencies; returns model, amplitudes."""
    X = [np.ones_like(t)]
    for f in freqs:
        X += [np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t)]
    X = np.array(X).T
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    amps = np.hypot(c[1::2], c[2::2])
    return X @ c, amps


def refine(t, y, f0, df):
    fine = np.linspace(f0 - df, f0 + df, 41)
    a = amp_spectrum(t, y, fine, method="cython")
    return fine[np.argmax(a)]


def local_noise(freq, amp, f0, half=1.0):
    m = (np.abs(freq - f0) <= half) & (freq >= FMIN)
    return np.median(amp[m])


def band_rms(freq, amp, lo, hi, T):
    """rms (same units as amp) of the signal in [lo, hi), from the oversampled spectrum:
    variance = sum A^2/2 over independent frequencies (spacing 1/T)."""
    m = (freq >= lo) & (freq < hi)
    df = freq[1] - freq[0]
    return np.sqrt(np.sum(amp[m] ** 2 / 2.0) * df * T)


def to_log_grid(freq, amp):
    """rms of A within log bins (interpolate where a bin holds no grid point)."""
    edges = np.sqrt(LOG_GRID[1:] * LOG_GRID[:-1])
    edges = np.r_[LOG_GRID[0] ** 2 / edges[0], edges, LOG_GRID[-1] ** 2 / edges[-1]]
    k = np.searchsorted(edges, freq) - 1
    out = np.interp(LOG_GRID, freq, amp)
    ok = (k >= 0) & (k < len(LOG_GRID))
    s = np.bincount(k[ok], weights=amp[ok] ** 2, minlength=len(LOG_GRID))
    n = np.bincount(k[ok], minlength=len(LOG_GRID))
    have = n > 0
    out[have] = np.sqrt(s[have] / n[have])
    return out


def smooth_1d(t, r, width=1.0):
    out = np.empty_like(r)
    for idx in segments(t, 0.5):
        tt, rr = t[idx], r[idx]
        lo = np.searchsorted(tt, tt - width / 2)
        hi = np.searchsorted(tt, tt + width / 2, side="right")
        cs = np.r_[0, np.cumsum(rr)]
        out[idx] = (cs[hi] - cs[lo]) / (hi - lo)
    return out


def sector_metrics(t, y):
    """t in days, y in ppt (already prepared). Returns (metrics dict, peaks list, log spectrum)."""
    T = t.max() - t.min()
    freq = freq_grid(t)
    df = freq[1] - freq[0]
    amp = amp_spectrum(t, y, freq)
    use = freq >= FMIN
    m = {"n_bins": len(t), "baseline_d": T}
    # strongest peak and its SNR (noise after removing it)
    i1 = np.argmax(np.where(use, amp, -1))
    f1 = refine(t, y, freq[i1], df)
    model, a = sinefit(t, y, [f1])
    res1 = y - model
    amp_r1 = amp_spectrum(t, res1, freq)
    m.update(f1=f1, A1=a[0], snr1=a[0] / local_noise(freq, amp_r1, f1))
    # prewhitening, up to 10 peaks with SNR >= 4
    peaks, fs, res, amp_r = [], [], y.copy(), amp.copy()
    for _ in range(NPEAK_MAX):
        i = np.argmax(np.where(use, amp_r, -1))
        fnew = refine(t, res, freq[i], df)
        model, a = sinefit(t, y, fs + [fnew])
        res_new = y - model
        amp_new = amp_spectrum(t, res_new, freq)
        snr = a[-1] / local_noise(freq, amp_new, fnew)
        if snr < SNR_MIN:
            break
        fs.append(fnew)
        res, amp_r = res_new, amp_new
        peaks.append(dict(freq=fnew, amp=a[-1], snr=snr))
    # refit amplitudes jointly for the final set
    if fs:
        _, a = sinefit(t, y, fs)
        for p, aa in zip(peaks, a):
            p["amp"] = aa
    for b, (lo, hi) in BANDS.items():
        m[b] = band_rms(freq, amp, lo, hi, T)
    m["n_peaks"] = len(peaks)
    m["red_noise"] = np.median(amp_r[(freq >= 0.5) & (freq <= 3.0)])
    m["flicker_rms"] = np.std(smooth_1d(t, res, 1.0))
    m["white_noise"] = np.median(amp[(freq >= 15) & (freq <= 24)])
    return m, peaks, to_log_grid(freq, amp)
