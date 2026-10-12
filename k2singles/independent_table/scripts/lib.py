"""Loaders, own TESScut photometry, and local-fit event search. Written from scratch for the
blind analysis in k2singles/independent_table. All times BTJD (= BJD_TDB - 2457000)."""
import os, glob, re, warnings
import numpy as np
from astropy.io import fits
warnings.filterwarnings("ignore")

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MAST = os.path.join(ROOT, "data", "mast")
CUT = os.path.join(ROOT, "data", "tesscut", "tesscut")
BKJD2BTJD = -2167.0

# Kepler/K2 quality bits we treat as hard-bad (lightkurve-like default + thruster firing)
K2_BAD = 1130799 & ~2048   # lightkurve default minus 'impulsive outlier' (can flag real in-eclipse points)
# TESS hard-bad bits (attitude tweak, safe mode, coarse/earth point, desat, manual exclude,
# impulsive outlier, straylight-2? not used)
TESS_BAD = 17087 & ~512    # lightkurve default minus 'impulsive outlier'


def _files(pat):
    return sorted(glob.glob(os.path.join(MAST, "**", pat), recursive=True))


def _seg(src, t, f, kind, camp=None, extra=None):
    t = np.asarray(t, float); f = np.asarray(f, float)
    m = np.isfinite(t) & np.isfinite(f) & (f != 0)
    t, f = t[m], f[m]
    if len(t) < 20:
        return None
    o = np.argsort(t); t, f = t[o], f[o]
    f = f / np.nanmedian(f)
    d = dict(src=src, t=t, f=f, kind=kind, camp=camp)
    if extra: d.update(extra)
    return d


def load_k2(epic, other_ids=()):
    segs = []
    ids = [epic] + list(other_ids)
    for i in ids:
        for fn in _files(f"ktwo{i}-c*_llc.fits"):
            c = re.search(r"-c(\d+)_", fn).group(1)
            h = fits.open(fn); d = h[1].data
            t = d["TIME"] + BKJD2BTJD; q = d["SAP_QUALITY"].astype(int)
            good = (q & K2_BAD) == 0
            s = _seg(f"K2 C{c} PDCSAP ({i})", t[good], d["PDCSAP_FLUX"][good], "K2pdc", c)
            if s: segs.append(s)
            # SAP: keep everything with finite flux (quality ignored) for inspection
            s = _seg(f"K2 C{c} SAP-all ({i})", t, d["SAP_FLUX"], "K2sap", c)
            if s: segs.append(s)
        for fn in _files(f"hlsp_everest_k2_llc_{i}-c*_lc.fits"):
            c = re.search(r"-c(\d+)_", fn).group(1)
            d = fits.open(fn)[1].data
            q = d["QUALITY"].astype(int)
            good = (q & K2_BAD) == 0
            s = _seg(f"K2 C{c} EVEREST ({i})", d["TIME"][good] + BKJD2BTJD, d["FLUX"][good], "K2ev", c)
            if s: segs.append(s)
        for fn in _files(f"hlsp_k2sff_k2_lightcurve_{i}-c*_llc.fits"):
            c = re.search(r"-c(\d+)_", fn).group(1)
            d = fits.open(fn)[1].data
            s = _seg(f"K2 C{c} K2SFF ({i})", d["T"] + BKJD2BTJD, d["FCOR"], "K2sff", c)
            if s: segs.append(s)
    # EVEREST / K2SFF / raw SAP: clip to the fine-point span covered by PDCSAP of the same
    # campaign (C0 EVEREST contains ~40 d of coarse-point data with flat filler values)
    span = {}
    for x in segs:
        if x["kind"] == "K2pdc":
            span[x["camp"]] = (x["t"][0] - 0.5, x["t"][-1] + 0.5)
    for x in segs:
        if x["kind"] in ("K2ev", "K2sff", "K2sap") and x["camp"] in span:
            a, b = span[x["camp"]]; m = (x["t"] >= a) & (x["t"] <= b)
            x["t"], x["f"] = x["t"][m], x["f"][m] / np.nanmedian(x["f"][m])
    return segs


def _btjd(t):
    t = np.asarray(t, float)
    return t - 2457000 if np.nanmedian(t) > 2e6 else t


def load_tess_pipeline(tic):
    segs = []
    t10 = f"{int(tic):016d}"
    files = set(_files(f"*{t10}*.fits"))
    for fn in _files("hlsp_tglc_*.fits") + _files("hlsp_t16_*.fits") + _files("hlsp_cdips_*.fits"):
        try:
            if int(fits.getheader(fn, 0).get("TICID", -1)) == int(tic): files.add(fn)
        except Exception:
            pass
    for fn in sorted(files):
        b = os.path.basename(fn).lower()
        try:
            h = fits.open(fn)
        except Exception:
            continue
        m = re.search(r"s(\d{4})", b); sec = int(m.group(1)) if m else None
        try:
            segs += _parse_tess(fn, b, h, sec)
        except Exception as ex:
            print("   parse fail", b, ex)
    return segs


def _parse_tess(fn, b, h, sec):
    segs = []
    if True:
        if b.startswith("tess") and b.endswith("_lc.fits"):  # SPOC 2-min / 20-s
            d = h[1].data; q = d["QUALITY"].astype(int); g = (q & TESS_BAD) == 0
            sec = int(h[0].header.get("SECTOR", sec or -1))
            s = _seg(f"TESS S{sec} SPOC-2min PDCSAP", d["TIME"][g], d["PDCSAP_FLUX"][g], "spoc", sec)
            if s: segs.append(s)
        elif "tess-spoc" in b:
            d = h[1].data; q = d["QUALITY"].astype(int); g = (q & TESS_BAD) == 0
            s = _seg(f"TESS S{sec} TESS-SPOC PDCSAP", d["TIME"][g], d["PDCSAP_FLUX"][g], "tspoc", sec)
            if s: segs.append(s)
        elif "_qlp_" in b:
            d = h[1].data; q = d["QUALITY"].astype(int); g = q == 0
            s = _seg(f"TESS S{sec} QLP SAP", d["TIME"][g], d["SAP_FLUX"][g], "qlp", sec)
            if s: segs.append(s)
            dc = "KSPSAP_FLUX" if "KSPSAP_FLUX" in d.columns.names else "DET_FLUX"
            s = _seg(f"TESS S{sec} QLP {dc}", d["TIME"][g], d[dc][g], "qlpk", sec)
            if s: segs.append(s)
        elif "_tglc_" in b:
            d = h[1].data
            g = (d["TESS_flags"] == 0) & (d["TGLC_flags"] == 0)
            s = _seg(f"TESS S{sec} TGLC aper", d["time"][g], d["cal_aper_flux"][g], "tglc", sec)
            if s: segs.append(s)
        elif "_tars_" in b:
            d = h[1].data
            s = _seg(f"TESS S{sec} TARS", d["TIME"], d["FLUX"], "tars", sec)
            if s: segs.append(s)
        elif "_cdips_" in b:
            d = h[1].data; sec = int(h[0].header.get("SECTOR", sec or -1))
            g = np.ones(len(d), bool)
            if "IRQ2" in d.columns.names: g = d["IRQ2"] == "G"
            col = "PCA2" if "PCA2" in d.columns.names else "IRM2"
            mag = np.asarray(d[col][g], float)
            s = _seg(f"TESS S{sec} CDIPS {col}", _btjd(d["TMID_BJD"][g]), 10 ** (-0.4 * (mag - np.nanmedian(mag))), "cdips", sec)
            if s: segs.append(s)
        elif "_t16_" in b:
            d = h[1].data; sec = int(h[0].header.get("SECTOR", sec or -1))
            col = "TFA2" if "TFA2" in d.columns.names else "IRM2"
            mag = np.asarray(d[col], float)
            s = _seg(f"TESS S{sec} T16 {col}", _btjd(d["TMID_BJD"]), 10 ** (-0.4 * (mag - np.nanmedian(mag))), "t16", sec)
            if s: segs.append(s)
        elif "eleanor" in b:
            d = h[1].data; g = d["QUALITY"] == 0
            s = _seg(f"TESS S{sec} eleanor-lite", d["TIME"][g], d["CORR_FLUX"][g], "eleanor", sec)
            if s: segs.append(s)
    return segs


def tesscut_files(ra, dec):
    return sorted(glob.glob(os.path.join(CUT, f"tess-s*_{ra:.6f}_{dec:.6f}_11x11_astrocut.fits")))


def own_tesscut_phot(fn):
    """Simple aperture photometry on an 11x11 TESScut cutout. Target at the WCS position.
    Background: per-cadence median of the faintest 40% of pixels outside r=3 px.
    Aperture: chosen among circular r=1,1.5,2,2.5 px apertures around the target pixel to
    minimise the 1-h-binned point-to-point scatter after a 1-d median detrend."""
    h = fits.open(fn)
    sec = int(h[0].header.get("SECTOR", re.search(r"tess-s(\d{4})", fn).group(1)))
    d = h[1].data
    t = d["TIME"]; F = d["FLUX"].astype(float); q = d["QUALITY"].astype(int)
    from astropy.wcs import WCS
    w = WCS(h[2].header)
    ra, dec = float(re.search(r"_(-?[\d.]+)_(-?[\d.]+)_11x11", fn).group(1)), float(re.search(r"_(-?[\d.]+)_(-?[\d.]+)_11x11", fn).group(2))
    x0, y0 = w.all_world2pix([[ra, dec]], 0)[0]
    ny, nx = F.shape[1:]
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = np.hypot(xx - x0, yy - y0)
    good = np.isfinite(t) & ((q & TESS_BAD) == 0) & np.all(np.isfinite(F.reshape(len(t), -1)), axis=1)
    t, F, q = t[good], F[good], q[good]
    outer = r > 3.5
    med_img = np.nanmedian(F, axis=0)
    thr = np.nanpercentile(med_img[outer], 40)
    bmask = outer & (med_img <= thr)
    bkg = np.nanmedian(F[:, bmask], axis=1)
    # drop strong scattered-light cadences: background far above its typical level
    bmed = np.nanmedian(bkg); bsd = robust_std(bkg)
    keep = bkg < bmed + 8 * bsd + 0.0
    t, F, q, bkg = t[keep], F[keep], q[keep], bkg[keep]
    Fs = F - bkg[:, None, None]
    best = None
    for rad in (1.0, 1.5, 2.0, 2.5):
        ap = r <= rad
        if ap.sum() < 1: continue
        fl = Fs[:, ap].sum(axis=1)
        if np.nanmedian(fl) <= 0: continue
        fn_ = fl / np.nanmedian(fl)
        tr = running_median(t, fn_, 1.0)
        res = fn_ / tr - 1
        sc = robust_std(bin_series(t, res, 1 / 24.)[1])
        if best is None or sc < best[0]:
            best = (sc, rad, fn_)
    if best is None or best[0] > 0.03:   # unusable (star off silicon / at edge / swamped)
        return None
    s = _seg(f"TESS S{sec} own-TESScut r={best[1]}", t, best[2], "own", sec,
             extra=dict(bkg=bkg / np.nanmedian(np.abs(bkg) + 1e-9)))
    return s


# ---------------------------------------------------------------- utilities
def robust_std(x):
    x = np.asarray(x); x = x[np.isfinite(x)]
    if len(x) < 3: return np.nan
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def running_median(t, f, w):
    out = np.empty_like(f)
    lo = np.searchsorted(t, t - w / 2); hi = np.searchsorted(t, t + w / 2)
    for i in range(len(t)):
        out[i] = np.median(f[lo[i]:hi[i]])
    return out


def bin_series(t, f, dt):
    if len(t) == 0: return np.array([]), np.array([])
    k = np.floor((t - t[0]) / dt).astype(int)
    u, inv = np.unique(k, return_inverse=True)
    s = np.bincount(inv, f); n = np.bincount(inv)
    tb = np.bincount(inv, t) / n
    return tb, s / n


def cadence(t):
    dt = np.diff(t)
    return np.median(dt) if len(dt) else np.nan


def local_fit(t, f, t0, dur, win=None, order=2, ingress_frac=0.15):
    """Fit a polynomial to out-of-event data around t0 and measure the mean depth in the
    central part of the event window. Returns dict with depth, err, coverage, n_in, n_out."""
    if win is None:
        win = max(1.5 * dur, 0.75)
    half = dur / 2
    sel = (t > t0 - half - win) & (t < t0 + half + win)
    tt, ff = t[sel], f[sel]
    cad = cadence(t)
    nexp_in = max(1, int(round(dur * (1 - 2 * ingress_frac) / cad)))
    if len(tt) < 10:
        return dict(depth=np.nan, err=np.nan, cov=0.0, cov_full=0.0, n_in=0, n_out=len(tt), nb=0, na=0)
    inn = np.abs(tt - t0) < half * (1 - 2 * ingress_frac)
    full = np.abs(tt - t0) < half
    out = np.abs(tt - t0) > half * 1.1
    nb = np.sum(out & (tt < t0)); na = np.sum(out & (tt > t0))
    cov = inn.sum() / nexp_in
    cov_full = full.sum() / max(1, int(round(dur / cad)))
    if out.sum() < 8 or inn.sum() < 2:
        return dict(depth=np.nan, err=np.nan, cov=cov, cov_full=cov_full, n_in=int(inn.sum()), n_out=int(out.sum()), nb=int(nb), na=int(na))
    x = (tt - t0)
    o = order if (nb > 5 and na > 5) else min(order, 1)
    if (nb < 3 or na < 3): o = 0
    c = np.polyfit(x[out], ff[out], o)
    model = np.polyval(c, x)
    r = ff - model
    sig = robust_std(r[out])
    # red-noise: compare scatter of means of chunks of length ~ in-transit size
    nin = inn.sum()
    chunk_scatter = []
    ro = r[out]
    if len(ro) > 3 * nin and nin > 1:
        k = len(ro) // nin
        cm = ro[:k * nin].reshape(k, nin).mean(axis=1)
        beta = max(1.0, robust_std(cm) / (sig / np.sqrt(nin))) if k >= 4 else 1.0
    else:
        beta = 1.0
    depth = -np.mean(r[inn])
    err = beta * sig / np.sqrt(nin)
    return dict(depth=depth, err=err, cov=cov, cov_full=cov_full, n_in=int(nin), n_out=int(out.sum()), nb=int(nb), na=int(na), sig=sig, beta=beta)


def scan(seg, dur, step=None, order=2):
    t, f = seg["t"], seg["f"]
    if step is None:
        step = dur / 8
    grid = np.arange(t[0] + dur / 2, t[-1] - dur / 2, step)
    rows = []
    for t0 in grid:
        # skip if no data in the central part
        lo, hi = np.searchsorted(t, [t0 - dur * 0.35, t0 + dur * 0.35])
        if hi - lo < 2: continue
        r = local_fit(t, f, t0, dur, order=order)
        r["t0"] = t0
        rows.append(r)
    return rows


# ---------------------------------------------------------------- trapezoid fit
def trapezoid(t, t0, depth, T14, T12, bl):
    """Trapezoid transit times baseline polynomial bl (callable)."""
    x = np.abs(t - t0)
    T12 = np.clip(T12, 1e-4, T14 / 2)
    m = np.ones_like(t)
    full = x <= T14 / 2 - T12
    ing = (x > T14 / 2 - T12) & (x < T14 / 2)
    m[full] = 1 - depth
    m[ing] = 1 - depth * (T14 / 2 - x[ing]) / T12
    return m


def fit_event(t, f, t0g, depg, durg, win=None, order=2, fix_shape=None):
    """Least-squares trapezoid + polynomial baseline. Returns params and 1-sigma errors
    (errors from the covariance, scaled by reduced chi2; t0 error also by red-noise beta)."""
    from scipy.optimize import least_squares
    if win is None:
        win = max(1.5 * durg, 0.75)
    sel = np.abs(t - t0g) < durg / 2 + win
    tt, ff = t[sel], f[sel]
    if len(tt) < 15:
        return None
    x = tt - t0g

    def model(p):
        t0, d, T14, T12 = p[:4]
        if fix_shape is not None:
            T14, T12 = fix_shape
        base = np.polyval(p[4:], x)
        return trapezoid(tt, t0, d, T14, T12, None) * base

    p0 = [t0g, depg, durg, durg * 0.15] + [0] * order + [1.0]
    lo = [t0g - durg, 0, durg * 0.3, durg * 0.01] + [-np.inf] * (order + 1)
    hi = [t0g + durg, 1, durg * 2.5, durg * 0.5] + [np.inf] * (order + 1)
    try:
        r = least_squares(lambda p: ff - model(p), p0, bounds=(lo, hi), x_scale="jac")
    except Exception:
        return None
    res = ff - model(r.x)
    dof = max(1, len(ff) - len(r.x))
    s2 = np.sum(res ** 2) / dof
    try:
        cov = np.linalg.inv(r.jac.T @ r.jac) * s2
        err = np.sqrt(np.diag(cov))
    except Exception:
        err = np.full(len(r.x), np.nan)
    # red-noise beta from binned residual scatter on ~ingress timescale
    nb = max(2, int(round(r.x[3] / cadence(tt))))
    beta = 1.0
    if len(res) > 5 * nb and nb > 1:
        k = len(res) // nb
        cm = res[:k * nb].reshape(k, nb).mean(axis=1)
        beta = max(1.0, robust_std(cm) / (robust_std(res) / np.sqrt(nb)))
    return dict(t0=r.x[0], depth=r.x[1], T14=r.x[2], T12=r.x[3], e_t0=err[0] * beta, e_depth=err[1] * beta,
                e_T14=err[2] * beta, rms=np.sqrt(s2), beta=beta, n=len(ff), x=r.x, tt=tt, ff=ff, model=model(r.x))
