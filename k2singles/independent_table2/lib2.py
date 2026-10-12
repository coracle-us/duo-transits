"""Own loader / photometry / detrending utilities (written from scratch for this blind analysis)."""
import numpy as np, glob, os, warnings
from astropy.io import fits
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')


def _f(epic, pat):
    return sorted(glob.glob(os.path.join(DATA, f'EPIC{epic}', '**', pat), recursive=True))


def _norm(t, f, e, ok, src, exptime):
    m = ok & np.isfinite(t) & np.isfinite(f)
    if e is None:
        e = np.full_like(f, np.nan)
    t, f, e = t[m], f[m], e[m]
    med = np.nanmedian(f)
    f = f / med
    e = e / med
    if not np.all(np.isfinite(e)) or np.nanmedian(e) <= 0:
        e = np.full_like(f, 1.4826 * np.nanmedian(np.abs(np.diff(f))) / np.sqrt(2))
    o = np.argsort(t)
    return dict(src=src, t=t[o], f=f[o], e=e[o], exptime=exptime)


def load_all(epic, include_cut=True):
    segs = []
    # K2 mission products, long + short cadence
    for fn in _f(epic, 'ktwo*_llc.fits') + _f(epic, 'ktwo*_slc.fits'):
        d = fits.open(fn)[1].data
        camp = os.path.basename(fn).split('-')[1].split('_')[0]
        cad = 'SC' if 'slc' in fn else 'LC'
        ex = 1 / 1440 if cad == 'SC' else 0.0204
        t = d['TIME'] - 2167.0
        q = d['SAP_QUALITY']
        ok = (q == 0)
        segs.append(_norm(t, d['PDCSAP_FLUX'].astype(float), d['PDCSAP_FLUX_ERR'].astype(float), ok, f'K2-{camp}-PDCSAP-{cad}', ex))
        segs.append(_norm(t, d['SAP_FLUX'].astype(float), d['SAP_FLUX_ERR'].astype(float), ok, f'K2-{camp}-SAP-{cad}', ex))
        segs.append(_norm(t, d['SAP_FLUX'].astype(float), d['SAP_FLUX_ERR'].astype(float), np.ones_like(ok), f'K2-{camp}-SAPallq-{cad}', ex))
    for fn in _f(epic, 'hlsp_everest*.fits'):
        d = fits.open(fn)[1].data
        camp = fn.split('-c')[-1][:2]
        t = d['TIME'] - 2167.0
        ok = d['QUALITY'] == 0
        segs.append(_norm(t, d['FCOR'].astype(float), None, ok, f'K2-c{camp}-EVEREST', 0.0204))
        # EVEREST flags deep in-eclipse points as outliers; keep a version without the quality cut
        segs.append(_norm(t, d['FCOR'].astype(float), None, np.ones(len(t), bool), f'K2-c{camp}-EVERESTallq', 0.0204))
    for fn in _f(epic, 'hlsp_k2sff*.fits'):
        d = fits.open(fn)[1].data
        camp = fn.split('-c')[-1][:2]
        t = d['T'] - 2167.0
        ok = np.ones(len(t), bool)
        segs.append(_norm(t, d['FCOR'].astype(float), None, ok, f'K2-c{camp}-K2SFF', 0.0204))
    for fn in _f(epic, 'hlsp_k2sc*.fits'):
        d = fits.open(fn)[1].data
        camp = fn.split('-c')[-1][:2]
        t = d['time'] - 2167.0
        ok = d['quality'] == 0
        fl = d['flux'] - d['trposi'] + np.nanmedian(d['trposi'])
        segs.append(_norm(t, fl, d['error'], ok, f'K2-c{camp}-K2SC', 0.0204))
    # TESS
    for fn in _f(epic, 'tess*-s_lc.fits'):
        d = fits.open(fn)[1].data
        s = int(os.path.basename(fn).split('-s')[1][:4])
        segs.append(_norm(d['TIME'], d['PDCSAP_FLUX'].astype(float), d['PDCSAP_FLUX_ERR'].astype(float), d['QUALITY'] == 0, f'TESS-S{s}-SPOC120', 120 / 86400))
        segs.append(_norm(d['TIME'], d['SAP_FLUX'].astype(float), d['SAP_FLUX_ERR'].astype(float), d['QUALITY'] == 0, f'TESS-S{s}-SPOC120SAP', 120 / 86400))
    for fn in _f(epic, 'hlsp_tess-spoc*.fits'):
        h = fits.open(fn)
        d = h[1].data
        s = int(fn.split('-s')[-1][:4])
        ex = h[1].header.get('TIMEDEL', 200 / 86400)
        segs.append(_norm(d['TIME'], d['PDCSAP_FLUX'].astype(float), d['PDCSAP_FLUX_ERR'].astype(float), d['QUALITY'] == 0, f'TESS-S{s}-TESSSPOC', ex))
    for fn in _f(epic, 'hlsp_qlp*.fits'):
        d = fits.open(fn)[1].data
        s = int(fn.split('_s')[-1][:4])
        ex = 600 / 86400 if s < 56 else 200 / 86400
        if 'DET_FLUX' in d.names:
            fk, ek = 'DET_FLUX', 'DET_FLUX_ERR'
        else:
            fk, ek = 'KSPSAP_FLUX', 'KSPSAP_FLUX_ERR'
        segs.append(_norm(d['TIME'], d[fk].astype(float), d[ek].astype(float), d['QUALITY'] == 0, f'TESS-S{s}-QLP', ex))
        segs.append(_norm(d['TIME'], d['SAP_FLUX'].astype(float), None, d['QUALITY'] == 0, f'TESS-S{s}-QLPSAP', ex))
    for fn in _f(epic, 'hlsp_tars*.fits'):
        d = fits.open(fn)[1].data
        s = int(fn.split('_s')[-1][:4])
        ex = 600 / 86400 if s < 56 else 200 / 86400
        segs.append(_norm(d['TIME'].astype(float), d['FLUX'].astype(float), None, np.ones(len(d), bool), f'TESS-S{s}-TARS', ex))
    if include_cut:
        for fn in _f(epic, '*astrocut.fits'):
            s = int(os.path.basename(fn).split('-s')[1][:4])
            r = cut_phot(fn)
            if r['ok'].sum() < 10:
                print(f'  [cut_phot] {os.path.basename(fn)}: target not on silicon / no good cadences -> skipped')
                continue
            ex =1800 / 86400 if s < 27 else (600 / 86400 if s < 56 else 200 / 86400)
            segs.append(_norm(r['t'], r['f'], None, r['ok'], f'TESS-S{s}-CUT', ex))
    return segs


def cut_phot(fn, rap=None, return_info=False):
    """Own aperture photometry on a TESScut cube. Target at the cube centre."""
    h = fits.open(fn)
    d = h[1].data
    t = d['TIME'].astype(float)
    fl = d['FLUX'].astype(float)
    q = d['QUALITY']
    ny, nx = fl.shape[1:]
    cy, cx = (ny - 1) / 2, (nx - 1) / 2
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(yy - cy, xx - cx)
    good = np.isfinite(t) & np.all(np.isfinite(fl[:, r <= 3]), axis=1)
    med = np.nanmedian(fl[good], axis=0)
    outer = r > 6
    bkgmask = outer & (med < np.nanpercentile(med[outer], 50))
    bkg = np.nanmedian(fl[:, bkgmask], axis=1)
    img = med - np.nanmedian(med[bkgmask])
    if rap is None:
        peak = np.nanmax(img[r <= 1.5])
        ap = (r <= 2.5) & (img > 0.1 * peak)
        if ap.sum() < 3:
            ap = r <= 1.5
    else:
        ap = r <= rap
    f = np.nansum(fl[:, ap] - bkg[:, None], axis=1)
    ok = good & (q == 0) & np.isfinite(f) & (f > 0)
    # reject strong scattered-light cadences (background far above its typical level)
    bm = np.nanmedian(bkg[ok]) if ok.any() else 0
    bmad = 1.4826 * np.nanmedian(np.abs(bkg[ok] - bm)) if ok.any() else 0
    ok &= bkg < bm + 10 * max(bmad, 1e-3) + 0.2 * abs(bm)
    # star must actually be on silicon (science columns 45-2092, rows 1-2048; cutout centre = ref pixel + 10)
    ccol = h[1].header.get('1CRV4P', 1000) + (nx - 1) // 2
    crow = h[1].header.get('2CRV4P', 1000) + (ny - 1) // 2
    if not (45 <= ccol <= 2092 and 1 <= crow <= 2048):
        ok[:] = False
    # and the centre pixels must be clearly above background
    if not (np.nanmax(img[r <= 1.5]) > 5 * 1.4826 * np.nanmedian(np.abs(img[bkgmask]))):
        ok[:] = False
    out = dict(t=t, f=f, ok=ok, bkg=bkg, napix=int(ap.sum()))
    if return_info:
        out['ap'] = ap
        out['img'] = img
        out['q'] = q
    return out


def detrend(t, f, win=1.5, mask=None):
    """Running median detrend with window `win` days, ignoring masked (in-event) points."""
    trend = np.full_like(f, np.nan)
    use = np.ones_like(f, bool) if mask is None else ~mask
    tu, fu = t[use], f[use]
    lo = np.searchsorted(tu, t - win / 2)
    hi = np.searchsorted(tu, t + win / 2)
    for i in range(len(t)):
        if hi[i] - lo[i] > 5:
            trend[i] = np.median(fu[lo[i]:hi[i]])
    bad = ~np.isfinite(trend)
    if bad.any() and (~bad).any():
        trend[bad] = np.interp(t[bad], t[~bad], trend[~bad])
    return f / trend, trend


def trapz(t, t0, depth, T14, T12):
    """Symmetric trapezoid: total duration T14, ingress duration T12 (days)."""
    T12 = np.clip(T12, 1e-4, T14 / 2)
    x = np.abs(t - t0)
    m = np.ones_like(t)
    full = x <= T14 / 2 - T12
    m[full] = 1 - depth
    ing = (x > T14 / 2 - T12) & (x < T14 / 2)
    m[ing] = 1 - depth * (T14 / 2 - x[ing]) / T12
    return m


def trapz_binned(t, t0, depth, T14, T12, exptime):
    if exptime < 0.004:
        return trapz(t, t0, depth, T14, T12)
    n = 9
    off = (np.arange(n) - (n - 1) / 2) / n * exptime
    return np.mean([trapz(t + o, t0, depth, T14, T12) for o in off], axis=0)
