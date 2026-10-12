"""My own loaders. Every function returns dict(t=BTJD, f=normalised flux, e=err, q=bool good) ."""
import glob, os
import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
BKJD_TO_BTJD = 2454833.0 - 2457000.0   # = -2167


def _find(pattern):
    m = sorted(glob.glob(os.path.join(DATA, "**", pattern), recursive=True))
    if not m:
        raise FileNotFoundError(pattern)
    return m[0]


def _norm(t, f, e, q):
    ok = q & np.isfinite(t) & np.isfinite(f)
    med = np.nanmedian(f[ok])
    return dict(t=t[ok], f=f[ok] / med, e=(e[ok] / med if e is not None else np.full(ok.sum(), np.nan)))


def k2_everest(c):
    f = _find(f"hlsp_everest_k2_llc_211924561-c{c:02d}_kepler_v2.0_lc.fits")
    d = fits.getdata(f, 1)
    t = d["TIME"] + BKJD_TO_BTJD
    q = d["QUALITY"] == 0
    # EVEREST FLUX is de-systematised flux; error column: FRAW_ERR
    return _norm(t, np.asarray(d["FLUX"], float), np.asarray(d["FRAW_ERR"], float), q)


def k2_k2sff(c):
    f = _find(f"hlsp_k2sff_k2_lightcurve_211924561-c{c:02d}_kepler_v1_llc.fits")
    d = fits.getdata(f, 1)  # BESTAPER
    t = d["T"] + BKJD_TO_BTJD
    q = np.ones(len(t), bool)
    return _norm(t, np.asarray(d["FCOR"], float), None, q)


def k2_k2sc(c=5):
    f = _find(f"hlsp_k2sc_k2_llc_211924561-c{c:02d}_kepler_v2_lc.fits")
    d = fits.getdata(f, 1)
    t = d["time"] + BKJD_TO_BTJD
    q = (d["quality"] == 0) & (d["mflags"] == 0) if "mflags" in d.columns.names else d["quality"] == 0
    # k2sc: flux is systematics-corrected; add back the time trend component median-free
    flux = np.asarray(d["flux"], float) - np.asarray(d["trposi"], float) + np.nanmedian(d["trposi"])
    return _norm(t, flux, np.asarray(d["error"], float), q)


def k2_sap(c, short=False):
    f = _find(f"ktwo211924561-c{c:02d}_{'slc' if short else 'llc'}.fits")
    d = fits.getdata(f, 1)
    t = d["TIME"] + BKJD_TO_BTJD
    q = d["SAP_QUALITY"] == 0
    return _norm(t, np.asarray(d["SAP_FLUX"], float), np.asarray(d["SAP_FLUX_ERR"], float), q)


def k2_pdc(c, short=False):
    f = _find(f"ktwo211924561-c{c:02d}_{'slc' if short else 'llc'}.fits")
    d = fits.getdata(f, 1)
    t = d["TIME"] + BKJD_TO_BTJD
    q = d["SAP_QUALITY"] == 0
    return _norm(t, np.asarray(d["PDCSAP_FLUX"], float), np.asarray(d["PDCSAP_FLUX_ERR"], float), q)


def tess_spoc2min(s, col="PDCSAP_FLUX", strict=True):
    f = _find(f"tess*-s{s:04d}-0000000238606535-*-s_lc.fits")
    d = fits.getdata(f, 1)
    q = d["QUALITY"] == 0 if strict else np.ones(len(d), bool)
    return _norm(np.asarray(d["TIME"], float), np.asarray(d[col], float), np.asarray(d[col + "_ERR"], float), q)


def tess_tessspoc(s, col="PDCSAP_FLUX"):
    f = _find(f"hlsp_tess-spoc_tess_phot_0000000238606535-s{s:04d}_tess_v1_lc.fits")
    d = fits.getdata(f, 1)
    q = d["QUALITY"] == 0
    return _norm(np.asarray(d["TIME"], float), np.asarray(d[col], float), np.asarray(d[col + "_ERR"], float), q)


def tess_qlp(s):
    f = _find(f"hlsp_qlp_tess_ffi_s{s:04d}-0000000238606535_tess_v01_llc.fits")
    d = fits.getdata(f, 1)
    col = "KSPSAP_FLUX" if "KSPSAP_FLUX" in d.columns.names else "DET_FLUX"
    q = d["QUALITY"] == 0
    return _norm(np.asarray(d["TIME"], float), np.asarray(d[col], float), np.asarray(d[col + "_ERR"], float), q)


def tess_qlp_sap(s):
    f = _find(f"hlsp_qlp_tess_ffi_s{s:04d}-0000000238606535_tess_v01_llc.fits")
    d = fits.getdata(f, 1)
    q = d["QUALITY"] == 0
    return _norm(np.asarray(d["TIME"], float), np.asarray(d["SAP_FLUX"], float), None, q)


def tess_tars(s):
    f = _find(f"hlsp_tars_tess_ffi_s{s:04d}-0000000238606535_tess_v01_lc.fits")
    d = fits.getdata(f, 1)
    return _norm(np.asarray(d["TIME"], float), np.asarray(d["FLUX"], float), None, np.ones(len(d), bool))


def tess_quality_raw(s, product="spoc2min"):
    """Return (time, quality) unfiltered, to locate momentum dumps (bit 6 = 32, bit 8 = 128)."""
    if product == "spoc2min":
        f = _find(f"tess*-s{s:04d}-0000000238606535-*-s_lc.fits")
    else:
        f = _find(f"hlsp_tess-spoc_tess_phot_0000000238606535-s{s:04d}_tess_v1_lc.fits")
    d = fits.getdata(f, 1)
    return np.asarray(d["TIME"], float), np.asarray(d["QUALITY"], int)


def flatten(t, f, window=1.0, mask=None, niter=3):
    """Robust running-median detrend (window in days). mask=True points are excluded from the trend."""
    m = np.ones(len(t), bool) if mask is None else ~mask
    trend = np.empty(len(t))
    for _ in range(niter):
        for i in range(len(t)):
            sel = m & (np.abs(t - t[i]) < window / 2)
            trend[i] = np.median(f[sel]) if sel.sum() > 5 else np.nan
        r = f / trend - 1
        s = 1.4826 * np.nanmedian(np.abs(r[m] - np.nanmedian(r[m])))
        m = m & ~(np.abs(r) > 4 * s)
    # fill gaps by interpolation
    bad = ~np.isfinite(trend)
    if bad.any() and (~bad).sum() > 2:
        trend[bad] = np.interp(t[bad], t[~bad], trend[~bad])
    return f / trend, trend
