"""Shared helpers: loading SPOC/QLP light curves and TESScut photometry, detrending,
transit model. All times BTJD."""
import os, glob, json, warnings
import numpy as np
from astropy.io import fits

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "out")

TARGETS = {
    # paper (Salinas et al.) single-transit fit; ra/dec from TIC v8.2
    275647095: dict(t0=2315.5171, depth=10298e-6, dur=8.50 / 24, ra=179.337620453092, dec=-45.5465975684301),
    241550976: dict(t0=2339.4447, depth=18366e-6, dur=5.81 / 24, ra=206.723556237836, dec=-51.5670489710803),
    107113048: dict(t0=2312.9442, depth=8227e-6, dur=7.51 / 24, ra=183.075785771868, dec=-33.2327098857705),
    230133616: dict(t0=2396.1307, depth=5703e-6, dur=9.20 / 24, ra=252.267345399536, dec=69.6063588727013),
}


def _clean(t, f, e, q):
    m = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (q == 0) & (e > 0)
    t, f, e = t[m], f[m], e[m]
    med = np.nanmedian(f)
    return t, f / med, e / med


def load_lc_files(tic):
    """Return list of dicts {sector, source, t, f, e, exptime} from SPOC/QLP HLSP files."""
    out = []
    pats = [f"*tess-spoc*{tic:016d}*s00*_lc.fits", f"*qlp*s0*{tic:016d}*_llc.fits"]
    files = []
    for p in pats:
        files += glob.glob(os.path.join(DATA, "lk", "**", p), recursive=True)
    for fn in sorted(files):
        with fits.open(fn) as h:
            hd0 = h[0].header
            d = h[1].data
            sec = int(hd0["SECTOR"])
            if "tess-spoc" in fn:
                src = "TESS-SPOC"
                t, f, e, q = d["TIME"], d["PDCSAP_FLUX"], d["PDCSAP_FLUX_ERR"], d["QUALITY"]
            else:
                src = "QLP"
                # SAP flux (not KSPSAP) so that our own detrending treats long transits gently
                t, f, q = d["TIME"], d["SAP_FLUX"], d["QUALITY"]
                e = d["KSPSAP_FLUX_ERR"] if "KSPSAP_FLUX_ERR" in d.names else np.full(len(t), np.nan)
                e = np.where(np.isfinite(e), e, np.nanmedian(np.abs(np.diff(f))) / 1.4)
            t, f, e = _clean(np.array(t, float), np.array(f, float), np.array(e, float), np.array(q))
            exptime = np.nanmedian(np.diff(t)) if len(t) > 2 else np.nan
            out.append(dict(sector=sec, source=src, t=t, f=f, e=e, exptime=exptime, file=os.path.basename(fn)))
    return out


def tesscut_files(tic):
    return sorted(glob.glob(os.path.join(DATA, "tesscut_phot", f"tic{tic}_s*.npz")))


def load_tesscut_phot(tic):
    out = []
    for fn in tesscut_files(tic):
        z = np.load(fn)
        sec = int(z["sector"])
        out.append(dict(sector=sec, source="TESScut", t=z["t"], f=z["f"], e=z["e"],
                        exptime=float(np.nanmedian(np.diff(z["t"]))), file=os.path.basename(fn)))
    return out


def detrend(t, f, window, mask=None):
    """Wotan biweight trend; points with mask=True are excluded from the trend fit and the
    trend is linearly interpolated across them."""
    from wotan import flatten
    if mask is None:
        mask = np.zeros(len(t), bool)
    g = ~mask
    trg = flatten(t[g], f[g], window_length=window, method="biweight", return_trend=True,
                  break_tolerance=0.5, edge_cutoff=0.0)[1]
    ok = np.isfinite(trg)
    tr = np.interp(t, t[g][ok], trg[ok])
    return f / tr, tr


def batman_model(t, t0, per, rp, a, inc, u=(0.4, 0.25), exptime=None):
    import batman
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a, p.inc = t0, per, rp, a, inc
    p.ecc, p.w, p.u, p.limb_dark = 0.0, 90.0, list(u), "quadratic"
    if exptime and exptime > 0.004:
        n = int(np.ceil(exptime / (2.0 / 1440)))
        m = batman.TransitModel(p, t, supersample_factor=max(n, 3), exp_time=exptime)
    else:
        m = batman.TransitModel(p, t)
    return m.light_curve(p)


def stellar():
    return {int(k): v for k, v in json.load(open(os.path.join(OUT, "stellar.json"))).items()}
