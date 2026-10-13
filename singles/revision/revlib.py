"""
Shared helpers for the AAS82207 revision checks (audit of the 34 not-recovered events, Q6, and
validation of our TESScut aperture photometry, Q4).

Nothing here writes to data/raw/: new downloads go to singles/revision/data/ (gitignored).
Existing caches in data/raw/ (light curves, TESScut cutouts) are read-only inputs.

The aperture photometry here is a line-for-line copy of singles/scripts/tesscut_photometry.py
(`photometry`), returning the raw aperture flux in e-/s plus diagnostics instead of writing a
normalised CSV into data/raw/. Keep the two in sync.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parents[2]
REV = Path(__file__).resolve().parent
DATA = REV / "data"
sys.path.insert(0, str(ROOT / "scripts"))
import alias_pipeline as ap  # noqa: E402

RAW_LC = ROOT / "data" / "raw" / "lightcurves"
RAW_CUT = ROOT / "data" / "raw" / "tesscut"
EPOCH_ERR = 0.01  # same as singles/scripts/pilot.py
SIZE = 15         # same cutout size as scripts/diff_image.py


# ----------------------------------------------------------------- TESScut cutouts
def get_cutout(tic, ra, dec, sector):
    """Existing cache (data/raw/tesscut, read only) first, else download to revision/data."""
    for d in (RAW_CUT / str(tic), DATA / "tesscut" / str(tic)):
        hits = list(d.glob(f"*s{sector:04d}*.fits")) if d.exists() else []
        if hits:
            return hits[0]
    import lightkurve as lk
    d = DATA / "tesscut" / str(tic)
    d.mkdir(parents=True, exist_ok=True)
    sr = lk.search_tesscut(f"{ra} {dec}", sector=sector)
    if len(sr) == 0:
        return None
    tpf = sr[0].download(cutout_size=SIZE, download_dir=str(d))
    return Path(tpf.path)


def own_photometry(tic, sector, ra, dec):
    """Same algorithm as tesscut_photometry.photometry. Returns (DataFrame, info).
    DataFrame: time, flux (e-/s, background-subtracted aperture sum), bkg (e-/s/pixel),
    quality (TESScut QUALITY), good (the mask tesscut_photometry keeps)."""
    path = get_cutout(tic, ra, dec, sector)
    if path is None:
        return None, None
    with fits.open(path) as h:
        tab, wcs = h[1].data, WCS(h[2].header)
    t, flux, q = tab["TIME"], tab["FLUX"].astype(float), tab["QUALITY"]
    tx, ty = wcs.all_world2pix(ra, dec, 0)
    yy, xx = np.mgrid[0:flux.shape[1], 0:flux.shape[2]]
    med = np.nanmedian(flux, axis=0)
    bkg_pix = med <= np.nanpercentile(med, 40)
    bkg = np.nanmedian(flux[:, bkg_pix], axis=1)
    apm = (np.hypot(xx - tx, yy - ty) <= 1.6) & (med > np.nanpercentile(med, 60))
    f = np.nansum(flux[:, apm] - bkg[:, None], axis=1)
    good = np.isfinite(t) & np.isfinite(f) & (q == 0) & (f > 0)
    df = pd.DataFrame({"time": np.asarray(t, float), "flux": f, "bkg": bkg, "quality": np.asarray(q, int),
                       "good": good})
    info = dict(path=str(path), n_ap=int(apm.sum()), tx=float(tx), ty=float(ty),
                ap_pix=[(int(a), int(b)) for a, b in zip(*np.where(apm))])
    return df, info


# ----------------------------------------------------------------- light-curve products
FLUX_COLS = {  # columns of the stored CSVs worth measuring, per product family
    "SPOC": ["pdcsap_flux", "sap_flux"],
    "TESS-SPOC": ["pdcsap_flux", "sap_flux"],
    "QLP": ["sap_flux", "kspsap_flux"],
    "TGLC": ["cal_psf_flux", "cal_aper_flux"],
    "GSFC-ELEANOR-LITE": ["corr_flux", "raw_flux", "pca_flux"],
}


def local_products(tic):
    """{(sector, author, exptime): path} for the stored CSVs in data/raw/lightcurves."""
    out = {}
    for f in (RAW_LC / str(tic)).glob("s*.csv.gz"):
        m = re.fullmatch(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", f.name)
        if m:
            out[(int(m.group(1)), m.group(2), int(m.group(3)))] = f
    return out


def pipeline_product(tic, sector):
    """The product the pipeline (alias_pipeline.load_lightcurves with pilot.PRODUCTS) would
    choose for this sector from the files on disk, or None."""
    avail = {}
    for (sec, author, exp), f in local_products(tic).items():
        if sec != sector:
            continue
        avail[f"SPOC_{exp}" if author == "SPOC" else author] = f
    for k in ap.PRODUCT_RANK + ["OWN"]:
        if k in avail:
            return k
    return None


def segments(t, gap=ap.SEG_GAP):
    return np.concatenate([[0], np.cumsum(np.diff(t) > gap)]) if len(t) else np.array([], int)


# ----------------------------------------------------------------- measurements
def measure(t, f, epoch, D, delta, qok=None, qlp_trim=False, scan=1.0):
    """Measure the catalogued event in one light curve.
    pipe_*: alias_pipeline.test_window at the catalogued epoch with the pipeline's +-0.03 d
            offset scan (EPOCH_ERR = 0.01 d), quality==0 only (as the pipeline), i.e. the same
            test as pilot.py's event check.
    all_*:  same test, but keeping quality-flagged cadences (finite flux only).
    scan_*: the most significant box dip of duration D within +-scan d of the epoch (step D/10).
    """
    t, f = np.asarray(t, float), np.asarray(f, float)
    fin = np.isfinite(t) & np.isfinite(f)
    out = {}
    for lab, m in (("pipe", fin & (qok if qok is not None else True)), ("all", fin)):
        tt, ff = t[m], f[m]
        o = np.argsort(tt)
        tt, ff = tt[o], ff[o]
        sg = segments(tt)
        if qlp_trim and lab == "pipe" and len(tt):
            s0 = pd.Series(tt).groupby(sg).transform("min").values
            keep = tt - s0 >= ap.QLP_TRIM
            tt, ff, sg = tt[keep], ff[keep], sg[keep]
        out[f"{lab}_n_in"] = int(np.sum(np.abs(tt - epoch) < D / 2))
        if len(tt) < 20:
            out[f"{lab}_verdict"] = "no-data"
            continue
        ff = ff / np.nanmedian(ff)
        sane = np.abs(ff - 1) < max(0.3, 5 * delta)  # drop gross outliers (e.g. eleanor 1e17 cadences)
        tt, ff, sg = tt[sane], ff[sane], sg[sane]
        r = ap.test_window(tt, ff, epoch, EPOCH_ERR, D, delta, sg)
        out[f"{lab}_verdict"] = r["verdict"]
        for k in ("depth", "snr", "offset", "cover", "beta"):
            if k in r:
                out[f"{lab}_{k}"] = float(r[k])
        if lab == "all":
            best = None
            for T in np.arange(epoch - scan, epoch + scan + 1e-9, D / 10):
                rr = ap.test_window(tt, ff, T, 0.0, D, delta, sg)
                if "snr" in rr and rr.get("cover", 0) >= 0.5 and (best is None or rr["snr"] > best["snr"]):
                    best = dict(rr, T=T)
            if best is not None:
                out.update(scan_T=best["T"], scan_dt_h=(best["T"] - epoch) * 24, scan_depth=best["depth"],
                           scan_snr=best["snr"])
    if "pipe_n_in" in out and "all_n_in" in out:
        out["n_flagged_in"] = out["all_n_in"] - out["pipe_n_in"]
    return out


def sector_scan(t, f, D, delta, q_ok=None, step_frac=0.25):
    """Single-event box search over one sector: returns the strongest dip (time, depth, SNR)
    and the number of distinct dips with depth >= 0.5 delta and SNR >= 5."""
    t, f = np.asarray(t, float), np.asarray(f, float)
    m = np.isfinite(t) & np.isfinite(f)
    if q_ok is not None:
        m &= q_ok
    t, f = t[m], f[m]
    o = np.argsort(t)
    t, f = t[o], f[o]
    if len(t) < 50:
        return {}
    f = f / np.nanmedian(f)
    sane = np.abs(f - 1) < max(0.3, 5 * delta)
    t, f = t[sane], f[sane]
    sg = segments(t)
    rows = []
    for T in np.arange(t.min(), t.max(), D * step_frac):
        j = np.searchsorted(t, T)
        if j == 0 or j >= len(t) or min(T - t[j - 1], t[j] - T) > D / 2:
            continue
        r = ap.test_window(t, f, T, 0.0, D, delta, sg)
        if "snr" in r and r.get("cover", 0) >= 0.8:
            rows.append((T, r["depth"], r["snr"]))
    if not rows:
        return {}
    g = pd.DataFrame(rows, columns=["T", "depth", "snr"])
    top = g.loc[g.snr.idxmax()]
    hits = g[(g.depth >= 0.5 * delta) & (g.snr >= 5)].sort_values("snr", ascending=False)
    dips = []
    for _, r in hits.iterrows():
        if all(abs(r["T"] - x) > D for x in dips):
            dips.append(r["T"])
    return dict(sec_top_T=float(top["T"]), sec_top_depth=float(top["depth"]), sec_top_snr=float(top["snr"]),
                sec_n_dips=len(dips), sec_dips=" ".join(f"{x:.3f}" for x in sorted(dips)[:12]))
