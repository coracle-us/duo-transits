"""Shared routines for the flare side quest.

Method summary (see ../README.md for the full write-up):
  * Data: SPOC 2-min PDCSAP light curves only (20-s used for vetting plots).
  * Quality: drop cadences whose SPOC flag hits BAD_BITS; KEEP impulsive-outlier
    (512) and cosmic-ray (64/1024) cadences, which can be real flare points.
  * Detrend: per continuous segment, iterative running median (6 h window,
    2 h if the star is strongly variable) with positive outliers masked
    between iterations; local noise = rolling robust (MAD) std of residuals.
  * Detection: Chang et al. (2015) criteria as implemented in AltaiPony
    (Ilin et al. 2021): N1=3, N2=2 (error-subtracted), N3=3 consecutive cadences.
  * ED = integral of relative flux excess (s); E_TESS = ED * L_TESS(Teff, R*)
    with L_TESS from a blackbody through a 600-1000 nm top-hat band.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter

ROOT = Path(__file__).resolve().parents[3]
SQ = ROOT / "sidequests" / "flares"
LCDIR = ROOT / "data" / "raw" / "lightcurves"
CATALOG = ROOT / "data" / "catalog" / "candidates.csv"
CACHE = SQ / "cache"

# SPOC quality bits (TESS Science Data Products Description Document)
Q_ATTITUDE_TWEAK = 1
Q_SAFE_MODE = 2
Q_COARSE_POINT = 4
Q_EARTH_POINT = 8
Q_ARGABRIGHTENING = 16
Q_DESAT = 32            # reaction-wheel desaturation = momentum dump
Q_APERTURE_COSMIC = 64
Q_MANUAL_EXCLUDE = 128
Q_DISCONTINUITY = 256
Q_IMPULSIVE_OUTLIER = 512
Q_COLLATERAL_COSMIC = 1024
Q_STRAYLIGHT = 2048
Q_SCATTERED_LIGHT = 4096
Q_PLANET_SEARCH_EXCLUDE = 8192
Q_BAD_CALIBRATION = 16384
Q_INSUFFICIENT_TARGETS = 32768

BAD_BITS = (Q_ATTITUDE_TWEAK | Q_SAFE_MODE | Q_COARSE_POINT | Q_EARTH_POINT
            | Q_ARGABRIGHTENING | Q_DESAT | Q_MANUAL_EXCLUDE | Q_STRAYLIGHT
            | Q_SCATTERED_LIGHT | Q_BAD_CALIBRATION | Q_INSUFFICIENT_TARGETS)
# Kept on purpose: 64, 256, 512, 1024, 8192 (flagged in the candidate table).

CADENCE_D = 120.0 / 86400.0
SEG_GAP_D = 0.5          # gap that splits a light curve into segments
EDGE_D = 0.5             # candidates within this of a segment edge are rejected
WIN_LONG_D = 0.25        # 6 h running median
WIN_SHORT_D = 1.0 / 12   # 2 h running median for strongly variable stars
SIG_WIN_D = 0.5          # window for the rolling robust noise estimate
N1, N2, N3 = 3.0, 2.0, 3
DUMP_PAD_D = 10.0 / 1440  # flare within 10 min of a momentum dump -> reject


# ----------------------------------------------------------------- catalog
def load_catalog() -> pd.DataFrame:
    c = pd.read_csv(CATALOG)
    return c.set_index("TIC", drop=False)


def spoc_files(tic: int, exptime: int = 120) -> dict[int, Path]:
    out = {}
    for p in sorted((LCDIR / str(tic)).glob(f"s*_SPOC_{exptime}s.csv.gz")):
        sec = int(p.name[1:4])
        if sec > 200:      # one mislabelled product exists in the archive
            continue
        out[sec] = p
    return out


# ------------------------------------------------------------- light curve
def load_sector(tic: int, sector: int, cat: pd.DataFrame, exptime: int = 120) -> dict:
    """Load one SPOC sector, apply quality mask and transit mask, find segments."""
    path = spoc_files(tic, exptime)[sector]
    d = pd.read_csv(path)
    t = d["time"].to_numpy(float)
    q = d["quality"].fillna(0).to_numpy(np.int64)
    f = d["pdcsap_flux"].to_numpy(float)
    e = d["flux_err"].to_numpy(float)
    tfin = np.isfinite(t)
    dumps = t[tfin & ((q & Q_DESAT) > 0)]

    row = cat.loc[tic]
    in_transit = np.zeros(len(t), bool)
    hw = 1.5 * float(row["Dur_h"]) / 24.0
    for tc in (row["Tc1"], row["Tc2"]):
        if np.isfinite(tc):
            in_transit |= np.abs(t - tc) < hw

    good = tfin & np.isfinite(f) & np.isfinite(e) & ((q & BAD_BITS) == 0) & ~in_transit
    # segments from the good cadences
    tg = t[good]
    seg_id = np.full(len(t), -1)
    if len(tg):
        brk = np.where(np.diff(tg) > SEG_GAP_D)[0]
        starts = np.r_[0, brk + 1]
        ends = np.r_[brk, len(tg) - 1]
        gidx = np.where(good)[0]
        segs = []
        for k, (a, b) in enumerate(zip(starts, ends)):
            seg_id[gidx[a:b + 1]] = k
            segs.append((tg[a], tg[b]))
    else:
        segs = []
    edge = np.zeros(len(t), bool)
    for (a, b) in segs:
        edge |= ((t >= a) & (t < a + EDGE_D)) | ((t <= b) & (t > b - EDGE_D))
    med = np.nanmedian(f[good]) if good.any() else np.nan
    return dict(tic=tic, sector=sector, t=t, f=f, e=e, q=q, good=good, seg_id=seg_id,
                segs=segs, edge=edge, dumps=dumps, in_transit=in_transit, med=med,
                sap=d["sap_flux"].to_numpy(float), bkg=d["sap_bkg"].to_numpy(float),
                c1=d["mom_centr1"].to_numpy(float), c2=d["mom_centr2"].to_numpy(float))


def _rolling_mad_std(x: np.ndarray, n: int) -> np.ndarray:
    s = pd.Series(x)
    med = s.rolling(n, center=True, min_periods=max(10, n // 4)).median()
    mad = (s - med).abs().rolling(n, center=True, min_periods=max(10, n // 4)).median()
    out = 1.4826 * mad.to_numpy()
    if np.all(~np.isfinite(out)):
        return np.full(len(x), 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x))))
    return pd.Series(out).bfill().ffill().to_numpy()


def _detrend_segment(t, f, win_d, n_iter=3):
    """Iterative masked running median. Returns trend, resid (relative), sigma."""
    n = len(t)
    w = max(5, int(round(win_d / CADENCE_D)) | 1)
    nsig = max(31, int(round(SIG_WIN_D / CADENCE_D)) | 1)
    mask = np.ones(n, bool)
    for _ in range(n_iter):
        tt, ff = t[mask], f[mask]
        if len(ff) < 2 * 5:
            mask[:] = True
            tt, ff = t, f
        tr_c = median_filter(ff, size=min(w, len(ff) - (1 - len(ff) % 2)), mode="nearest")
        trend = np.interp(t, tt, tr_c)
        resid = f / trend - 1.0
        sig_c = _rolling_mad_std(resid[mask], min(nsig, mask.sum()))
        sigma = np.interp(t, t[mask], sig_c)
        hi = resid > 2.5 * sigma
        # dilate: 2 cadences before, 5 after each high point (flare tails)
        newmask = np.ones(n, bool)
        for k in np.where(hi)[0]:
            newmask[max(0, k - 2):k + 6] = False
        if np.array_equal(newmask, mask):
            break
        mask = newmask
    return trend, resid, sigma, mask


def detrend(lc: dict) -> dict:
    """Detrend each segment; choose a 6 h window unless the star is strongly
    variable on short timescales, in which case use 2 h."""
    t, f, good = lc["t"], lc["f"], lc["good"]
    resid = np.full(len(t), np.nan)
    sigma = np.full(len(t), np.nan)
    trend = np.full(len(t), np.nan)
    win_used = []
    for k in range(len(lc["segs"])):
        idx = np.where(lc["seg_id"] == k)[0]
        if len(idx) < 30:
            continue
        ts, fs = t[idx], f[idx]
        tr, r, s, m = _detrend_segment(ts, fs, WIN_LONG_D)
        # white-noise level from point-to-point scatter
        p2p = 1.4826 * np.nanmedian(np.abs(np.diff(r[m]) - np.nanmedian(np.diff(r[m])))) / np.sqrt(2)
        rob = 1.4826 * np.nanmedian(np.abs(r[m] - np.nanmedian(r[m])))
        win = WIN_LONG_D
        if rob > 1.3 * p2p:
            tr2, r2, s2, m2 = _detrend_segment(ts, fs, WIN_SHORT_D)
            rob2 = 1.4826 * np.nanmedian(np.abs(r2[m2] - np.nanmedian(r2[m2])))
            if rob2 < rob:
                tr, r, s, m, win = tr2, r2, s2, m2, WIN_SHORT_D
        trend[idx], resid[idx], sigma[idx] = tr, r, s
        win_used.append(win)
    lc = dict(lc)
    lc.update(trend=trend, resid=resid, sigma=sigma, erel=lc["e"] / np.where(np.isfinite(trend), trend, np.nan),
              win_used=win_used)
    return lc


# ---------------------------------------------------------------- detection
def find_flares(lc: dict, n1=N1, n2=N2, n3=N3) -> pd.DataFrame:
    t, r, s, er = lc["t"], lc["resid"], lc["sigma"], lc["erel"]
    good = lc["good"] & np.isfinite(r) & np.isfinite(s)
    gidx = np.where(good)[0]
    if len(gidx) == 0:
        return pd.DataFrame()
    tg, rg, sg, eg = t[gidx], r[gidx], s[gidx], er[gidx]
    cand = (rg > 0) & (rg / sg > n1) & ((rg - eg) / sg > n2)
    # runs of consecutive cadences (no time gap) that pass
    contig = np.r_[False, np.diff(tg) < 1.5 * CADENCE_D]
    runs = []
    k = 0
    n = len(tg)
    while k < n:
        if cand[k]:
            j = k
            while j + 1 < n and cand[j + 1] and contig[j + 1]:
                j += 1
            if j - k + 1 >= n3:
                runs.append([k, j])
            k = j + 1
        else:
            k += 1
    # merge runs separated by <= 3 cadences
    merged = []
    for a, b in runs:
        if merged and a - merged[-1][1] <= 4 and tg[a] - tg[merged[-1][1]] < 4.5 * CADENCE_D:
            merged[-1][1] = b
        else:
            merged.append([a, b])
    ext = []
    for a, b in merged:
        # extend start back while > 1 sigma (max 5), end forward while > 0.5 sigma (max 2 h)
        st = a
        while st - 1 >= 0 and a - (st - 1) <= 5 and contig[st] and rg[st - 1] > 1.0 * sg[st - 1]:
            st -= 1
        en = b
        while en + 1 < n and (en + 1 - b) <= 60 and contig[en + 1] and rg[en + 1] > 0.5 * sg[en + 1]:
            en += 1
        if ext and st <= ext[-1][1] + 1:   # overlapping extents -> one event
            ext[-1][1] = max(en, ext[-1][1])
            ext[-1][2] += b - a + 1
        else:
            ext.append([st, en, b - a + 1])
    rows = []
    for st, en, nab in ext:
        seg = slice(st, en + 1)
        ipk = st + int(np.argmax(rg[seg]))
        ed = float(np.sum(rg[seg]) * 120.0)
        ed_err = float(np.sqrt(np.sum(sg[seg] ** 2)) * 120.0)
        rows.append(dict(
            tic=lc["tic"], sector=lc["sector"],
            t_start=tg[st], t_peak=tg[ipk], t_stop=tg[en],
            n_pts=en - st + 1, n_above=nab,
            amp=float(rg[ipk]), amp_sig=float(rg[ipk] / sg[ipk]), sigma=float(sg[ipk]),
            dur_min=float((tg[en] - tg[st]) * 1440 + 2.0),
            rise_min=float((tg[ipk] - tg[st]) * 1440), decay_min=float((tg[en] - tg[ipk]) * 1440),
            ed_s=ed, ed_err_s=ed_err,
            i0=int(gidx[st]), i1=int(gidx[en]), ipk=int(gidx[ipk]),
        ))
    return pd.DataFrame(rows)


def auto_flags(fl: pd.DataFrame, lc: dict) -> pd.DataFrame:
    """Automatic vetting flags, applied identically to real and injected events."""
    if fl.empty:
        return fl
    fl = fl.copy()
    q = lc["q"]
    near_dump, edge, shape_ok, n512, ncr = [], [], [], [], []
    for _, r in fl.iterrows():
        a, b = r.t_start - DUMP_PAD_D, r.t_stop + DUMP_PAD_D
        near_dump.append(bool(np.any((lc["dumps"] >= a) & (lc["dumps"] <= b))))
        edge.append(bool(lc["edge"][int(r.ipk)]))
        # fast rise / slower decay: peak not after the midpoint of the event
        shape_ok.append(bool(r.decay_min >= r.rise_min))
        qq = q[int(r.i0):int(r.i1) + 1]
        n512.append(int(np.sum((qq & Q_IMPULSIVE_OUTLIER) > 0)))
        ncr.append(int(np.sum((qq & (Q_APERTURE_COSMIC | Q_COLLATERAL_COSMIC)) > 0)))
    fl["near_dump"] = near_dump
    fl["edge"] = edge
    fl["shape_ok"] = shape_ok
    fl["n_flag512"] = n512
    fl["n_flag_cr"] = ncr
    return fl


def obs_time_days(lc: dict) -> float:
    """Good, non-edge, out-of-transit 2-min cadences, in days."""
    ok = lc["good"] & ~lc["edge"]
    if "resid" in lc:
        ok &= np.isfinite(lc["resid"])
    return float(ok.sum() * CADENCE_D)


# ------------------------------------------------------------------ energy
H, C, KB, SIGMA_SB = 6.62607015e-27, 2.99792458e10, 1.380649e-16, 5.670374e-5  # cgs
RSUN = 6.957e10
BAND_NM = (600.0, 1000.0)


def band_flux(teff: float, band=BAND_NM, n=400) -> float:
    """pi * integral B_lambda(T) dlambda over a top-hat band [erg s^-1 cm^-2]."""
    lam = np.linspace(band[0], band[1], n) * 1e-7  # cm
    b = 2 * H * C ** 2 / lam ** 5 / np.expm1(H * C / (lam * KB * teff))
    return float(np.pi * np.trapezoid(b, lam))


def l_tess(teff: float, rstar_rsun: float) -> float:
    return 4 * np.pi * (rstar_rsun * RSUN) ** 2 * band_flux(teff)


def bol_correction(t_flare=9000.0) -> float:
    """E_bol / E_TESS for a flare radiating as a T_flare blackbody (Shibayama et al. 2013)."""
    return SIGMA_SB * t_flare ** 4 / band_flux(t_flare)


# --------------------------------------------------------------- templates
def davenport2014(t, tpeak, fwhm, amp):
    """Davenport et al. (2014) analytic flare template; t, tpeak, fwhm in days."""
    x = (np.asarray(t) - tpeak) / fwhm
    out = np.zeros_like(x, dtype=float)
    r = (x > -1) & (x <= 0)
    out[r] = 1 + 1.941 * x[r] - 0.175 * x[r] ** 2 - 2.246 * x[r] ** 3 - 1.125 * x[r] ** 4
    d = x > 0
    out[d] = 0.6890 * np.exp(-1.600 * x[d]) + 0.3030 * np.exp(-0.2783 * x[d])
    return amp * out


def davenport_ed_per_amp_fwhm() -> float:
    x = np.linspace(-1, 200, 2_000_001)
    return float(np.trapezoid(davenport2014(x, 0.0, 1.0, 1.0), x))  # ~1.827


def binned_template(t, tpeak, fwhm, amp, exp_d=CADENCE_D, nsub=15):
    off = (np.arange(nsub) + 0.5) / nsub * exp_d - exp_d / 2
    return davenport2014(np.asarray(t)[:, None] + off[None, :], tpeak, fwhm, amp).mean(axis=1)
