"""
AAS82207 revision, editor Q4: validate our own TESScut aperture photometry
(singles/scripts/tesscut_photometry.py, called by pilot.add_own_photometry; reproduced in
revlib.own_photometry) against SPOC 2-min light curves of the same sector.

Targets: 11 stars, T = 8.5-15.5: the original validation case (TIC 147035057, S33, grazing
EB), eight known transiting planets (TOI CP/KP, depth 5-47 ppt) observed at 2-min in a recent
sector (S91-S107: the 200-s FFI regime in which the own photometry is actually used, sectors >= 100
in production), one non-transiting star from our sample (TIC 65290698, S105) and one faint one
(TIC 1461133, S93, T = 15.5, beyond the range of interest).

Comparisons per star (SPOC 2-min binned onto each FFI exposure, quality == 0 for both):
 (i)   absolute flux: median own/SAP_FLUX and own/PDCSAP_FLUX; aperture sizes; CROWDSAP,
       FLFRCSAP; own flux vs the TESS zero point (15,000 e-/s at T = 10, TESS Instrument Handbook);
       background per pixel (own vs SAP_BKG / N_ap)
 (ii)  normalised light curves: point-to-point scatter, 1-h binned scatter after a 1-d running
       median (transits masked), Pearson correlation of the raw 1-h binned curves
 (iii) transits: batman fit (rp, a/R*, b, T0 shift; P fixed from the TOI table; quadratic LD fixed;
       exposure-time supersampling) to own, SAP and PDCSAP; depth and per-transit mid-times
 (iv)  systematics: raw (undetrended) 1-h scatter, rms of own - SAP (1-h bins), cadence retention
Output: singles/revision/photometry_validation.csv, photometry_validation_transits.csv,
        plots/phot_<TIC>_s<sector>.png
"""
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")
REV = Path(__file__).resolve().parent
sys.path.insert(0, str(REV))
import revlib as rl  # noqa: E402

ROOT = rl.ROOT
PLOTS = REV / "plots"
LCDIR = REV / "data" / "lc"

# TIC, sector, label; ephemeris (BJD T0, P d, duration h, depth ppm) from data/raw/exofop_toi.csv
TARGETS = [
    (147035057, 33, "CTOI 147035057.01 (grazing EB; original validation case)"),
    (166836920, 105, "TOI-267.01 (CP)"),
    (144065872, 106, "TOI-105.01 (CP)"),
    (149603524, 98, "TOI-102.01 (CP)"),
    (183537452, 103, "TOI-192.01 (CP)"),
    (65290698, 105, "CTOI 65290698.01 host (no transit in S105)"),
    (294301883, 91, "TOI-774.01 (CP)"),
    (46432937, 107, "TOI-7716.01 (CP/KP)"),
    (289793076, 95, "TOI-108.01 (CP)"),
    (271898990, 92, "TOI-7501.01 (CP/KP)"),
    (1461133, 93, "CTOI 1461133 host (faint, no transit)"),
]
# 147035057: our own ephemeris (singles/results/vetting/147035057_ephemeris.json)
OWN_EPHEM = {147035057: dict(T0=1515.0661, P=37.31497, dur_h=4.53, depth_ppm=3300)}


def coords_and_ephem(tic):
    toi = pd.read_csv(ROOT / "data" / "raw" / "exofop_toi.csv")
    r = toi[toi["TIC ID"] == tic]
    if len(r):
        r = r.iloc[0]
        c = SkyCoord(r["RA"], r["Dec"], unit=(u.hourangle, u.deg))
        eph = dict(T0=r["Epoch (BJD)"] - 2457000, P=r["Period (days)"], dur_h=r["Duration (hours)"],
                   depth_ppm=r["Depth (ppm)"], Tmag=r["TESS Mag"])
        return c.ra.deg, c.dec.deg, eph
    cs = pd.read_csv(ROOT / "singles" / "data" / "ctoi_singles.csv")
    r = cs[cs.TIC == tic].iloc[0]
    eph = OWN_EPHEM.get(tic, {})
    return r.RA, r.Dec, dict(eph, Tmag=r.Tmag)


def spoc_lc(tic, sector):
    import lightkurve as lk
    d = LCDIR / str(tic)
    d.mkdir(parents=True, exist_ok=True)
    hits = list(d.glob(f"tess*-s{sector:04d}-*{tic:016d}*-s_lc.fits"))
    if not hits:
        sr = lk.search_lightcurve(f"TIC {tic}", author="SPOC", exptime=120, sector=sector)
        if len(sr) == 0:
            return None
        sr[0].download(download_dir=str(d), quality_bitmask="none")
        hits = list(d.rglob(f"tess*-s{sector:04d}-*{tic:016d}*-s_lc.fits"))
    p = hits[0]
    with fits.open(p) as h:
        t = h[1].data
        hd = h[1].header
        apim = h[2].data
        df = pd.DataFrame({"time": np.asarray(t["TIME"], float), "sap": np.asarray(t["SAP_FLUX"], float),
                           "pdc": np.asarray(t["PDCSAP_FLUX"], float), "bkg": np.asarray(t["SAP_BKG"], float),
                           "q": np.asarray(t["QUALITY"], int)})
        info = dict(CROWDSAP=hd.get("CROWDSAP"), FLFRCSAP=hd.get("FLFRCSAP"),
                    n_ap_spoc=int(np.sum((apim & 2) > 0)), spoc_file=p.name)
    return df, info


def bin_to(t_ffi, cad, ts, ys):
    """Average 2-min samples (ts sorted) within each FFI exposure; NaN if < 60% coverage."""
    lo = np.searchsorted(ts, t_ffi - cad / 2)
    hi = np.searchsorted(ts, t_ffi + cad / 2)
    need = 0.6 * cad / (120 / 86400)
    cs = np.concatenate([[0], np.cumsum(ys)])
    n = hi - lo
    out = np.where(n >= need, (cs[hi] - cs[lo]) / np.maximum(n, 1), np.nan)
    return out


def running_median(t, y, w=1.0):
    out = np.empty_like(y)
    for i, ti in enumerate(t):
        m = np.abs(t - ti) < w / 2
        out[i] = np.nanmedian(y[m])
    return out


def hour_bins(t, y):
    b = np.floor((t - t.min()) * 24).astype(int)
    s = pd.Series(y).groupby(b)
    c = s.count()
    return s.mean()[c >= max(2, int(0.5 * c.max()))]


def p2p(y):
    return 1.4826 * np.nanmedian(np.abs(np.diff(y))) / np.sqrt(2)


# ----------------------------------------------------------------- transit fitting
def batman_model(t, T0, P, rp, a, b, exp_s):
    import batman
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp, pm.a = T0, P, rp, a
    pm.inc = np.degrees(np.arccos(min(b / a, 1.0)))
    pm.ecc, pm.w, pm.u, pm.limb_dark = 0.0, 90.0, [0.4, 0.25], "quadratic"
    ss = max(1, int(round(exp_s / 120)))
    m = batman.TransitModel(pm, t, supersample_factor=ss, exp_time=exp_s / 86400) if ss > 1 else batman.TransitModel(pm, t)
    return m.light_curve(pm)


def windows(t, f, T0, P, D):
    """Locally detrended (linear, out-of-transit) data around each transit with baseline both sides."""
    ks = np.arange(np.ceil((t.min() - T0) / P), np.floor((t.max() - T0) / P) + 1)
    tt, ff, kk = [], [], []
    for k in ks:
        Tk = T0 + k * P
        w = np.abs(t - Tk) < max(1.5 * D, D / 2 + 0.25)
        oot = w & (np.abs(t - Tk) > 0.6 * D)
        if w.sum() < 10 or (oot & (t < Tk)).sum() < 3 or (oot & (t > Tk)).sum() < 3:
            continue
        if (np.abs(t - Tk) < 0.3 * D).sum() < 2:
            continue
        c = np.polyfit(t[oot] - Tk, f[oot], 1)
        tt.append(t[w])
        ff.append(f[w] / np.polyval(c, t[w] - Tk))
        kk.append(np.full(w.sum(), int(k)))
    if not tt:
        return None
    return np.concatenate(tt), np.concatenate(ff), np.concatenate(kk)


def fit_transits(t, f, eph, exp_s):
    T0, P, D = eph["T0"], eph["P"], eph["dur_h"] / 24
    W = windows(t, f, T0, P, D)
    if W is None:
        return None
    tw, fw, kw = W
    d0 = eph["depth_ppm"] * 1e-6
    a0 = max(3.0, P / (np.pi * D))

    def res(p):
        dt, rp, a, b = p
        return batman_model(tw, T0 + dt, P, rp, a, b, exp_s) - fw

    best = None
    for b0 in (0.2, 0.6, 0.9):
        try:
            r = least_squares(res, [0.0, np.sqrt(d0), a0, b0], bounds=([-0.05, 0.005, 1.5, 0.0], [0.05, 0.6, 200, 1.5]))
        except Exception:  # noqa: BLE001
            continue
        if best is None or r.cost < best.cost:
            best = r
    dt, rp, a, b = best.x
    model = batman_model(tw, T0 + dt, P, rp, a, b, exp_s)
    sig = np.std(fw - model)
    # depth = maximum model depth at native cadence (no smearing)
    tf = np.linspace(T0 + dt - D, T0 + dt + D, 2001)
    depth = 1 - batman_model(tf, T0 + dt, P, rp, a, b, 120).min()
    # per-transit mid-times with the shape fixed
    mids = {}
    for k in np.unique(kw):
        m = kw == k
        Tk = T0 + dt + k * P

        def r1(p):
            return batman_model(tw[m], T0 + dt + p[0], P, rp, a, b, exp_s) - fw[m]
        try:
            rr = least_squares(r1, [0.0], bounds=([-0.05], [0.05]))
            J = rr.jac
            cov = np.linalg.inv(J.T @ J) * np.var(rr.fun) if np.linalg.matrix_rank(J) else np.nan
            mids[int(k)] = (Tk + rr.x[0], float(np.sqrt(cov[0, 0])) if np.ndim(cov) else np.nan)
        except Exception:  # noqa: BLE001
            pass
    return dict(dt_min=dt * 1440, rp=rp, a=a, b=b, depth_ppm=depth * 1e6, rms_ppm=sig * 1e6,
                n_transits=len(np.unique(kw)), mids=mids)


# ----------------------------------------------------------------- one star
def validate(tic, sector, label):
    ra, dec, eph = coords_and_ephem(tic)
    own, info = rl.own_photometry(tic, sector, ra, dec)
    sp = spoc_lc(tic, sector)
    if own is None or sp is None:
        return None, []
    sdf, sinfo = sp
    cad = np.nanmedian(np.diff(own.time[own.good]))
    exp_s = cad * 86400
    o = own[own.good].reset_index(drop=True)
    s = sdf[(sdf.q == 0) & np.isfinite(sdf.sap) & np.isfinite(sdf.pdc) & np.isfinite(sdf.time)].sort_values("time")
    ts = s.time.to_numpy()
    sap_b = bin_to(o.time.to_numpy(), cad, ts, s.sap.to_numpy())
    pdc_b = bin_to(o.time.to_numpy(), cad, ts, s.pdc.to_numpy())
    bkg_b = bin_to(o.time.to_numpy(), cad, ts, s.bkg.to_numpy())
    both = np.isfinite(sap_b) & np.isfinite(pdc_b)
    t = o.time.to_numpy()[both]
    fo, fs, fp = o.flux.to_numpy()[both], sap_b[both], pdc_b[both]
    bo, bs = o.bkg.to_numpy()[both], bkg_b[both]
    row = dict(TIC=tic, sector=sector, label=label, Tmag=round(eph["Tmag"], 2), ffi_cadence_s=int(round(exp_s)),
               n_ap_own=info["n_ap"], n_ap_spoc=sinfo["n_ap_spoc"], CROWDSAP=sinfo["CROWDSAP"],
               FLFRCSAP=sinfo["FLFRCSAP"], n_common=int(both.sum()),
               own_cadences_kept=round(float(own.good.mean()), 3),
               spoc_cadences_kept=round(float(((sdf.q == 0) & np.isfinite(sdf.pdc)).mean()), 3))
    # (i) absolute flux
    row["own_flux_e_s"] = round(float(np.median(fo)), 1)
    row["sap_flux_e_s"] = round(float(np.median(fs)), 1)
    row["pdc_flux_e_s"] = round(float(np.median(fp)), 1)
    row["ratio_own_sap"] = round(float(np.median(fo / fs)), 4)
    row["ratio_own_pdc"] = round(float(np.median(fo / fp)), 4)
    row["ratio_own_sap_p16_p84"] = f"{np.percentile(fo / fs, 16):.4f}-{np.percentile(fo / fs, 84):.4f}"
    zp = 15000 * 10 ** (-0.4 * (eph["Tmag"] - 10))
    row["own_over_zp"] = round(float(np.median(fo)) / zp, 3)
    row["pdc_over_zp"] = round(float(np.median(fp)) / zp, 3)
    row["bkg_own_per_px"] = round(float(np.median(bo)), 2)
    row["bkg_spoc_per_px"] = round(float(np.median(bs)) / sinfo["n_ap_spoc"], 2)
    row["bkg_corr"] = round(float(np.corrcoef(bo, bs)[0, 1]), 3)
    # transit mask for the noise statistics
    mask = np.zeros(len(t), bool)
    if "P" in eph and np.isfinite(eph.get("P", np.nan)):
        ph = ((t - eph["T0"]) / eph["P"] + 0.5) % 1 - 0.5
        mask = np.abs(ph * eph["P"]) < eph["dur_h"] / 24
    stats = {}
    for lab, y in (("own", fo), ("sap", fs), ("pdc", fp)):
        yn = y / np.median(y)
        rm = running_median(t[~mask], yn[~mask], 1.0)
        det = yn[~mask] / rm
        stats[lab] = yn
        row[f"p2p_ppm_{lab}"] = round(p2p(det) * 1e6, 0)
        row[f"rms1h_det_ppm_{lab}"] = round(float(np.std(hour_bins(t[~mask], det))) * 1e6, 0)
        row[f"rms1h_raw_ppm_{lab}"] = round(float(np.std(hour_bins(t[~mask], yn[~mask]))) * 1e6, 0)
    ho = hour_bins(t[~mask], stats["own"][~mask])
    hs = hour_bins(t[~mask], stats["sap"][~mask])
    hp = hour_bins(t[~mask], stats["pdc"][~mask])
    j = ho.index.intersection(hs.index).intersection(hp.index)
    row["corr1h_own_sap"] = round(float(np.corrcoef(ho[j], hs[j])[0, 1]), 3)
    row["corr1h_own_pdc"] = round(float(np.corrcoef(ho[j], hp[j])[0, 1]), 3)
    row["rms1h_own_minus_sap_ppm"] = round(float(np.std(ho[j] - hs[j])) * 1e6, 0)
    row["rms1h_own_minus_pdc_ppm"] = round(float(np.std(ho[j] - hp[j])) * 1e6, 0)
    # level steps: split at gaps > 1 h (orbit gaps, momentum dumps, flagged stretches); the event
    # test fits a free offset per segment (gap > 0.5 d) only, so also report steps between those
    tt_ = t[~mask]
    ratio = stats["own"][~mask] / stats["sap"][~mask]
    for lab, gap in (("1h", 1 / 24), ("05d", 0.5)):
        seg = np.concatenate([[0], np.cumsum(np.diff(tt_) > gap)])
        lev = pd.Series(ratio).groupby(seg).median()
        nper = pd.Series(ratio).groupby(seg).size()
        lev = lev[nper >= 20]
        row[f"n_seg_{lab}"] = int(len(lev))
        row[f"seg_level_spread_ppt_{lab}"] = round(float((lev.max() - lev.min()) * 1e3), 1) if len(lev) > 1 else 0.0
        on = stats["own"][~mask] / pd.Series(stats["own"][~mask]).groupby(seg).transform("median").to_numpy()
        sn = stats["sap"][~mask] / pd.Series(stats["sap"][~mask]).groupby(seg).transform("median").to_numpy()
        hos, hss = hour_bins(tt_, on), hour_bins(tt_, sn)
        jj = hos.index.intersection(hss.index)
        row[f"corr1h_own_sap_segnorm_{lab}"] = round(float(np.corrcoef(hos[jj], hss[jj])[0, 1]), 3)
        row[f"rms1h_own_minus_sap_segnorm_ppm_{lab}"] = round(float(np.std(hos[jj] - hss[jj])) * 1e6, 0)
    # (iii) transits: own at FFI cadence, SPOC at native 2-min (q == 0)
    trows = []
    if "P" in eph:
        fits_ = {}
        for lab, tt, yy, ex in (("own", o.time.to_numpy(), o.flux.to_numpy(), exp_s),
                                ("sap", ts, s.sap.to_numpy(), 120), ("pdc", ts, s.pdc.to_numpy(), 120)):
            fits_[lab] = fit_transits(tt, yy / np.median(yy), eph, ex)
        for lab, fr in fits_.items():
            if fr:
                row[f"depth_ppm_{lab}"] = round(fr["depth_ppm"], 0)
                row[f"dt0_min_{lab}"] = round(fr["dt_min"], 2)
                row[f"b_{lab}"] = round(fr["b"], 2)
                row[f"n_tr_{lab}"] = fr["n_transits"]
        if fits_.get("own") and fits_.get("sap"):
            row["depth_ratio_own_sap"] = round(fits_["own"]["depth_ppm"] / fits_["sap"]["depth_ppm"], 3)
        if fits_.get("own") and fits_.get("pdc"):
            row["depth_ratio_own_pdc"] = round(fits_["own"]["depth_ppm"] / fits_["pdc"]["depth_ppm"], 3)
            ks = sorted(set(fits_["own"]["mids"]) & set(fits_["pdc"]["mids"]))
            dd = []
            for k in ks:
                (to, eo), (tp, ep_) = fits_["own"]["mids"][k], fits_["pdc"]["mids"][k]
                dd.append((to - tp) * 1440)
                trows.append(dict(TIC=tic, sector=sector, epoch_k=k, T_own=to, T_own_err_min=eo * 1440,
                                  T_pdc=tp, T_pdc_err_min=ep_ * 1440, dT_min=(to - tp) * 1440))
            if dd:
                row["n_midtimes"] = len(dd)
                row["midtime_diff_mean_min"] = round(float(np.mean(dd)), 2)
                row["midtime_diff_rms_min"] = round(float(np.std(dd)), 2)
                row["midtime_err_own_med_min"] = round(float(np.median([r["T_own_err_min"] for r in trows if r["TIC"] == tic])), 2)
    plot(tic, sector, label, t, stats, o, s, eph, row)
    return row, trows


def plot(tic, sector, label, t, stats, o, s, eph, row):
    fig, ax = plt.subplots(2, 1, figsize=(10, 6))
    ho = hour_bins(t, stats["own"])
    hs = hour_bins(t, stats["sap"])
    hp = hour_bins(t, stats["pdc"])
    tb = t.min() + (ho.index.to_numpy() + 0.5) / 24
    ax[0].plot(t, (stats["own"] - 1) * 1e3, ".", ms=1, color="0.7")
    ax[0].plot(tb, (ho.values - 1) * 1e3, "-", lw=1, label=f"own TESScut ({row['n_ap_own']} px)")
    ax[0].plot(t.min() + (hs.index.to_numpy() + 0.5) / 24, (hs.values - 1) * 1e3 - 0, "-", lw=1,
               label=f"SPOC SAP ({row['n_ap_spoc']} px)")
    ax[0].plot(t.min() + (hp.index.to_numpy() + 0.5) / 24, (hp.values - 1) * 1e3, "-", lw=1, label="SPOC PDCSAP")
    lo = np.nanpercentile((stats["own"] - 1) * 1e3, 0.5)
    hi = np.nanpercentile((stats["own"] - 1) * 1e3, 99.5)
    ax[0].set_ylim(lo - 0.3 * (hi - lo), hi + 0.3 * (hi - lo))
    ax[0].legend(fontsize=7)
    ax[0].set_ylabel("ppt (1-h bins)")
    ax[0].set_title(f"TIC {tic} S{sector} {label}: own/SAP = {row['ratio_own_sap']}, own/PDC = {row['ratio_own_pdc']}",
                    fontsize=8)
    if "P" in eph:
        P, T0 = eph["P"], eph["T0"] + row.get("dt0_min_pdc", 0) / 1440
        for lab, tt, yy, c in (("PDCSAP 2-min", s.time.to_numpy(), s.pdc.to_numpy(), "tab:green"),
                               ("own", o.time.to_numpy(), o.flux.to_numpy(), "tab:blue")):
            W = windows(tt, yy / np.median(yy), T0, P, eph["dur_h"] / 24)
            if W is None:
                continue
            ph = ((W[0] - T0) / P + 0.5) % 1 - 0.5
            ax[1].plot(ph * P * 24, (W[1] - 1) * 1e3, ".", ms=1.5 if lab == "own" else 0.8, color=c, alpha=0.6,
                       label=lab)
        ax[1].set_xlim(-2 * eph["dur_h"], 2 * eph["dur_h"])
        ax[1].set_xlabel("hours from mid-transit")
        ax[1].set_ylabel("ppt")
        ax[1].legend(fontsize=7)
        ax[1].set_title(f"depth own {row.get('depth_ppm_own')} / SAP {row.get('depth_ppm_sap')} / "
                        f"PDC {row.get('depth_ppm_pdc')} ppm; mid-time own-PDC {row.get('midtime_diff_mean_min')} "
                        f"+- {row.get('midtime_diff_rms_min')} min (rms)", fontsize=8)
    else:
        ax[1].plot(t, (stats["own"] - stats["sap"]) * 1e3, ".", ms=1)
        ax[1].set_ylabel("own - SAP (ppt)")
        ax[1].set_xlabel("BTJD")
    fig.tight_layout()
    fig.savefig(PLOTS / f"phot_{tic}_s{sector}.png", dpi=100)
    plt.close(fig)


def main():
    sel = [int(x) for x in sys.argv[1:]]
    rows, trows = [], []
    for tic, sec, label in TARGETS:
        if sel and tic not in sel:
            continue
        print(f"TIC {tic} S{sec}", flush=True)
        try:
            r, tr = validate(tic, sec, label)
        except Exception:  # noqa: BLE001
            import traceback
            print(traceback.format_exc())
            continue
        if r:
            rows.append(r)
            trows += tr
            print({k: v for k, v in r.items() if k != "label"}, flush=True)
    out = pd.DataFrame(rows)
    f = REV / "photometry_validation.csv"
    if sel and f.exists():
        old = pd.read_csv(f)
        out = pd.concat([old[~old.TIC.isin(sel)], out], ignore_index=True)
    out.to_csv(f, index=False)
    f2 = REV / "photometry_validation_transits.csv"
    tro = pd.DataFrame(trows)
    if sel and f2.exists():
        old = pd.read_csv(f2)
        tro = pd.concat([old[~old.TIC.isin(sel)], tro], ignore_index=True)
    tro.to_csv(f2, index=False)


if __name__ == "__main__":
    main()
