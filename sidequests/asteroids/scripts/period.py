"""
Rotation period from a moving-target light curve (output of mtphot.py).

 1. Frames kept: flag == 0 (quality 0, enough static samples, not noisy) and not star-contaminated:
    star_ap (static-scene flux inside the aperture) < STAR_FRAC x the asteroid's median flux.
 2. Binning: BIN_MIN-minute bins in light-time-corrected time, median of the frames (>= 3 per bin),
    error = 1.2533 x std / sqrt(n). Flux -> mag; reduced mag Hred = m - 5 log10(r Delta).
 3. Bins > 5 robust sigma from a running median (window 2 h... set by the clip) are dropped.
 4. Phase-angle trend: Hred = c + beta*alpha (linear), fitted simultaneously with the Fourier
    series below; for the periodogram the linear-alpha fit alone is removed first.
 5. Lomb-Scargle (astropy) on 0.05-24 c/d. The strongest peak gives P_LS. Fourier fits of order
    NHARM at P_LS and 2*P_LS (period free, least squares). The rotation period is chosen
    per object with the 'odd-harmonic test': at 2*P_LS, an F-test of the full series against
    the even-harmonics-only series (which is the same as a fit at P_LS). Significant odd
    harmonics (p < 0.001) mean the two halves of the 2*P_LS light curve differ -> double-peaked,
    P_rot = 2 P_LS. Otherwise both are reported and the choice is made by eye (default 2 P_LS
    for amplitudes > 0.1 mag, the usual convention, e.g. Harris et al. 2014).
 6. Period error: least-squares covariance scaled by sqrt(reduced chi2), x2 for the red noise
    (beta factor from residuals binned to 3 h).
Plots: plots/<id>_s<sector>.png (raw + binned light curve, periodogram, folds at P_LS and 2 P_LS).
Usage: python period.py ID SECTOR [P_known_h]
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy import stats
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parents[1]
LCDIR = HERE / "data" / "lightcurves"
PLOTS = HERE / "plots"
ZP = 20.44
STAR_FRAC = 3.0
BIN_MIN = 30.0
NHARM = 4
FMIN, FMAX = 0.05, 24.0


def load(sid, sector):
    df = pd.read_csv(LCDIR / f"{sid}_s{sector:03d}.csv.gz")
    fmed = np.nanmedian(df.flux[df.flag == 0])
    df["contam"] = df.star_ap > STAR_FRAC * max(fmed, 1.0)
    keep = (df.flag == 0) & ~df.contam & np.isfinite(df.flux)
    return df, keep, fmed


def bin_lc(df, keep, bin_min=BIN_MIN):
    d = df[keep]
    tb = np.floor((d.t_lt - d.t_lt.min()) / (bin_min / 1440)).astype(int)
    rows = []
    for k, g in d.groupby(tb):
        if len(g) < 3:
            continue
        f = g.flux.values
        m = np.median(f)
        s = 1.4826 * np.median(np.abs(f - m))
        if s == 0:
            s = np.std(f)
        ok = np.abs(f - m) < 5 * s if s > 0 else np.ones(len(f), bool)
        f = f[ok]
        rows.append(dict(t=g.t_lt.values[ok].mean(), f=np.median(f),
                         ef=1.2533 * np.std(f) / np.sqrt(len(f)) if len(f) > 1 else np.nan,
                         n=len(f), r=g.r.mean(), delta=g.delta.mean(), alpha=g.alpha.mean()))
    b = pd.DataFrame(rows)
    b = b[(b.f > 0) & np.isfinite(b.ef) & (b.ef > 0)]
    b["mag"] = ZP - 2.5 * np.log10(b.f)
    b["emag"] = 1.0857 * b.ef / b.f
    b["hred"] = b.mag - 5 * np.log10(b.r * b.delta)
    return b.reset_index(drop=True)


def design(t, alpha, P, nh, t0, odd=True):
    cols = [np.ones_like(t), alpha - alpha.mean()]
    ph = 2 * np.pi * (t - t0) / P
    for k in range(1, nh + 1):
        if not odd and k % 2 == 1:
            continue
        cols += [np.cos(k * ph), np.sin(k * ph)]
    return np.column_stack(cols)


def linfit(b, P, nh, odd=True):
    t0 = b.t.min()
    A = design(b.t.values, b.alpha.values, P, nh, t0, odd)
    w = 1 / b.emag.values
    c, *_ = np.linalg.lstsq(A * w[:, None], b.hred.values * w, rcond=None)
    model = A @ c
    chi2 = np.sum(((b.hred.values - model) * w) ** 2)
    return c, model, chi2, A.shape[1]


def fit_period(b, P0, nh):
    """Nonlinear least squares with the period free."""
    t0 = b.t.min()
    w = 1 / b.emag.values

    def res(p):
        A = design(b.t.values, b.alpha.values, p[0], nh, t0)
        c, *_ = np.linalg.lstsq(A * w[:, None], b.hred.values * w, rcond=None)
        return (b.hred.values - A @ c) * w

    r = least_squares(res, [P0], x_scale=[P0 * 1e-3], diff_step=[1e-6])
    J = r.jac
    dof = len(b) - (2 + 2 * nh + 1)
    rchi2 = np.sum(r.fun ** 2) / dof
    sig = np.sqrt(np.linalg.inv(J.T @ J)[0, 0] * rchi2) if np.all(np.isfinite(J)) else np.nan
    return r.x[0], sig, rchi2


def amplitude(c, P, nh, n=2000):
    ph = np.linspace(0, 1, n)
    y = np.zeros(n)
    for k in range(1, nh + 1):
        y += c[2 * k] * np.cos(2 * np.pi * k * ph) + c[2 * k + 1] * np.sin(2 * np.pi * k * ph)
    return y.max() - y.min()


def red_beta(t, resid, emag, binh=3.0):
    """Time-averaging beta (Winn et al. 2008 style): scatter of the residuals averaged in binh-hour
    bins over the white-noise expectation std(resid)/sqrt(n). (First version divided by the
    propagated per-bin errors instead, which conflated underestimated errors with red noise;
    changed during validation on (624) Hektor.)"""
    tb = np.floor((t - t.min()) / (binh / 24)).astype(int)
    g = pd.DataFrame(dict(tb=tb, r=resid)).groupby("tb")
    m = g.r.mean()
    n = g.r.size()
    ok = n >= 3
    if ok.sum() < 5:
        return 1.0
    exp = np.std(resid) / np.sqrt(n[ok])
    return float(max(1.0, np.std(m[ok]) / np.median(exp)))


def analyse(sid, sector, P_known_h=None, label=None, plot=True, bin_min=BIN_MIN):
    df, keep, fmed = load(sid, sector)
    b = bin_lc(df, keep, bin_min)
    out = dict(id=sid, sector=sector, n_frames=len(df), n_keep=int(keep.sum()),
               frac_contam=round(float(df.contam[df.flag == 0].mean()), 3), n_bins=len(b),
               median_flux=round(fmed, 2), V_med=round(float(df.V.median()), 2),
               span_d=round(float(b.t.max() - b.t.min()), 2) if len(b) else 0)
    if len(b) < 30:
        out["status"] = "too few bins"
        return out
    # 5-sigma clip on residuals of a 6-h running median of binned points
    med = b.hred.rolling(13, center=True, min_periods=3).median()
    res = b.hred - med
    s = 1.4826 * np.nanmedian(np.abs(res))
    b = b[np.abs(res) < 6 * s].reset_index(drop=True)
    out["n_bins"] = len(b)
    out["bin_scatter_mag"] = round(float(1.4826 * np.median(np.abs(np.diff(b.hred))) / np.sqrt(2)), 4)
    # remove linear alpha trend for the periodogram
    A0 = np.column_stack([np.ones(len(b)), b.alpha - b.alpha.mean()])
    c0, *_ = np.linalg.lstsq(A0 / b.emag.values[:, None], b.hred / b.emag, rcond=None)
    y = b.hred - A0 @ c0
    freq = np.linspace(FMIN, min(FMAX, 0.5 * 1440 / bin_min), 200000)
    ls = LombScargle(b.t, y, b.emag)
    pw = ls.power(freq)
    i = np.argmax(pw)
    f1 = freq[i]
    P_ls = 1 / f1
    out["P_LS_h"] = round(P_ls * 24, 5)
    out["LS_power"] = round(float(pw[i]), 4)
    out["LS_FAP"] = float(ls.false_alarm_probability(pw[i], minimum_frequency=FMIN,
                                                     maximum_frequency=freq.max()))
    # Fourier fits at P_LS and 2 P_LS
    res = {}
    for tag, P0 in [("1x", P_ls), ("2x", 2 * P_ls)]:
        P, sP, rchi2 = fit_period(b, P0, NHARM)
        c, model, chi2, npar = linfit(b, P, NHARM)
        res[tag] = dict(P=P, sP=sP, rchi2=rchi2, c=c, chi2=chi2, npar=npar, amp=amplitude(c, P, NHARM),
                        resid=b.hred.values - model)
    # odd-harmonic F-test at 2P
    P2 = res["2x"]["P"]
    c_e, _, chi2_e, np_e = linfit(b, P2, NHARM, odd=False)
    chi2_f, np_f = res["2x"]["chi2"], res["2x"]["npar"]
    dof_f = len(b) - np_f - 1
    F = ((chi2_e - chi2_f) / (np_f - np_e)) / (chi2_f / dof_f)
    p_odd = float(stats.f.sf(F, np_f - np_e, dof_f))
    beta = red_beta(b.t.values, res["2x"]["resid"], b.emag.values)
    # F-test with red noise: inflate chi2 by beta^2
    F_r = F / beta ** 2
    p_odd_red = float(stats.f.sf(F_r, np_f - np_e, dof_f))
    out.update(P1_h=round(res["1x"]["P"] * 24, 5), amp1=round(res["1x"]["amp"], 3),
               P2_h=round(P2 * 24, 5), sP2_h=round(res["2x"]["sP"] * 24 * beta, 5),
               amp2=round(res["2x"]["amp"], 3), rchi2_2=round(res["2x"]["rchi2"], 2), beta=round(beta, 2),
               beta_alpha_mag_per_deg=round(float(res["2x"]["c"][1]), 4),
               F_odd=round(float(F), 2), p_odd=p_odd, p_odd_red=p_odd_red)
    out["auto_choice"] = "2P_LS (odd harmonics significant)" if p_odd_red < 1e-3 else "ambiguous (P_LS or 2P_LS)"
    if P_known_h:
        out["P_known_h"] = P_known_h
    if plot:
        make_plot(sid, sector, df, keep, b, freq, pw, res, out, P_known_h, label)
    return out


def make_plot(sid, sector, df, keep, b, freq, pw, res, out, P_known_h, label):
    PLOTS.mkdir(exist_ok=True)
    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.1])
    ax = fig.add_subplot(gs[0, :])
    fm = out["median_flux"]
    ok = keep.values
    ax.plot(df.t_lt[ok], df.flux[ok] / fm, ",", color="0.6", alpha=0.5)
    ax.plot(b.t, 10 ** (-0.4 * (b.hred - np.median(b.hred))), ".", ms=3, color="k")
    lo, hi = np.nanpercentile(df.flux[ok] / fm, [1, 99])
    ax.set_ylim(min(lo, 0.3), max(hi, 1.7))
    bad = df.contam & (df.flag == 0)
    ax.plot(df.t_lt[bad], np.full(bad.sum(), ax.get_ylim()[0] + 0.02), "|", color="r", ms=4)
    ax.set_xlabel("BTJD (light-time corrected)")
    ax.set_ylabel("rel. flux")
    ax.set_title(f"{label or sid}  S{sector}  V~{out['V_med']}  kept {out['n_keep']}/{out['n_frames']} frames "
                 f"(red ticks: star-contaminated)  bins {BIN_MIN:.0f} min: scatter {out.get('bin_scatter_mag')} mag",
                 fontsize=9)
    ax = fig.add_subplot(gs[1, :])
    ax.plot(freq, pw, lw=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("frequency (c/d)")
    ax.set_ylabel("LS power")
    ax.axvline(1 / (out["P_LS_h"] / 24), color="r", ls=":", lw=0.8)
    if P_known_h:
        for m in (1, 2):
            ax.axvline(m * 24 / P_known_h, color="g", ls="--", lw=0.8)
    ax.set_title(f"P_LS = {out['P_LS_h']:.4f} h (FAP {out['LS_FAP']:.1e}); green: known P and P/2" if P_known_h
                 else f"P_LS = {out['P_LS_h']:.4f} h (FAP {out['LS_FAP']:.1e})", fontsize=9)
    for j, tag in enumerate(["1x", "2x"]):
        r = res[tag]
        ax = fig.add_subplot(gs[2, j])
        P = r["P"]
        ph = ((b.t - b.t.min()) / P) % 1
        yd = b.hred - r["c"][0] - r["c"][1] * (b.alpha - b.alpha.mean())
        ax.errorbar(ph, yd, b.emag, fmt=".", ms=2, color="k", ecolor="0.8", lw=0.5)
        ax.errorbar(ph + 1, yd, b.emag, fmt=".", ms=2, color="0.5", ecolor="0.9", lw=0.5)
        xx = np.linspace(0, 2, 500)
        yy = np.zeros_like(xx)
        for k in range(1, NHARM + 1):
            yy += r["c"][2 * k] * np.cos(2 * np.pi * k * xx) + r["c"][2 * k + 1] * np.sin(2 * np.pi * k * xx)
        ax.plot(xx, yy, "r-", lw=1)
        ax.invert_yaxis()
        lo, hi = np.percentile(yd, [0.5, 99.5])
        pad = 0.2 * (hi - lo)
        ax.set_ylim(hi + pad, lo - pad)
        ax.set_xlabel("phase")
        ax.set_ylabel("reduced mag - trend")
        extra = f"; odd-harm p={out['p_odd_red']:.1e} (red-noise scaled)" if tag == "2x" else ""
        ax.set_title(f"fold at {P*24:.4f} h, amp {r['amp']:.3f} mag{extra}", fontsize=9)
    fig.tight_layout()
    fig.savefig(PLOTS / f"{sid}_s{sector:03d}.png", dpi=80)
    plt.close(fig)


if __name__ == "__main__":
    sid, sector = sys.argv[1], int(sys.argv[2])
    Pk = float(sys.argv[3]) if len(sys.argv) > 3 else None
    o = analyse(sid, sector, Pk)
    for k, v in o.items():
        print(f"{k}: {v}")
