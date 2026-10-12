"""Deep K2-C0 EBs with VSX periods: EPIC 202072917 (TIC 308255000), 202071902 (TIC 57051340), 202126877 (TIC 86491361).
Independent check with the TESS light curves in data/raw (best product per sector, as the pipeline uses):
  1. find every eclipse-like dip (sliding box, depth > 0.25 x catalogue depth, SNR > 10) and time it with a trapezoid fit;
  2. test the VSX period: predicted eclipses in TESS data without a dip;
  3. test the pipeline period: residuals, odd/even depths, search for a secondary eclipse at any phase;
  4. for 202126877 also test the 'pair 7.25 d apart repeating every 21.78 d' hypothesis directly.
Writes vsx_ebs.txt and vsx_ebs_<TIC>.png."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import alias_pipeline as ap  # noqa: E402

TARGETS = [  # EPIC, TIC, K2 epoch (BTJD), depth, dur_h, VSX P, pipeline P(s)
    (202072917, 308255000, -213.3739, 0.248395, 7.3536, 3.842372, [34.746]),
    (202071902, 57051340, -208.4293, 0.091224, 5.3952, 8.69725, [22.904]),
    (202126877, 86491361, -211.5659, 0.061834, 5.8848, 3.877092, [21.775, 21.952, 43.550]),
]
out = open(HERE / "vsx_ebs.txt", "w")


def say(*a):
    print(*a); print(*a, file=out)


def trap(t, t0, depth, T14, Tin):
    Tin = np.clip(Tin, 1e-3, T14 / 2)
    x = np.abs(t - t0); m = np.zeros_like(t)
    full = x <= T14 / 2 - Tin; ramp = (x > T14 / 2 - Tin) & (x < T14 / 2)
    m[full] = depth; m[ramp] = depth * (T14 / 2 - x[ramp]) / Tin
    return 1 - m


def fit_dip(t, f, t_seed, D):
    sel = np.abs(t - t_seed) < D / 2 + 0.6
    tt, ff = t[sel], f[sel]
    if len(tt) < 15:
        return None
    tm = tt.mean()

    def model(p, tt):
        return trap(tt, *p[:4]) * (p[4] + p[5] * (tt - tm))

    p0 = [t_seed, 0.1, D, 0.2 * D, 1.0, 0.0]
    lo = [t_seed - D / 2, 0, 0.3 * D, 0.005, 0.5, -1]; hi = [t_seed + D / 2, 1, 2.5 * D, 1.25 * D, 1.5, 1]
    try:
        r = least_squares(lambda p: ff - model(p, tt), p0, bounds=(lo, hi), x_scale="jac")
    except Exception:
        return None
    res = ff - model(r.x, tt)
    J = r.jac; cov = np.linalg.pinv(J.T @ J) * np.var(res) * len(res) / max(1, len(res) - 6)
    return dict(t0=r.x[0], t0_err=np.sqrt(cov[0, 0]), depth=r.x[1], depth_err=np.sqrt(cov[1, 1]), T14=r.x[2], Tin=r.x[3], n=len(tt))


def find_dips(lc, D, delta):
    """Sliding box (width 0.7 D) with local linear baseline; returns clustered dip seeds (time, depth, snr, sector)."""
    t, f, sec = lc.time.to_numpy(), lc.flux.to_numpy(), lc.sector.to_numpy()
    hits = []
    for T in np.arange(t.min(), t.max(), D / 8):
        inn = np.abs(t - T) < 0.35 * D
        oot = (np.abs(t - T) > 0.5 * D + 0.05) & (np.abs(t - T) < 0.5 * D + 0.8)
        if inn.sum() < 5 or (oot & (t < T)).sum() < 8 or (oot & (t > T)).sum() < 8:
            continue
        c = np.polyfit(t[oot] - T, f[oot], 1)
        r = f / np.polyval(c, t - T) - 1
        sig = 1.4826 * np.median(np.abs(r[oot] - np.median(r[oot]))) + 1e-6
        d = -np.mean(r[inn]); s = sig / np.sqrt(inn.sum())
        if d > 0.25 * delta and d / s > 10:
            hits.append((T, d, d / s, int(sec[np.searchsorted(t, T) - 1])))
    # cluster within D
    hits.sort()
    cl = []
    for h in hits:
        if cl and h[0] - cl[-1][-1][0] < D:
            cl[-1].append(h)
        else:
            cl.append([h])
    return [max(c, key=lambda x: x[1]) for c in cl]


def coverage(lc, Tc, D):
    """Is a full eclipse at Tc covered by data (>=60% of the +-0.35 D core and baseline both sides)?"""
    t = lc.time.to_numpy()
    cad = np.median(np.diff(t[np.abs(t - Tc) < 3])) if (np.abs(t - Tc) < 3).sum() > 5 else np.nan
    inn = np.abs(t - Tc) < 0.35 * D
    oot_l = (t < Tc - 0.5 * D - 0.05) & (t > Tc - 0.5 * D - 0.8); oot_r = (t > Tc + 0.5 * D + 0.05) & (t < Tc + 0.5 * D + 0.8)
    if not np.isfinite(cad) or inn.sum() * cad < 0.6 * 0.7 * D or oot_l.sum() < 8 or oot_r.sum() < 8:
        return None
    f = lc.flux.to_numpy()
    c = np.polyfit(t[oot_l | oot_r] - Tc, f[oot_l | oot_r], 1)
    r = f / np.polyval(c, t - Tc) - 1
    return -np.mean(r[inn]), 1.4826 * np.median(np.abs(r[oot_l | oot_r])) / np.sqrt(inn.sum())


for epic, tic, epoch, delta, dur_h, p_vsx, p_pipe in TARGETS:
    D = dur_h / 24
    lc = ap.load_lightcurves(tic, products=ap.PRODUCT_RANK + ["OWN"])
    say(f"\n================ EPIC {epic} = TIC {tic}: K2 C0 event {epoch}, depth {delta*1e3:.0f} ppt, D {dur_h:.1f} h")
    say("TESS sectors/products:", {int(s): str(lc[lc.sector == s]["product"].iloc[0]) for s in sorted(lc.sector.unique())})
    dips = find_dips(lc, D, delta)
    fits = []
    for T, d, snr, sec in dips:
        r = fit_dip(lc.time.to_numpy(), lc.flux.to_numpy(), T, D)
        if r:
            fits.append(dict(sector=sec, **{k: r[k] for k in ("t0", "t0_err", "depth", "depth_err", "T14", "Tin")}, box_depth=d, box_snr=snr))
    df = pd.DataFrame(fits)
    say(f"Dips found in TESS ({len(df)}):")
    if len(df):
        say(df.assign(T14_h=df.T14 * 24, Tin_h=df.Tin * 24).drop(columns=["T14", "Tin"]).round(4).to_string(index=False))
    # VSX period test
    say(f"\nVSX period {p_vsx} d: predicted eclipses (from the K2 epoch) that fall in covered TESS data:")
    t_all = lc.time.to_numpy()
    n_cov = n_flat = 0
    for k in range(int((t_all.min() - epoch) / p_vsx) - 1, int((t_all.max() - epoch) / p_vsx) + 2):
        Tc = epoch + k * p_vsx
        if Tc < t_all.min() or Tc > t_all.max():
            continue
        cv = coverage(lc, Tc, D)
        if cv is None:
            continue
        n_cov += 1
        if cv[0] < 0.25 * delta:
            n_flat += 1
    say(f"   {n_cov} predicted eclipses covered; {n_flat} show NO eclipse (depth < 25% of the K2 depth). "
        f"Note: the K2 epoch is {abs(t_all.min()-epoch)/p_vsx:.0f} cycles away, so phase drift is possible; the "
        f"decisive test is the number of eclipses per sector: expected ~{27/p_vsx:.0f} per sector at P={p_vsx}, found {len(df)} in {lc.sector.nunique()} sectors.")
    # pipeline period(s)
    for P in p_pipe:
        say(f"\nPipeline period {P} d, linear ephemeris from the K2 epoch:")
        if len(df) == 0:
            continue
        cyc = np.round((df.t0 - epoch) / P)
        oc = df.t0 - (epoch + cyc * P)
        # refit P, T0 with all dips whose |O-C| < D/3
        use = np.abs(oc) < D / 3
        if use.sum() >= 1:
            A = np.vstack([np.ones(use.sum() + 1), np.concatenate([[0], cyc[use]])]).T
            y = np.concatenate([[epoch], df.t0[use]])
            T0f, Pf = np.linalg.lstsq(A, y, rcond=None)[0]
            cyc2 = np.round((df.t0 - T0f) / Pf); oc2 = df.t0 - (T0f + cyc2 * Pf)
            say(f"   refit: P = {Pf:.5f} d, T0 = {T0f:.4f}; dips on ephemeris: {int((np.abs(oc2) < D/3).sum())}/{len(df)}")
            for i, r in df.iterrows():
                say(f"     t0 {r.t0:.4f} (S{r.sector}) cycle {int(cyc2[i])} O-C {oc2[i]*24:+.2f} h depth {r.depth*1e3:.1f} ppt T14 {r.T14*24:.1f} h"
                    f" {'ON' if abs(oc2[i]) < D/3 else 'OFF'}")
            on = df[np.abs(oc2) < D / 3]
            if len(on) >= 2:
                odd = on[cyc2[on.index] % 2 == 1]; even = on[cyc2[on.index] % 2 == 0]
                say(f"   odd/even: even depths {np.round(even.depth.values*1e3,1)}, odd depths {np.round(odd.depth.values*1e3,1)} ppt;"
                    f" even T14 {np.round(even.T14.values*24,1)} h, odd T14 {np.round(odd.T14.values*24,1)} h")
            # predicted eclipses at Pf that are covered but missing
            miss = []
            for k in range(int((t_all.min() - T0f) / Pf) - 1, int((t_all.max() - T0f) / Pf) + 2):
                Tc = T0f + k * Pf
                if Tc < t_all.min() or Tc > t_all.max():
                    continue
                cv = coverage(lc, Tc, D)
                if cv is not None and cv[0] < 0.25 * delta:
                    miss.append((k, round(Tc, 2), round(cv[0] * 1e3, 1)))
            say(f"   covered predicted eclipses WITHOUT a dip: {miss if miss else 'none'}")
            # secondary search: fold all data at Pf, scan phase for a dip of width 0.7 D, excluding +-D of phase 0
            t, f, sec = lc.time.to_numpy(), lc.flux.to_numpy(), lc.sector.to_numpy()
            fd = np.empty_like(f)
            for s in np.unique(sec):
                m = sec == s
                tt, ff = t[m], f[m]
                inpri = np.abs(((tt - T0f) / Pf + 0.5) % 1 - 0.5) * Pf < D
                base = np.full_like(ff, np.nan)
                for i in range(len(tt)):
                    w = (np.abs(tt - tt[i]) < 1.0) & ~inpri
                    if w.sum() > 10:
                        base[i] = np.median(ff[w])
                fd[m] = ff / base
            ok = np.isfinite(fd)
            ph = ((t[ok] - T0f) / Pf) % 1
            best = None
            scan = []
            for pc in np.arange(0, 1, 0.002):
                dph = np.abs(((ph - pc) + 0.5) % 1 - 0.5) * Pf
                if pc * Pf < D or (1 - pc) * Pf < D:
                    continue
                inn = dph < 0.35 * D
                if inn.sum() < 10:
                    continue
                d = 1 - np.mean(fd[ok][inn]); s = np.std(fd[ok][~inn & (dph < 2)]) / np.sqrt(inn.sum()) if (~inn & (dph < 2)).sum() > 10 else np.nan
                scan.append((pc, d, s, inn.sum()))
                if best is None or d > best[1]:
                    best = (pc, d, s, inn.sum())
            if best:
                say(f"   secondary-eclipse scan at P={Pf:.4f}: deepest phase {best[0]:.3f} (dt={best[0]*Pf:.2f} d after primary): depth {best[1]*1e3:.2f} +- {best[2]*1e3:.2f} ppt (n={best[3]});"
                    f" phase-0.5 depth: {[round(x[1]*1e3,2) for x in scan if abs(x[0]-0.5) < 0.0015]} ppt")
                # at P/2 would phase 0.5 need to be a primary of equal depth -> report
            fig, axs = plt.subplots(2, 1, figsize=(11, 7))
            axs[0].plot(t, f, ".", ms=1.5, color="0.4")
            for _, r in df.iterrows():
                axs[0].axvline(r.t0, color="r", lw=0.5)
            for k in range(int((t_all.min() - T0f) / Pf) - 1, int((t_all.max() - T0f) / Pf) + 2):
                axs[0].axvline(T0f + k * Pf, color="b", lw=0.5, ls="--")
            axs[0].set_title(f"TIC {tic}: red = dips found, blue dashed = P={Pf:.4f} ephemeris")
            axs[1].plot(ph, fd[ok], ".", ms=1.5); axs[1].set_xlabel(f"phase (P={Pf:.4f})"); axs[1].set_ylim(1 - 1.5 * delta, 1.02)
            fig.tight_layout(); fig.savefig(HERE / f"vsx_ebs_{tic}_{P}.png", dpi=90); plt.close(fig)
    if epic == 202126877 and len(df) >= 2:
        say("\nPair structure test (202126877): all dip-time differences < 30 d:")
        tt = np.sort(df.t0.values)
        for i in range(len(tt)):
            for j in range(i + 1, len(tt)):
                if tt[j] - tt[i] < 30:
                    say(f"   {tt[i]:.3f} -> {tt[j]:.3f}: {tt[j]-tt[i]:.3f} d")
out.close()
