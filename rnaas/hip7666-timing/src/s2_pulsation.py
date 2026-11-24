"""Pulsation phases and time delays (RULES.md section 4).

Needs results/prepared_main.pkl and results/template_S58.pkl from s1_eclipses.py.
Writes results/pulsation_delays.csv (one row per orbit-long segment and mode, plus a
common-delay row "all"), results/pulsation_orbit_check.csv and results/pulsation_summary.txt.
"""
import csv
import pickle

import numpy as np
from scipy.optimize import least_squares

import hip

P = pickle.load(open(hip.OUTDIR / "prepared_main.pkl", "rb"))
TPL = pickle.load(open(hip.OUTDIR / "template_S58.pkl", "rb"))
NAMES = list(hip.FREQS)
FN = np.array([hip.FREQS[n] for n in NAMES])
MODES = hip.TIMING_MODES
A_LTT_INNER = K_A_S = hip.K_A * hip.P_LEE * 86400 / (2 * np.pi * hip.C_KMS)   # a_A sin i / c in s


def residual(d):
    """Out-of-eclipse pulsation signal: f/trend - binary template (RULES 2)."""
    shape = hip.smear(lambda tt: TPL(hip.phase(tt)), d["t"], d["exptime"])
    return d["f"] / d["trend"] - shape


def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def per_mode(t, r, coef, modes):
    """Fit amplitude+phase of `modes` with the other frequencies subtracted at sector values."""
    others = [i for i, n in enumerate(NAMES) if n not in modes]
    sub = hip.design(t, FN[others]) @ np.concatenate([coef[2 * i:2 * i + 2] for i in others])
    fm = [hip.FREQS[m] for m in modes]
    c, cov, rms = hip.fit_sines(t, r - sub, fm)
    out = {}
    for j, m in enumerate(modes):
        a, b = c[2 * j], c[2 * j + 1]
        A, ph = hip.amp_phase(a, b)
        va, vb, cab = cov[2 * j, 2 * j], cov[2 * j + 1, 2 * j + 1], cov[2 * j, 2 * j + 1]
        sph = np.sqrt(max(b ** 2 * va + a ** 2 * vb - 2 * a * b * cab, 0)) / A ** 2
        out[m] = (A, ph, sph)
    return out


def common_tau(t, r, coef):
    """Fit r = k * M(t - tau) + c with the whole sector pulsation model M (all 14 modes)."""
    def model(p):
        return p[1] * hip.sines_model(t - p[0], FN, coef) + p[2]
    fun = lambda p: model(p) - r
    sol = least_squares(fun, [0.0, 1.0, 0.0], x_scale=[1e-4, 0.1, 1e-4])
    J = sol.jac
    s2 = np.sum(sol.fun ** 2) / max(len(t) - 3, 1)
    cov = np.linalg.pinv(J.T @ J) * s2
    return sol.x[0], np.sqrt(cov[0, 0]), sol.x[1]


def segments(t, length=None):
    """Orbit-long segments from phase -0.25 to 0.75 (default), else fixed-length in time."""
    if length is None:
        n = np.floor((t - hip.T0_LEE) / hip.P_LEE + 0.25).astype(int)
    else:
        n = np.floor((t - t.min()) / length).astype(int)
    return [np.where(n == k)[0] for k in np.unique(n)]


def main():
    rows, summary = [], []
    for s, d in P.items():
        t, r, coef, ex = d["t"], residual(d), d["coef"], d["exptime"]
        oot = ~hip.in_eclipse(t)
        sector_ph = per_mode(t[oot], r[oot], coef, MODES)
        nexp = hip.P_LEE * (1 - 4 * hip.ECL_HALF) * 86400 / ex
        for idx in segments(t):
            idx = idx[oot[idx]]
            if len(idx) < 0.6 * nexp:
                continue
            tm = np.mean(t[idx])
            seg = per_mode(t[idx], r[idx], coef, MODES)
            for m in MODES:
                A, ph, sph = seg[m]
                dphi = wrap(ph - sector_ph[m][1])
                f = hip.FREQS[m]
                rows.append([s, m, f"{tm:.4f}", len(idx), f"{A*1e3:.4f}",
                             f"{-dphi / (2*np.pi*f) * 86400:.2f}", f"{sph / (2*np.pi*f) * 86400:.2f}"])
            tau, stau, k = common_tau(t[idx], r[idx], coef)
            rows.append([s, "all", f"{tm:.4f}", len(idx), f"{k:.3f}", f"{tau*86400:.2f}", f"{stau*86400:.2f}"])
    with open(hip.OUTDIR / "pulsation_delays.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sector", "mode", "t_mid", "n", "amp_ppt_or_scale", "tau_s", "err_s"])
        w.writerows(rows)

    # within-sector: remove a linear term per sector and mode, then chi2 of the rest (RULES 4)
    import pandas as pd
    df = pd.read_csv(hip.OUTDIR / "pulsation_delays.csv")
    df["tau_detr"] = np.nan
    for (s, m), g in df.groupby(["sector", "mode"]):
        if len(g) < 3:
            continue
        w8 = 1 / g.err_s ** 2
        X = np.column_stack([np.ones(len(g)), g.t_mid - g.t_mid.mean()])
        c = np.linalg.solve((X * w8.values[:, None]).T @ X, (X * w8.values[:, None]).T @ g.tau_s.values)
        res = g.tau_s.values - X @ c
        df.loc[g.index, "tau_detr"] = res
        chi2 = np.sum(res ** 2 * w8)
        summary.append(f"S{s} {m}: {len(g)} segments, median err {g.err_s.median():.1f} s, "
                       f"slope {c[1]:+.2f} s/d, chi2 about line {chi2:.1f} for {len(g)-2} dof, "
                       f"rms {np.std(res):.1f} s")
    df.to_csv(hip.OUTDIR / "pulsation_delays.csv", index=False)

    # S58-S85 phase connection test (RULES 4): chi2 vs frequency with common amplitude+phase
    summary.append("")
    summary.append("S58-S85 phase connection (common amplitude and phase, other modes subtracted per sector):")
    data = {}
    for s in (58, 85):
        d = P[s]
        t, r, coef = d["t"], residual(d), d["coef"]
        oot = ~hip.in_eclipse(t)
        data[s] = (t[oot], r[oot], coef)
    gap = np.mean(data[85][0]) - np.mean(data[58][0])
    for m in MODES:
        f0 = hip.FREQS[m]
        grid = np.arange(f0 - 0.003, f0 + 0.003, 2e-6)
        tt, yy = [], []
        for s in (58, 85):
            t, r, coef = data[s]
            others = [i for i, n in enumerate(NAMES) if n != m]
            sub = hip.design(t, FN[others]) @ np.concatenate([coef[2*i:2*i+2] for i in others])
            tt.append(t); yy.append(r - sub - np.mean(r - sub))
        tt, yy = np.concatenate(tt), np.concatenate(yy)
        chi = np.array([np.sum((yy - hip.design(tt, [f]) @ np.linalg.lstsq(hip.design(tt, [f]), yy, rcond=None)[0])**2)
                        for f in grid])
        s2 = chi.min() / (len(tt) - 3)
        chi = chi / s2
        # local minima
        mins = np.where((chi[1:-1] < chi[:-2]) & (chi[1:-1] < chi[2:]))[0] + 1
        order = mins[np.argsort(chi[mins])]
        best, second = order[0], order[1]
        summary.append(f"  {m}: best f {grid[best]:.6f} (Lee {f0}), next alias {grid[second]:.6f}, "
                       f"dchi2 {chi[second]-chi[best]:.1f}; alias spacing 1/gap = {1/gap:.6f}")
    summary.append("  Note: with two sectors and a free frequency, a constant delay between them is "
                   "absorbed by the frequency, so even a unique alias gives no delay measurement.")

    # validation: inner-orbit light time of the pulsator in S58 (0.3-d segments, common tau)
    d = P[58]
    t, r, coef = d["t"], residual(d), d["coef"]
    oot = ~hip.in_eclipse(t)
    chk = []
    for idx in segments(t, length=0.3):
        idx = idx[oot[idx]]
        if len(idx) < 0.6 * 0.3 * 86400 / d["exptime"]:
            continue
        tau, stau, k = common_tau(t[idx], r[idx], coef)
        chk.append([np.mean(t[idx]), hip.phase(np.mean(t[idx])) % 1.0, tau * 86400, stau * 86400])
    chk = np.array(chk)
    np.savetxt(hip.OUTDIR / "pulsation_orbit_check.csv", chk, delimiter=",",
               header="t_mid,orb_phase,tau_s,err_s", comments="", fmt="%.5f")
    ph, y, e = chk[:, 1], chk[:, 2], chk[:, 3]
    X = np.column_stack([np.ones_like(ph), np.cos(2*np.pi*ph), np.sin(2*np.pi*ph),
                         chk[:, 0] - chk[:, 0].mean()])
    W = 1 / e**2
    cov = np.linalg.inv((X * W[:, None]).T @ X)
    c = cov @ (X * W[:, None]).T @ y
    chi2 = np.sum((y - X @ c)**2 * W)
    sc = np.sqrt(max(chi2 / (len(y) - 4), 1))
    summary.append("")
    summary.append(f"Inner-orbit check (S58, {len(y)} segments of 0.3 d, common tau, linear term free):")
    summary.append(f"  tau = c0 + a cos(2 pi phi) + b sin(2 pi phi): a = {c[1]:+.1f} +/- {np.sqrt(cov[1,1])*sc:.1f} s, "
                   f"b = {c[2]:+.1f} +/- {np.sqrt(cov[2,2])*sc:.1f} s (errors scaled by {sc:.2f})")
    summary.append(f"  prediction if the primary pulsates: a = -{A_LTT_INNER:.1f} s, b = 0 "
                   f"(secondary: a = +{hip.K_B*hip.P_LEE*86400/(2*np.pi*hip.C_KMS):.1f} s)")
    txt = "\n".join(summary)
    print(txt)
    (hip.OUTDIR / "pulsation_summary.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()
