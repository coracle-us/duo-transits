"""Eclipse O-C models, BIC choice and the LTTE / tertiary-mass limit (RULES.md sections 3, 5, 6).

Usage: python src/s3_models.py [main|qlp]
Writes results/oc_<kind>.csv, results/models_<kind>.txt and (main only) results/ltte_limit_<kind>.csv.
"""
import sys

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import chi2 as chi2dist

import hip

KIND = sys.argv[1] if len(sys.argv) > 1 else "main"
G, MSUN, CLIGHT, DAY = 6.674e-11, 1.989e30, 2.998e8, 86400.0
MJUP = 1.898e27
out = []


def say(s=""):
    print(s)
    out.append(s)


df = pd.read_csv(hip.OUTDIR / f"eclipse_times_{KIND}.csv")
df["cyc"] = df.E + np.where(df.type == "s", 0.5, 0.0)
df["sec"] = (df.type == "s").astype(float)
df["err_raw"] = df.err_d

# ---- error scaling per sector and type (RULES 3)
say(f"Eclipse times ({KIND}): {len(df)} used; per sector and type, scatter about a line:")
for (s, ty), g in df.groupby(["sector", "type"]):
    w = 1 / g.err_d ** 2
    X = np.column_stack([np.ones(len(g)), g.cyc - g.cyc.mean()])
    c = np.linalg.solve((X * w.values[:, None]).T @ X, (X * w.values[:, None]).T @ g.t_btjd.values)
    res = g.t_btjd.values - X @ c
    chi = np.sum(res ** 2 * w)
    dof = len(g) - 2
    fac = np.sqrt(chi / dof) if dof > 0 and chi / dof > 1 else 1.0
    df.loc[g.index, "err_d"] = g.err_d * fac
    say(f"  S{s} {ty}: n={len(g)}, median boot err {g.err_raw.median()*DAY:.1f} s, "
        f"rms about line {np.std(res)*DAY:.1f} s, chi2/dof {chi:.1f}/{dof}, error scale {fac:.2f}")

t, e, cyc, sec = df.t_btjd.values, df.err_d.values, df.cyc.values, df.sec.values
N = len(t)
W = 1 / e ** 2
c0 = np.round(np.average(cyc, weights=W))


def wls(X, y):
    A = (X * W[:, None]).T @ X
    cov = np.linalg.inv(A)
    beta = cov @ (X * W[:, None]).T @ y
    r = y - X @ beta
    return beta, cov, r, float(np.sum(r ** 2 * W))


x = cyc - c0
X0 = np.column_stack([np.ones(N), x, sec])
X1 = np.column_stack([X0, x ** 2])
b0, cov0, r0, chi0 = wls(X0, t)
b1, cov1, r1, chi1 = wls(X1, t)
bic0, bic1 = chi0 + 3 * np.log(N), chi1 + 4 * np.log(N)

say()
say(f"M0 linear: T0 = {b0[0]:.6f} +/- {np.sqrt(cov0[0,0]):.6f} BTJD (cycle {c0:.0f}), "
    f"P = {b0[1]:.8f} +/- {np.sqrt(cov0[1,1]):.8f} d, secondary offset = {b0[2]*DAY:+.1f} +/- "
    f"{np.sqrt(cov0[2,2])*DAY:.1f} s; chi2 = {chi0:.1f} for {N-3} dof, BIC {bic0:.1f}")
say(f"   P(TESS) - P(Lee Eq.1) = {(b0[1]-hip.P_LEE)*DAY:+.3f} +/- {np.hypot(np.sqrt(cov0[1,1]), 1.8e-7)*DAY:.3f} s")
say(f"M1 quadratic: Q = {b1[3]:.3e} +/- {np.sqrt(cov1[3,3]):.3e} d/cycle^2 "
    f"(dP/dE = 2Q = {2*b1[3]:.3e}); chi2 = {chi1:.1f}, BIC {bic1:.1f}, dBIC vs M0 {bic1-bic0:+.1f}")

# M2: circular LTTE on a P3 grid
P3 = np.geomspace(10, 3000, 400)
best = (np.inf, None)
for p3 in P3:
    ph = 2 * np.pi * t / p3
    X2 = np.column_stack([X0, np.cos(ph), np.sin(ph)])
    try:
        b2, cov2, r2, chi2v = wls(X2, t)
    except np.linalg.LinAlgError:
        continue
    if chi2v < best[0]:
        best = (chi2v, p3, b2, cov2)
chi2b, p3b, b2, cov2 = best
bic2 = chi2b + 6 * np.log(N)
A2 = np.hypot(b2[3], b2[4]) * DAY
say(f"M2 LTTE (best on grid): P3 = {p3b:.1f} d, A = {A2:.1f} s; chi2 = {chi2b:.1f}, BIC {bic2:.1f}, "
    f"dBIC vs M0 {bic2-bic0:+.1f}")
pref = min((bic0, "M0"), (bic1, "M1"), (bic2, "M2"))[1]
sig = "M0" if min(bic1, bic2) - bic0 > -10 else pref
say(f"Preferred by BIC with the dBIC < -10 rule: {sig}")

# primary and secondary separately
say()
for ty in ("p", "s"):
    m = df.type == ty
    Xs = np.column_stack([np.ones(m.sum()), x[m]])
    Ws = W[m]
    A = (Xs * Ws[:, None]).T @ Xs
    cv = np.linalg.inv(A)
    bb = cv @ (Xs * Ws[:, None]).T @ t[m]
    rr = t[m] - Xs @ bb
    say(f"{ty} only: P = {bb[1]:.8f} +/- {np.sqrt(cv[1,1]):.8f} d, chi2 {np.sum(rr**2*Ws):.1f}/{m.sum()-2}")

# O-C table and per-sector summaries
df["oc_lee_s"] = (t - (hip.T0_LEE + cyc * hip.P_LEE)) * DAY
df["oc_own_s"] = (t - (b0[0] + x * b0[1] + sec * b0[2])) * DAY
df["err_s"] = e * DAY
df.to_csv(hip.OUTDIR / f"oc_{KIND}.csv", index=False)
say()
say("Per sector, weighted mean O-C (s) [against Lee Eq. 1 | against own M0 linear]:")
for s, g in df.groupby("sector"):
    for ty in ("p", "s"):
        h = g[g.type == ty]
        w = 1 / h.err_s ** 2
        say(f"  S{s} {ty}: n={len(h)}, Lee {np.average(h.oc_lee_s, weights=w):+.1f} +/- {1/np.sqrt(w.sum()):.1f}, "
            f"own {np.average(h.oc_own_s, weights=w):+.1f} +/- {1/np.sqrt(w.sum()):.1f}")
    hp, hs = g[g.type == "p"], g[g.type == "s"]
    wp, ws = 1 / hp.err_s ** 2, 1 / hs.err_s ** 2
    d = np.average(hs.oc_lee_s, weights=ws) - np.average(hp.oc_lee_s, weights=wp)
    say(f"  S{s} secondary - primary: {d:+.1f} +/- {np.sqrt(1/wp.sum()+1/ws.sum()):.1f} s")
romer = hip.P_LEE * DAY * (hip.K_B - hip.K_A) / (np.pi * hip.C_KMS)
say(f"  inner-orbit Romer delay predicted (Kaplan 2010): {romer:+.1f} s")

if KIND in ("main", "mixq18"):
    # ---- injection-recovery limit (RULES 6)
    thr = chi2dist.ppf(0.9973, N - 3)
    R = np.eye(N) - X0 @ np.linalg.inv((X0 * W[:, None]).T @ X0) @ (X0 * W[:, None]).T
    rng = np.random.default_rng(7)
    phis = rng.uniform(0, 2 * np.pi, 1000)
    Agrid = np.geomspace(1, 20000, 300) / DAY
    rows = []
    for p3 in np.geomspace(10, 3000, 600):
        rc, rs = R @ np.cos(2 * np.pi * t / p3), R @ np.sin(2 * np.pi * t / p3)
        A95 = np.nan
        for A in Agrid:
            res = r0[None, :] + A * (np.cos(phis)[:, None] * rc[None, :] + np.sin(phis)[:, None] * rs[None, :])
            chi = np.sum(res ** 2 * W[None, :], axis=1)
            if np.mean(chi > thr) >= 0.95:
                A95 = A * DAY
                break
        if np.isfinite(A95):
            fm = 4 * np.pi ** 2 * (A95 * CLIGHT) ** 3 / (G * (p3 * DAY) ** 2) / MSUN
            m3 = brentq(lambda m: m ** 3 / (hip.M12 + m) ** 2 - fm, 1e-8, 1e4)
        else:
            m3 = np.nan
        rows.append([p3, A95, m3])
    lim = pd.DataFrame(rows, columns=["P3_d", "A95_s", "m3min_msun"])
    lim.to_csv(hip.OUTDIR / f"ltte_limit_{KIND}.csv", index=False)
    say()
    say(f"LTTE limit (95% of phases give M0 chi2 > {thr:.1f}, the 99.73% point for {N-3} dof; "
        f"observed M0 chi2 {chi0:.1f}):")
    for p in (10, 20, 30, 50, 100, 200, 300, 500, 700, 1000, 1500, 2000, 3000):
        rr = lim.iloc[(lim.P3_d - p).abs().argmin()]
        say(f"  P3 {rr.P3_d:7.1f} d: A < {rr.A95_s:8.1f} s, m3 sin i3 >~ excluded above "
            f"{rr.m3min_msun:7.3f} Msun ({rr.m3min_msun*MSUN/MJUP:7.1f} MJup)")

(hip.OUTDIR / f"models_{KIND}.txt").write_text("\n".join(out) + "\n")
