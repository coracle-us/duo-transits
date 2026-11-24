"""Step 3: M0 linear, M1 quadratic, M2 linear + circular LTTE (P3 grid 10-3000 d) on the
eclipse times (primary + secondary, free secondary offset). Also per-type periods,
per-sector secondary-primary offsets and O-C summaries.
Writes independent/out/models.txt, independent/out/m2_periodogram.csv, independent/out/oc.csv
"""
import os, sys
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
import hiplib as H

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
DAY_S = 86400.0


TIMES = os.environ.get("TIMES_CSV", os.path.join(HERE, "eclipse_times.csv"))
TAG = os.environ.get("OUTTAG", "")


def read_times(path=None):
    path = TIMES if path is None else path
    import csv
    rows = list(csv.DictReader(open(path)))
    sec = np.array([int(r["sector"]) for r in rows])
    typ = np.array([r["type"] for r in rows])
    cyc = np.array([float(r["cycle"]) for r in rows])
    t = np.array([float(r["t_btjd"]) for r in rows])
    e = np.array([float(r["err_d"]) for r in rows])
    e0 = np.array([float(r["err_d_unscaled"]) for r in rows])
    return sec, typ, cyc, t, e, e0


E_REF = 3420.0   # reference cycle (first S58 primary) to decorrelate T0 and P


def design(cyc, is_sec, t, model="M0", P3=None):
    x = cyc - E_REF
    cols = [np.ones_like(x), x, is_sec.astype(float)]
    names = ["T0", "P", "dsec"]
    if model == "M1":
        cols.append(x ** 2); names.append("Q")
    if model == "M2":
        arg = 2 * np.pi * t / P3
        cols += [np.sin(arg), np.cos(arg)]; names += ["a", "b"]
    return np.vstack(cols).T, names


def wls(X, y, e):
    A = X / e[:, None]
    b = y / e
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    cov = np.linalg.inv(A.T @ A)
    chi2 = float(np.sum((b - A @ c) ** 2))
    return c, cov, chi2


def p3_grid(tmin=10.0, tmax=3000.0, n=20000):
    f = np.linspace(1.0 / tmax, 1.0 / tmin, n)
    return 1.0 / f[::-1]


def main():
    sec, typ, cyc, t, e, e0 = read_times()
    is_sec = typ == "secondary"
    N = len(t)
    out = open(os.path.join(OUT, TAG + "models.txt"), "w")

    def say(*a):
        s = " ".join(str(x) for x in a); print(s); out.write(s + "\n")

    say("N eclipses = %d (primaries %d, secondaries %d)" % (N, (~is_sec).sum(), is_sec.sum()))
    for s in (18, 58, 85):
        say("  S%d: P=%d S=%d" % (s, np.sum((sec == s) & ~is_sec), np.sum((sec == s) & is_sec)))

    res = {}
    # M0
    X, names = design(cyc, is_sec, t, "M0")
    c0, cov0, chi0 = wls(X, t, e)
    k0 = 3
    res["M0"] = (chi0, N - k0, chi0 + k0 * np.log(N))
    sig0 = np.sqrt(np.diag(cov0))
    chi2nu0 = chi0 / (N - k0)
    say("\nM0 linear: chi2=%.3f dof=%d chi2nu=%.3f BIC=%.3f  p(chi2>obs)=%.4f" %
        (chi0, N - k0, chi2nu0, res["M0"][2], stats.chi2.sf(chi0, N - k0)))
    say("  T0(E=%.0f) = %.7f +- %.7f BTJD" % (E_REF, c0[0], sig0[0]))
    say("  P = %.9f +- %.9f d  (+- %.4f s); with sqrt(chi2nu) inflation: +- %.9f d" %
        (c0[1], sig0[1], sig0[1] * DAY_S, sig0[1] * max(1, np.sqrt(chi2nu0))))
    say("  P - P_Lee = %.3f +- %.3f s" % ((c0[1] - H.P_LEE) * DAY_S, sig0[1] * DAY_S))
    say("  secondary offset from phase 0.5 = %.2f +- %.2f s (inflated: +- %.2f s)" %
        (c0[2] * DAY_S, sig0[2] * DAY_S, sig0[2] * DAY_S * max(1, np.sqrt(chi2nu0))))
    # M1
    X, names = design(cyc, is_sec, t, "M1")
    c1, cov1, chi1 = wls(X, t, e)
    k1 = 4
    res["M1"] = (chi1, N - k1, chi1 + k1 * np.log(N))
    s1 = np.sqrt(np.diag(cov1))
    say("\nM1 quadratic: chi2=%.3f dof=%d BIC=%.3f dBIC(M1-M0)=%.3f" %
        (chi1, N - k1, res["M1"][2], res["M1"][2] - res["M0"][2]))
    say("  Q = %.4e +- %.4e d/cycle^2 ; dP/dE = 2Q = %.4e +- %.4e d/cycle ; dP/dt = %.4e +- %.4e (d/d)" %
        (c1[3], s1[3], 2 * c1[3], 2 * s1[3], 2 * c1[3] / c1[1], 2 * s1[3] / c1[1]))
    say("  P(E_ref) = %.8f +- %.8f ; dsec = %.2f +- %.2f s" % (c1[1], s1[1], c1[2] * DAY_S, s1[2] * DAY_S))
    # M2 grid
    grid = p3_grid()
    chis = np.empty(len(grid)); amps = np.empty(len(grid))
    for i, P3 in enumerate(grid):
        X, _ = design(cyc, is_sec, t, "M2", P3)
        c, cov, chi = wls(X, t, e)
        chis[i] = chi; amps[i] = np.hypot(c[3], c[4]) * DAY_S
    ib = int(np.argmin(chis))
    k2 = 6
    res["M2"] = (chis[ib], N - k2, chis[ib] + k2 * np.log(N))
    X, _ = design(cyc, is_sec, t, "M2", grid[ib])
    c2, cov2, _ = wls(X, t, e)
    say("\nM2 linear+circular LTTE: best P3=%.3f d, A=%.2f s, chi2=%.3f dof=%d BIC=%.3f dBIC(M2-M0)=%.3f" %
        (grid[ib], amps[ib], chis[ib], N - k2, res["M2"][2], res["M2"][2] - res["M0"][2]))
    say("  (k counts T0, P, dsec, a, b and P3 = 6; grid of %d P3 values uniform in frequency, 10-3000 d)" % len(grid))
    # next-best local minima
    order = np.argsort(chis)
    seen = []
    for j in order:
        if all(abs(np.log(grid[j] / g)) > 0.05 for g in seen):
            seen.append(grid[j])
            say("    local min: P3=%.2f d chi2=%.3f A=%.2f s" % (grid[j], chis[j], amps[j]))
        if len(seen) >= 6:
            break
    with open(os.path.join(OUT, TAG + "m2_periodogram.csv"), "w") as fh:
        fh.write("P3_d,chi2,A_s\n")
        for P3, ch, a in zip(grid, chis, amps):
            fh.write("%.5f,%.5f,%.4f\n" % (P3, ch, a))
    say("\nSummary table:")
    say("  model  chi2     dof  BIC      dBIC_vs_M0")
    for m in ("M0", "M1", "M2"):
        say("  %s  %8.3f  %3d  %8.3f  %+8.3f" % (m, res[m][0], res[m][1], res[m][2], res[m][2] - res["M0"][2]))
    pref = [m for m in ("M1", "M2") if res[m][2] - res["M0"][2] < -10]
    say("  preferred over M0 (dBIC < -10):", pref if pref else "none -> M0")

    # per-type periods
    say("\nPer-type linear fits (all sectors):")
    for lab, m in (("primary", ~is_sec), ("secondary", is_sec)):
        x = cyc[m] - E_REF
        X = np.vstack([np.ones_like(x), x]).T
        c, cov, chi = wls(X, t[m], e[m])
        s = np.sqrt(np.diag(cov)); nu = m.sum() - 2
        say("  %s: N=%d P=%.9f +- %.9f d (+- %.4f s) chi2=%.2f dof=%d; inflated err %.9f" %
            (lab, m.sum(), c[1], s[1], s[1] * DAY_S, chi, nu, s[1] * max(1, np.sqrt(chi / nu))))
    # unscaled-error M0 for reference
    X, _ = design(cyc, is_sec, t, "M0")
    cu, covu, chiu = wls(X, t, e0)
    say("\nM0 with unscaled errors: chi2=%.3f dof=%d P=%.8f +- %.8f dsec=%.2f +- %.2f s" %
        (chiu, N - 3, cu[1], np.sqrt(covu[1, 1]), cu[2] * DAY_S, np.sqrt(covu[2, 2]) * DAY_S))

    # per-sector secondary-primary offset, P fixed at M0 value
    say("\nPer-sector secondary - primary offset (P fixed at M0; offset = secondary O-C minus primary O-C):")
    Pm0 = c0[1]
    for s in (18, 58, 85):
        m = sec == s
        y = t[m] - Pm0 * (cyc[m] - E_REF)
        X = np.vstack([np.ones(m.sum()), is_sec[m].astype(float)]).T
        c, cov, chi = wls(X, y, e[m])
        say("  S%d: dsec = %.2f +- %.2f s (chi2=%.2f, dof=%d)" %
            (s, c[1] * DAY_S, np.sqrt(cov[1, 1]) * DAY_S, chi, m.sum() - 2))
    say("  predicted inner-orbit Romer delay P(K_B-K_A)/(pi c) = %.2f s" %
        (H.P_LEE * DAY_S * (136.17 - 93.16) / (np.pi * 299792.458)))

    # O-C vs Lee and vs M0
    say("\nSector-mean O-C (weighted), seconds:")
    oc_lee = (t - (H.T0_LEE + H.P_LEE * cyc)) * DAY_S
    model0 = design(cyc, is_sec, t, "M0")[0] @ c0
    oc_m0 = (t - model0) * DAY_S
    for s in (18, 58, 85):
        for lab, mt in (("primary", ~is_sec), ("secondary", is_sec)):
            m = (sec == s) & mt
            w = 1 / e[m] ** 2
            mu = np.sum(w * oc_lee[m]) / w.sum(); se = 1 / np.sqrt(w.sum()) * DAY_S
            mu2 = np.sum(w * oc_m0[m]) / w.sum()
            say("  S%d %-9s N=%2d  O-C(Lee)=%+7.2f +- %5.2f   O-C(M0)=%+7.2f" % (s, lab, m.sum(), mu, se, mu2))
    with open(os.path.join(OUT, TAG + "oc.csv"), "w") as fh:
        fh.write("sector,type,cycle,t_btjd,err_s,oc_lee_s,oc_m0_s\n")
        for i in range(N):
            fh.write("%d,%s,%.1f,%.7f,%.2f,%.2f,%.2f\n" % (sec[i], typ[i], cyc[i], t[i], e[i] * DAY_S,
                                                           oc_lee[i], oc_m0[i]))
    np.savez(os.path.join(OUT, TAG + "m0_fit.npz"), c=c0, cov=cov0, chi2=chi0)


if __name__ == "__main__":
    main()
