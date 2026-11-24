"""Step 4: injection-recovery LTTE limit (RULES 6) on the real eclipse times and scaled errors.

For each P3 (log grid 10-3000 d) and 1000 random phases, a circular LTTE
A sin(2 pi t / P3 + phi) is added to the real times and M0 (T0, P, dsec) is refitted.
Detection: chi2_M0 > chi2.ppf(0.9973, dof). Because the injected signal enters the M0
residuals linearly, chi2(A) is a quadratic in A for each phase, so the threshold amplitude
A*(phi) is solved exactly (equivalent to an infinitely fine A grid). A95 = smallest A
detected in >= 95% of phases = the 950th smallest A*(phi). A brute-force check on an A grid
is run for a few P3 values.

Minimum mass: (m3 sin i3)^3/(M12+m3)^2 = 4 pi^2 (A c)^3 / (G P3^2), sin i3 = 1, M12 = 3.068.
Writes independent/ltte_limit.csv and independent/out/ltte_summary.txt
"""
import os, sys
import numpy as np
from scipy import stats
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(__file__))
from s3_models import read_times, design, wls, TAG

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
DAY_S = 86400.0
NPH = 1000
G = 6.67430e-11; C = 299792458.0; MSUN = 1.98840987e30
M12 = 3.068


def m3_min(A_s, P3_d):
    fm = 4 * np.pi ** 2 * (A_s * C) ** 3 / (G * (P3_d * DAY_S) ** 2) / MSUN  # Msun
    g = lambda m: m ** 3 / (M12 + m) ** 2 - fm
    return brentq(g, 1e-8, 1e4)


def main():
    sec, typ, cyc, t, e, e0 = read_times()
    is_sec = typ == "secondary"
    N = len(t)
    X, _ = design(cyc, is_sec, t, "M0")
    Xw = X / e[:, None]
    # projector onto residual space in whitened coordinates
    Pm = np.eye(N) - Xw @ np.linalg.pinv(Xw)
    r0 = Pm @ (t / e)
    chi0 = r0 @ r0
    dof = N - X.shape[1]
    thr = stats.chi2.ppf(0.9973, dof)
    rng = np.random.default_rng(20261124)
    grid = np.logspace(np.log10(10.0), np.log10(3000.0), 600)
    rows = []
    for P3 in grid:
        arg = 2 * np.pi * t / P3
        us = Pm @ (np.sin(arg) / e)   # per unit A (days)
        uc = Pm @ (np.cos(arg) / e)
        phi = rng.uniform(0, 2 * np.pi, NPH)
        U = np.outer(np.cos(phi), us) + np.outer(np.sin(phi), uc)   # sin(arg+phi)
        a = np.sum(U * U, axis=1)
        b = U @ r0
        disc = b ** 2 + a * (thr - chi0)
        Astar = (-b + np.sqrt(disc)) / a       # days; chi0 < thr so a positive root exists
        Astar_s = np.sort(Astar) * DAY_S
        A95 = Astar_s[int(np.ceil(0.95 * NPH)) - 1]
        rows.append((P3, A95, m3_min(A95, P3)))
    rows = np.array(rows)
    with open(os.path.join(HERE if not TAG else OUT, TAG + "ltte_limit.csv"), "w") as fh:
        fh.write("P3_d,A95_s,m3min_msun\n")
        for P3, A, m in rows:
            fh.write("%.4f,%.4f,%.5f\n" % (P3, A, m))

    # brute-force check
    chk = []
    for P3 in (15.0, 120.0, 700.0, 2500.0):
        arg = 2 * np.pi * t / P3
        phi = np.random.default_rng(1).uniform(0, 2 * np.pi, NPH)
        Agrid = np.logspace(-0.5, 4, 400) / DAY_S
        frac = []
        for A in Agrid:
            det = 0
            for p in phi:
                y = t + A * np.sin(arg + p)
                _, _, chi = wls(X, y, e)
                det += chi > thr
            frac.append(det / NPH)
            if det / NPH >= 0.95:
                break
        frac = np.array(frac)
        Abf = Agrid[len(frac) - 1] * DAY_S
        # analytic for same phases
        us = Pm @ (np.sin(arg) / e); uc = Pm @ (np.cos(arg) / e)
        U = np.outer(np.cos(phi), us) + np.outer(np.sin(phi), uc)
        a = np.sum(U * U, 1); b = U @ r0
        As = np.sort((-b + np.sqrt(b ** 2 + a * (thr - chi0))) / a) * DAY_S
        chk.append((P3, Abf, As[949]))

    P3 = rows[:, 0]; m = rows[:, 2]; A = rows[:, 1]
    excl = (P3 >= 330) & (P3 <= 400)
    sel = (P3 >= 12) & ~excl
    out = open(os.path.join(OUT, TAG + "ltte_summary.txt"), "w")

    def say(*a):
        s = " ".join(str(x) for x in a); print(s); out.write(s + "\n")
    say("N=%d dof=%d chi2_M0(real)=%.3f threshold(99.73%%)=%.3f" % (N, dof, chi0, thr))
    say("P3 grid: %d log-spaced points 10-3000 d; %d random phases per P3" % (len(grid), NPH))
    say("median m3min, P3>=12 d excl 330-400 d: %.4f Msun (median A95 %.2f s) [n=%d]" %
        (np.median(m[sel]), np.median(A[sel]), sel.sum()))
    for lo, hi in ((12, 100), (100, 1000), (1000, 3000)):
        s2 = sel & (P3 >= lo) & (P3 <= hi)
        say("  median %4d-%4d d: m3min=%.4f Msun, A95=%.2f s [n=%d]" % (lo, hi, np.median(m[s2]), np.median(A[s2]), s2.sum()))
    j = np.argmax(np.where(sel, m, -1))
    say("max m3min (P3>=12, excl 330-400): %.4f Msun at P3=%.2f d (A95=%.2f s)" % (m[j], P3[j], A[j]))
    j2 = np.argmax(np.where(~excl, m, -1))
    say("max m3min (all P3 10-3000, excl 330-400): %.4f Msun at P3=%.2f d (A95=%.2f s)" % (m[j2], P3[j2], A[j2]))
    jj = np.argmax(m)
    say("max m3min incl 330-400: %.4f Msun at P3=%.2f d" % (m[jj], P3[jj]))
    say("min m3min: %.4f Msun at P3=%.2f d" % (m.min(), P3[np.argmin(m)]))
    say("brute-force check (P3, A95 brute on 400-pt log A grid, A95 analytic, same phases):")
    for c in chk:
        say("  P3=%7.1f  brute=%.2f s  analytic=%.2f s" % c)


if __name__ == "__main__":
    main()
