"""Does a circular orbit fit the transit shape when the stellar density is held at the Gaia value?
Joint shape fit (same segments and model as s05) twice: (1) free a/R*, b; (2) a/R* fixed from the
Gaia-based density (out/radius_<tic>.txt median), b free. Reports chi2 of each and the difference.
Usage: s15_rho_fixed.py TIC"""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import least_squares
from common import TARGETS, OUT
from s05_fit_times import gather, model

tic = int(sys.argv[1])
dur = TARGETS[tic]["dur"]
eph = json.load(open(os.path.join(OUT, f"ephem_{tic}.json")))
P, T0 = eph["P"], eph["T0_E0"]
rho = float(re.search(r"rho = ([\d.]+)", open(os.path.join(OUT, f"radius_{tic}.txt")).read()).group(1))
a_fix = (rho * 1.41 * 6.674e-8 * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
segs = gather(tic, P, T0)
use = [s for s in segs if s["timing"] and s["cov"] > 0.7 and s["ing"] > 0.5 and s["egr"] > 0.5]
k0 = np.sqrt(TARGETS[tic]["depth"])


def fit(fix_a):
    def resid(x):
        rp, b = x[0], x[1]
        a = a_fix if fix_a else x[2]
        o = 2 if fix_a else 3
        return np.concatenate([(s["f"] - model((rp, a, b), x[o + i], s["t"], s["cad"], P) * (1 + x[o + len(use) + i]))
                               / s["e"] for i, s in enumerate(use)])
    best = None
    for b0 in [0.05, 0.3, 0.6, 0.85]:
        x0 = [k0, b0] + ([] if fix_a else [P / (np.pi * dur)]) + [s["tpred"] for s in use] + [0.0] * len(use)
        lo = [0.01, 0.0] + ([] if fix_a else [1.5]) + [s["tpred"] - 0.3 * dur for s in use] + [-0.01] * len(use)
        hi = [0.5, 1.2] + ([] if fix_a else [500]) + [s["tpred"] + 0.3 * dur for s in use] + [0.01] * len(use)
        r = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac")
        if best is None or r.cost < best.cost:
            best = r
    return best


free, fixd = fit(False), fit(True)
n = sum(len(s["t"]) for s in use)
c_free, c_fix = 2 * free.cost, 2 * fixd.cost
L = [f"TIC {tic}: {len(use)} full transits, {n} points; Gaia rho = {rho:.2f} rho_sun -> a/R* = {a_fix:.2f} (circular)",
     f"  free:  Rp/R* = {free.x[0]:.4f}  b = {free.x[1]:.2f}  a/R* = {free.x[2]:.2f}  chi2 = {c_free:.1f}",
     f"  fixed: Rp/R* = {fixd.x[0]:.4f}  b = {fixd.x[1]:.2f}  a/R* = {a_fix:.2f} (fixed)  chi2 = {c_fix:.1f}",
     f"  delta chi2 (fixed - free) = {c_fix - c_free:.1f} for 1 fewer parameter; "
     f"chi2/dof(free) = {c_free / (n - len(free.x)):.2f}  "
     f"scaled delta chi2 = {(c_fix - c_free) / (c_free / (n - len(free.x))):.1f}"]
txt = "\n".join(L)
open(os.path.join(OUT, f"rhofixed_{tic}.txt"), "w").write(txt + "\n")
print(txt)
