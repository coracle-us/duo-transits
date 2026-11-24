"""Impact-parameter profile of the joint transit-shape fit: chi2 with b held fixed on a grid, a/R* and
Rp/R* free (same segments and model as s05/s15). Gives the upper limit on the transit impact
parameter, which decides whether an eccentric orbit still produces a secondary eclipse.
Usage: s18_b_profile.py TIC"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import least_squares
from common import TARGETS, OUT
from s05_fit_times import gather, model

tic = int(sys.argv[1])
dur = TARGETS[tic]["dur"]
eph = json.load(open(os.path.join(OUT, f"ephem_{tic}.json")))
P, T0 = eph["P"], eph["T0_E0"]
segs = gather(tic, P, T0)
use = [s for s in segs if s["timing"] and s["cov"] > 0.7 and s["ing"] > 0.5 and s["egr"] > 0.5]
k0 = np.sqrt(TARGETS[tic]["depth"])
n = sum(len(s["t"]) for s in use)


def fit(b):
    def resid(x):
        rp, a = x[0], x[1]
        return np.concatenate([(s["f"] - model((rp, a, b), x[2 + i], s["t"], s["cad"], P) * (1 + x[2 + len(use) + i]))
                               / s["e"] for i, s in enumerate(use)])
    a0 = P / (np.pi * dur) * np.sqrt(max(1 - b**2, 0.05))
    x0 = [k0, a0] + [s["tpred"] for s in use] + [0.0] * len(use)
    lo = [0.01, 1.5] + [s["tpred"] - 0.3 * dur for s in use] + [-0.01] * len(use)
    hi = [0.5, 500] + [s["tpred"] + 0.3 * dur for s in use] + [0.01] * len(use)
    r = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac")
    return 2 * r.cost, r.x[0], r.x[1]


grid = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
res = [(b,) + fit(b) for b in grid]
cmin = min(r[1] for r in res)
scale = cmin / (n - (3 + 2 * len(use)))
L = [f"TIC {tic}: b profile, {len(use)} full transits, {n} points, chi2/dof at best = {scale:.2f}"]
for b, c, rp, a in res:
    L.append(f"  b = {b:.2f}: Rp/R* = {rp:.4f}  a/R* = {a:6.2f}  chi2 = {c:.1f}  "
             f"scaled dchi2 = {(c - cmin) / scale:6.1f}")
# 2-sigma (dchi2 = 4) upper limit by linear interpolation
d = np.array([(c - cmin) / scale for _, c, _, _ in res])
bs = np.array(grid)
i = np.argmax(d > 4) if (d > 4).any() else None
if i:
    bup = np.interp(4, d[i - 1:i + 1], bs[i - 1:i + 1])
    L.append(f"  b < {bup:.2f} at 2 sigma (scaled dchi2 = 4)")
else:
    L.append("  no 2-sigma upper limit on b within the grid")
txt = "\n".join(L)
open(os.path.join(OUT, f"bprofile_{tic}.txt"), "w").write(txt + "\n")
print(txt)
