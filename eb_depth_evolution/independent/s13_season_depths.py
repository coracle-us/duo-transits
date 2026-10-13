"""Step 13: per-season model depths and geometric thresholds from the TESS posterior.

Depth definition: D = 1 - F_min/F_out, the minimum of the best-fit physical model at mid-eclipse
(instantaneous, i.e. without cadence smearing), for that season's inclination.
 * 'obs'  : as seen in our r = 2.0 px TESScut aperture (l3 = nominal 0.035)
 * 'undil': corrected to zero third light
Also: primary/secondary minimum sky separations (impact parameters in units of a),
i_onset(P) (primary disappears: d_min = r1 + r2), i_onset(S), i_total(P) (d_min = |r2 - r1|).
usage: s13_season_depths.py [tag]
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy.optimize import brentq
import ebmodel as M

H = os.path.dirname(os.path.abspath(__file__))
tag = sys.argv[1] if len(sys.argv) > 1 else "nominal"
pj = json.load(open(os.path.join(H, "fit_tess_%s.json" % tag)))
ch = np.load(os.path.join(H, "data", "chain_%s.npz" % tag))
names = list(ch["names"]); flat = ch["flat"]
fixed = json.loads(str(ch["fixed"]))
kw = pj["_fitter_kw"]; u1, u2, l3 = kw["u1"], kw["u2"], kw["l3"]
P = fixed["P"]


def dmin(inc, e, w, which):
    fc = (np.pi / 2 if which == "P" else 3 * np.pi / 2) - w
    f = fc + np.linspace(-0.3, 0.3, 4001)
    r = (1 - e ** 2) / (1 + e * np.cos(f))
    s = np.sin(w + f)
    return np.min(r * np.sqrt(1 - np.sin(np.radians(inc)) ** 2 * s ** 2))


def depth(inc, e, w, r1, r2, J, which, l3v):
    T0 = 0.0
    tc = 0.0 if which == "P" else M.secondary_offset(P, e, w)
    t = tc + np.linspace(-0.6, 0.6, 6001)
    m = M.lc(t, P, T0, e, w, inc, r1, r2, J, u1[0], u1[1], u2[0], u2[1], l3=l3v)
    return 1 - m.min()


def thresholds(e, w, r1, r2):
    out = {}
    g = lambda inc, which, target: dmin(inc, e, w, which) - target
    for which in ("P", "S"):
        out["ion_" + which] = brentq(g, 60, 89.999999, args=(which, r1 + r2)) if g(60, which, r1 + r2) > 0 else np.nan
    dd = abs(r2 - r1)
    out["itot_P"] = brentq(g, 60, 90, args=("P", dd)) if g(90, "P", dd) < 0 else 90.0
    out["itot_S"] = brentq(g, 60, 90, args=("S", dd)) if g(90, "S", dd) < 0 else 90.0
    return out


rng = np.random.default_rng(3)
idx = rng.choice(len(flat), 1500, replace=False)
rows = []
for j in idx:
    p = dict(fixed); p.update(zip(names, flat[j]))
    row = dict(r1=p["r1"], r2=p["r2"], J=p["J"], esw=p["esw"])
    for s in (2021, 2023, 2026):
        ecw = p["ecw%d" % s]
        e = np.hypot(ecw, p["esw"]); w = np.arctan2(p["esw"], ecw)
        inc = p["i%d" % s]
        row["i%d" % s] = inc
        for which in ("P", "S"):
            row["D%s%d_obs" % (which, s)] = depth(inc, e, w, p["r1"], p["r2"], p["J"], which, l3)
            row["D%s%d_undil" % (which, s)] = depth(inc, e, w, p["r1"], p["r2"], p["J"], which, 0.0)
            row["b%s%d" % (which, s)] = dmin(inc, e, w, which)
        if s == 2023:
            th = thresholds(e, w, p["r1"], p["r2"])
            row.update(th)
            row["e"] = e; row["w_deg"] = np.degrees(w)
    rows.append(row)
R = pd.DataFrame(rows)
R.to_csv(os.path.join(H, "data", "season_depth_samples_%s.csv" % tag), index=False)
summ = {}
for c in R.columns:
    q = np.nanpercentile(R[c], [2.5, 16, 50, 84, 97.5])
    summ[c] = q.tolist()
    print("%-14s %10.5f  +%.5f -%.5f   [95%%: %.5f .. %.5f]" % (c, q[2], q[3] - q[2], q[2] - q[1], q[0], q[4]))
json.dump(summ, open(os.path.join(H, "season_depths_%s.json" % tag), "w"), indent=1)
print("fraction r2 > r1: %.3f" % np.mean(R.r2 > R.r1))
