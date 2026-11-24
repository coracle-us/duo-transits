"""Synthetic checks of the eclipse timer (sign, integration, no bias).
1. Pure synthetic: S58 template (symmetrised eclipse shapes) placed at Lee + dP (primary) and
   Lee + 0.5P + dS (secondary), integrated over 1800 s / 200 s / 120 s, white noise added.
2. Real S85 prepared data with all time stamps shifted by +30 s.
"""
import os, sys, pickle, io
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import hiplib as H
from s2_times import time_sector

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
tmpl = H.Template(pickle.load(open(os.path.join(OUT, "template_s58.pkl"), "rb")))
log = io.StringIO()
rng = np.random.default_rng(11)
dP, dS = 40.0 / 86400, -20.0 / 86400

for texp, noise in ((200.0, 0.9e-3), (1800.0, 1.5e-3)):
    cad = texp / 86400
    t = np.arange(2882.3, 2960.0, cad)
    nsub = 31
    offs = (np.arange(nsub) - (nsub - 1) / 2) / nsub * cad
    y = np.zeros_like(t)
    for o in offs:
        tt = t + o
        ph = H.phase(tt)
        f = np.full_like(t, np.nan)
        # piecewise: use eclipse shapes inside |dphi|<0.07, full template elsewhere
        f = tmpl.full(ph)
        dpp = (tt - (H.T0_LEE + H.P_LEE * np.round((tt - H.T0_LEE) / H.P_LEE)) - dP) / H.P_LEE
        mP = np.abs(dpp) < 0.07
        f[mP] = tmpl.eclipse("P", dpp[mP])
        Es = np.round((tt - H.T0_LEE) / H.P_LEE - 0.5)
        dss = (tt - (H.T0_LEE + H.P_LEE * (Es + 0.5)) - dS) / H.P_LEE
        mS = np.abs(dss) < 0.07
        f[mS] = tmpl.eclipse("S", dss[mS])
        y += f / nsub
    y += rng.normal(0, noise, len(t))
    e = np.full_like(t, noise)
    rows, dropped = time_sector(58, tmpl, log, data=(t, y, e), texp=texp, nboot=100)
    for k, d in (("P", dP), ("S", dS)):
        sh = np.array([r["shift_s"] for r in rows if r["type"] == k])
        er = np.array([r["err_raw"] for r in rows if r["type"] == k]) * 86400
        print("texp=%5.0f %s: injected %+.1f s, recovered mean %+.2f s, scatter %.2f s, mean boot err %.2f s, N=%d"
              % (texp, k, d * 86400, sh.mean(), sh.std(ddof=1), er.mean(), len(sh)))

