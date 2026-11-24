"""Check the per-eclipse optimizer: noise-free recovery and lm fit vs a dense shift scan."""
import os, sys, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import hiplib as H
from s2_times import fit_one
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
tmpl = H.Template(pickle.load(open(os.path.join(OUT, "template_s58.pkl"), "rb")))
rng = np.random.default_rng(5)
for texp in (120.0, 200.0):
    cad = texp/86400; model = H.eclipse_model_factory(tmpl, "P", texp)
    dt = np.arange(-H.FIT_HALF*H.P_LEE, H.FIT_HALF*H.P_LEE, cad)
    truth = 40/86400
    y0 = model(dt, truth, 1.0, 1.0, 0.0)
    p, _ = fit_one(dt, y0, np.ones_like(dt), model, [0, 1, 1, 0])
    print("texp %g noise-free: recovered %.3f s" % (texp, p[0]*86400))
    d_lm, d_scan = [], []
    for k in range(40):
        y = y0 + rng.normal(0, 0.9e-3, len(dt))
        p, _ = fit_one(dt, y, np.ones_like(dt), model, [0, 1, 1, 0])
        grid = np.linspace(-20, 100, 1201)/86400
        chi = []
        for s in grid:
            vals = model(dt, s, 1, 0, 0) - 0  # = scale*(T-1) with c0=0
            A = np.vstack([np.ones_like(dt), dt, vals + 0*dt]).T
            # model = c0 + c1 dt + scale*(T-1); vals currently = 1*(T-1)
            c, *_ = np.linalg.lstsq(A, y, rcond=None)
            chi.append(np.sum((y - A@c)**2))
        d_lm.append(p[0]*86400); d_scan.append(grid[np.argmin(chi)]*86400)
    d_lm, d_scan = np.array(d_lm), np.array(d_scan)
    print("  lm: mean %.2f sd %.2f | scan: mean %.2f sd %.2f | max |lm-scan| %.2f s" %
          (d_lm.mean(), d_lm.std(), d_scan.mean(), d_scan.std(), np.abs(d_lm-d_scan).max()))
