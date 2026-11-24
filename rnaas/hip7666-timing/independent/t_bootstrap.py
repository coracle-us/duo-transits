"""Compare the block-bootstrap error with the true scatter on white-noise synthetic eclipses."""
import os, sys, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import hiplib as H
from s2_times import fit_one
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
tmpl = H.Template(pickle.load(open(os.path.join(OUT, "template_s58.pkl"), "rb")))
rng = np.random.default_rng(9)
for texp, kind in ((200.0, "P"), (120.0, "P"), (200.0, "S")):
    cad = texp/86400; model = H.eclipse_model_factory(tmpl, kind, texp)
    dt = np.arange(-H.FIT_HALF*H.P_LEE, H.FIT_HALF*H.P_LEE, cad); n = len(dt)
    y0 = model(dt, 0.0, 1.0, 1.0, 0.0)
    fits, berr = [], {1: [], 22 if texp == 200 else 36: []}
    for k in range(30):
        y = y0 + rng.normal(0, 0.9e-3, n)
        p, _ = fit_one(dt, y, np.ones(n), model, [0, 1, 1, 0])
        fits.append(p[0]*86400)
        best = model(dt, *p); res = y - best
        for L in berr:
            nb = int(np.ceil(n/L)); sh = []
            for b in range(100):
                st = rng.integers(0, n-L+1, nb)
                idx = (st[:, None]+np.arange(L)).ravel()[:n]
                pb, _ = fit_one(dt, best+res[idx], np.ones(n), model, p)
                sh.append(pb[0]*86400)
            berr[L].append(np.std(sh, ddof=1))
    print("texp %g %s n=%d: true scatter %.2f s; " % (texp, kind, n, np.std(fits, ddof=1)) +
          "; ".join("boot L=%d mean err %.2f s" % (L, np.mean(v)) for L, v in berr.items()))
