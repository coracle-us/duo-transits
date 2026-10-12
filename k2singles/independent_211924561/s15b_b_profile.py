"""Profile chi2 vs impact parameter b (free k, a/R*, mid-times), with a/R* starts that give T14~3 h,
plus the circular-orbit density posterior from the free-shape 4-event chain scaled to each surviving P,
and the minimum eccentricity (omega=90 deg) needed to reconcile it with the Gaia FLAME density."""
import os, json, sys
import numpy as np
from scipy.optimize import minimize
import lcio, tfit
import s07_fit_events as F

HERE = lcio.HERE
ev = ["K2C5_K2SFF", "K2C16_EVEREST", "S46_SPOC2min_PDC", "S72_SPOC2min_PDC"]
Ds = [F.prep(e) for e in ev]
fr = json.load(open(os.path.join(HERE, "fit_results.json")))
tc0 = np.array([x[0] for x in fr["joint4_primary"]["tc"]])
chi2 = lambda k, b, aR, tcs: -2 * sum(F.lnlike_event(D, tc, k, b, aR)[0] for D, tc in zip(Ds, tcs))

prof = []
for b in [0.0, 0.3, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15]:
    best = None
    for k0 in [0.045, 0.06, 0.1, 0.15, 0.25]:
        if b >= 1 + k0 - 0.01:
            continue
        for T in [2.6, 3.0]:
            aR0 = tfit.PREF / (np.pi * T / 24) * np.sqrt((1 + k0) ** 2 - b ** 2)
            x0 = np.r_[k0, np.log(aR0), tc0]
            fn = lambda x: chi2(x[0], b, np.exp(x[1]), x[2:]) if (0.005 < x[0] < 0.6 and b < 1 + x[0]) else 1e12
            r = minimize(fn, x0, method="Nelder-Mead", options=dict(maxiter=3000, xatol=1e-7, fatol=1e-3))
            if best is None or r.fun < best.fun:
                best = r
    k, aR = best.x[0], np.exp(best.x[1])
    prof.append(dict(b=b, chi2=float(best.fun), k=float(k), T14_h=float(tfit.t14(tfit.PREF, aR, k, b) * 24), T23_h=float(tfit.t23(tfit.PREF, aR, k, b) * 24),
                     rho_circ_at_P100=float(tfit.rho_star(tfit.PREF, aR))))
cmin = min(p["chi2"] for p in prof)
for p in prof:
    p["dchi2"] = p["chi2"] - cmin
    print("b=%.2f dchi2=%7.2f k=%.4f T14=%.2f h T23=%.2f h rho_circ(P=100 d)=%.3f" % (p["b"], p["dchi2"], p["k"], p["T14_h"], p["T23_h"], p["rho_circ_at_P100"]))

# density posterior from the free-shape chain
ch = np.load(os.path.join(HERE, "cache", "chain_joint4_primary.npy"))
rho100 = tfit.rho_star(tfit.PREF, np.exp(ch[:, 2]))
RHO = json.load(open(os.path.join(HERE, "density_check.json")))["rho_star_gcc"]
res = dict(b_profile=prof, rho_star_gaia=RHO, per_P={})
for P in (140.49401, 70.24701):
    r = rho100 * P / tfit.PREF
    g = (r / RHO) ** (1 / 3)
    emin = (g ** 2 - 1) / (g ** 2 + 1)
    pc = lambda x: [float(np.percentile(x, q)) for q in (2.5, 16, 50, 84, 97.5)]
    res["per_P"][str(P)] = dict(rho_circ_pct=pc(r), ratio_pct=pc(r / RHO), e_min_pct=pc(np.clip(emin, 0, None)),
                                frac_rho_below_gaia=float(np.mean(r < RHO)))
    print(P, json.dumps(res["per_P"][str(P)]))
json.dump(res, open(os.path.join(HERE, "b_profile.json"), "w"), indent=1)
