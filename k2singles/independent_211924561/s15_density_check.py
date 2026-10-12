"""Shape + stellar-density consistency for each surviving period.

(a) Free-shape profile: chi2 of the 4-event joint data vs impact parameter b on a fixed grid
    (k, a/R*, mid-times free) -> tells whether a flat-bottomed (low-b) transit is allowed.
(b) For each surviving P: MCMC with a/R* tied to P and the Gaia DR3 FLAME star
    (R=1.739, M=1.424 -> rho=0.382 g/cc, 15% Gaussian prior on rho), circular orbit; free k, b, mid-times.
    Report b, k posteriors and the chi2 penalty relative to the free-shape fit.
(c) Implied circular-orbit stellar density from the free-shape posterior, scaled to each P.
"""
import os, json, sys
import numpy as np
import emcee
from scipy.optimize import minimize
import lcio, tfit
sys.argv = [sys.argv[0], "none"]
import s07_fit_events as F

HERE = lcio.HERE
G = 6.674e-8; MSUN = 1.989e33; RSUN = 6.957e10
R, M = 1.739, 1.424
RHO = M * MSUN / (4 / 3 * np.pi * (R * RSUN) ** 3)
SRHO = 0.15 * RHO
ev = ["K2C5_K2SFF", "K2C16_EVEREST", "S46_SPOC2min_PDC", "S72_SPOC2min_PDC"]
Ds = [F.prep(e) for e in ev]
alias = json.load(open(os.path.join(HERE, "alias_table.json")))
survP = [p for p in alias["periods"] if p["status"] == "survives"]
fr = json.load(open(os.path.join(HERE, "fit_results.json")))
tc0 = np.array([x[0] for x in fr["joint4_primary"]["tc"]])


def chi2_total(k, b, aR, tcs):
    return -2 * sum(F.lnlike_event(D, tc, k, b, aR)[0] for D, tc in zip(Ds, tcs))


def aR_of(P, rho):
    return (G * rho * (P * 86400.) ** 2 / (3 * np.pi)) ** (1 / 3)


out = dict(rho_star_gcc=RHO, srho=SRHO)
# (a) profile in b
best_free = None
prof = []
for b in [0.0, 0.3, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2]:
    bestb = None
    for k0 in [0.05, 0.15]:
        if b > 1 + k0 - 0.01:
            continue
        for la in [np.log(60), np.log(110)]:
            x0 = np.r_[k0, la, tc0]
            fn = lambda x: chi2_total(x[0], b, np.exp(x[1]), x[2:]) if (0.005 < x[0] < 0.6 and b < 1 + x[0] and 0.4 < x[1] < 7) else 1e12
            r = minimize(fn, x0, method="Nelder-Mead", options=dict(maxiter=2500, xatol=1e-6, fatol=1e-3))
            if bestb is None or r.fun < bestb.fun:
                bestb = r
    if bestb is None:
        continue
    k, la = bestb.x[0], bestb.x[1]
    prof.append(dict(b=b, chi2=float(bestb.fun), k=float(k), aR_Pref=float(np.exp(la)), T14_h=float(tfit.t14(tfit.PREF, np.exp(la), k, b) * 24),
                     T23_h=float(tfit.t23(tfit.PREF, np.exp(la), k, b) * 24), rho_at_Pref=float(tfit.rho_star(tfit.PREF, np.exp(la)))))
    print("b=%.2f chi2=%.2f k=%.4f T14=%.2fh T23=%.2fh rho(P=100)=%.4f" % (b, bestb.fun, k, prof[-1]["T14_h"], prof[-1]["T23_h"], prof[-1]["rho_at_Pref"]))
cmin = min(p["chi2"] for p in prof)
for p in prof:
    p["dchi2"] = p["chi2"] - cmin
    # implied density (circular) at each surviving P: rho scales linearly with P at fixed shape
    p["rho_implied_at_P"] = {f"{s['P']:.3f}": p["rho_at_Pref"] * s["P"] / tfit.PREF for s in survP}
out["b_profile"] = prof
ntot = sum(len(D["t"]) for D in Ds)
print("total points", ntot)

# (b) density-tied MCMC per surviving P
rng = np.random.default_rng(3)
out["per_period"] = []
for s in survP:
    P = s["P"]

    def lnp(th):
        k, b, lrho = th[:3]
        tcs = th[3:]
        rho = np.exp(lrho)
        if not (0.005 < k < 0.8 and 0 <= b < 1 + k):
            return -np.inf
        aR_true = aR_of(P, rho)
        aR_ref = aR_true * tfit.PREF / P   # same transit chord velocity at the reference period
        if aR_ref < 1.5:
            return -np.inf
        lp = -0.5 * ((rho - RHO) / SRHO) ** 2
        return lp - 0.5 * chi2_total(k, b, aR_ref, tcs)

    best = None
    for b0 in [0.98, 1.08]:
        for k0 in [0.08, 0.15]:
            x0 = np.r_[k0, b0, np.log(RHO), tc0]
            if b0 >= 1 + k0:
                continue
            r = minimize(lambda x: -lnp(x), x0, method="Nelder-Mead", options=dict(maxiter=6000, xatol=1e-7, fatol=1e-4))
            if best is None or r.fun < best.fun:
                best = r
    nw, nd = 40, len(best.x)
    pos = best.x + 1e-4 * rng.standard_normal((nw, nd))
    sam = emcee.EnsembleSampler(nw, nd, lnp)
    sam.run_mcmc(pos, 3000, progress=False)
    c = sam.get_chain(discard=1200, thin=10, flat=True)
    q = lambda x: [float(np.percentile(x, 50)), float(np.percentile(x, 50) - np.percentile(x, 16)), float(np.percentile(x, 84) - np.percentile(x, 50))]
    chi2_best = 2 * best.fun - (0)  # includes the rho prior term (~0 at best)
    rec = dict(P=P, aR=aR_of(P, RHO), T14_central_h=float(tfit.t14(P, aR_of(P, RHO), 0.045, 0.0) * 24),
               b=q(c[:, 1]), k=q(c[:, 0]), rho=q(np.exp(c[:, 2])), b_2p5=float(np.percentile(c[:, 1], 2.5)),
               k_97p5=float(np.percentile(c[:, 0], 97.5)), chi2_best=float(chi2_best), dchi2_vs_free=float(chi2_best - cmin),
               Rp_RJ=q(c[:, 0] * R * 9.731), grazing_frac=float(np.mean(c[:, 1] > 1 - c[:, 0])))
    out["per_period"].append(rec)
    print(json.dumps(rec))
json.dump(out, open(os.path.join(HERE, "density_check.json"), "w"), indent=1)
