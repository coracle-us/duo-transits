"""MCMC fits: (a) each event alone, in two independent reductions each; (b) joint fit of the three events.

Shape is parametrised with a reference period PREF=100 d (k, b, a/R* at PREF); T14, T23 and depth are
period-independent outputs. Polynomial (quadratic) baselines are profiled out per event.
LD fixed: Kepler (0.31, 0.31), TESS (0.23, 0.31) for an F star (Teff~6300 K, logg~4.1).
"""
import os, sys, json
import numpy as np
from astropy.io import fits
import emcee
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lcio, tfit

rng = np.random.default_rng(42)
HW = 0.7
U_K = (0.31, 0.31); U_T = (0.23, 0.31)


def window(d, tc, hw=HW):
    w = np.abs(d["t"] - tc) < hw
    return d["t"][w], d["f"][w]


def everest_noq(c):
    d = fits.getdata(lcio._find(f"hlsp_everest_k2_llc_211924561-c{c:02d}_kepler_v2.0_lc.fits"), 1)
    q = (d["QUALITY"] & (1 << 20)) == 0
    return lcio._norm(d["TIME"] - 2167.0, np.asarray(d["FLUX"], float), None, q)


def k2sap_hard(c):
    d = fits.getdata(lcio._find(f"ktwo211924561-c{c:02d}_llc.fits"), 1)
    hard = 1 | 2 | 4 | 8 | 32 | 128 | (1 << 20)
    q = (d["SAP_QUALITY"] & hard) == 0
    return lcio._norm(d["TIME"] - 2167.0, np.asarray(d["SAP_FLUX"], float), None, q)


def sap_tess(s):
    return lcio.tess_spoc2min(s, col="SAP_FLUX")


SETS = {
    "K2C5_K2SFF": dict(load=lambda: lcio.k2_k2sff(5), tc=179.90, u=U_K, exp=0.020434),
    "K2C5_EVEREST": dict(load=lambda: everest_noq(5), tc=179.90, u=U_K, exp=0.020434),
    "K2C16_EVEREST": dict(load=lambda: everest_noq(16), tc=1163.36, u=U_K, exp=0.020434),
    "K2C16_SAP": dict(load=lambda: k2sap_hard(16), tc=1163.36, u=U_K, exp=0.020434),
    "S46_SPOC2min_PDC": dict(load=lambda: lcio.tess_spoc2min(46), tc=2568.30, u=U_T, exp=None),
    "S46_QLP": dict(load=lambda: lcio.tess_qlp(46), tc=2568.30, u=U_T, exp=600 / 86400.),
    "S72_SPOC2min_PDC": dict(load=lambda: lcio.tess_spoc2min(72), tc=3270.75, u=U_T, exp=None),
    "S72_QLP": dict(load=lambda: lcio.tess_qlp(72), tc=3270.75, u=U_T, exp=200 / 86400.),
}


def prep(name):
    s = SETS[name]
    d = s["load"]()
    t, f = window(d, s["tc"])
    o = np.argsort(t); t, f = t[o], f[o]
    # white noise from out-of-transit points (|dt|>0.15) after quadratic baseline
    oot = np.abs(t - s["tc"]) > 0.15
    c = np.polyfit(t[oot] - s["tc"], f[oot], 2)
    r = f[oot] - np.polyval(c, t[oot] - s["tc"])
    sig = 1.4826 * np.median(np.abs(r - np.median(r)))
    # clip >5 sigma outliers out of transit only
    keep = ~(oot & (np.abs(f - np.polyval(c, t - s["tc"])) > 5 * sig))
    return dict(t=t[keep], f=f[keep], sig=sig, u=s["u"], exp=s["exp"], tc0=s["tc"])


def lnlike_event(D, tc, k, b, aR):
    m = tfit.model(D["t"], tc, k, b, aR, D["u"], D["exp"])
    base = tfit.profile_baseline(D["t"], D["f"], m, D["tc0"])
    r = D["f"] - m * base
    return -0.5 * np.sum((r / D["sig"]) ** 2), m * base


def lnprior_shape(k, b, laR):
    if not (0.005 < k < 0.25 and 0 <= b < 1 + k and np.log(1.5) < laR < np.log(1000)):
        return -np.inf
    return 0.0


def run(events, joint_name, nwalk=40, nstep=6000, burn=2000):
    Ds = [prep(e) for e in events]
    ne = len(Ds)

    def lnprob(th):
        k, b, laR = th[:3]
        tcs = th[3:]
        lp = lnprior_shape(k, b, laR)
        if not np.isfinite(lp):
            return -np.inf
        aR = np.exp(laR)
        if np.any(np.abs(tcs - np.array([D["tc0"] for D in Ds])) > 0.15):
            return -np.inf
        return lp + sum(lnlike_event(D, tc, k, b, aR)[0] for D, tc in zip(Ds, tcs))

    p0 = np.array([0.042, 0.4, np.log(60.)] + [D["tc0"] for D in Ds])
    # quick optimisation from a few starts
    from scipy.optimize import minimize
    best = None
    for b0 in [0.1, 0.5, 0.8, 0.95]:
        for la in [np.log(30), np.log(60), np.log(120)]:
            x0 = p0.copy(); x0[1] = b0; x0[2] = la
            res = minimize(lambda x: -lnprob(x), x0, method="Nelder-Mead", options=dict(maxiter=6000, xatol=1e-6, fatol=1e-4))
            if best is None or res.fun < best.fun:
                best = res
    ndim = len(p0)
    pos = best.x + 1e-4 * rng.standard_normal((nwalk, ndim)) * np.r_[0.05, 1, 1, np.ones(ne) * 0.02]
    pos[:, 1] = np.abs(pos[:, 1])
    sampler = emcee.EnsembleSampler(nwalk, ndim, lnprob)
    sampler.run_mcmc(pos, nstep, progress=False)
    ch = sampler.get_chain(discard=burn, thin=10, flat=True)
    k, b, aR = ch[:, 0], ch[:, 1], np.exp(ch[:, 2])
    T14 = tfit.t14(tfit.PREF, aR, k, b) * 24
    T23 = tfit.t23(tfit.PREF, aR, k, b) * 24
    # depth of the model at mid-transit (no smearing) for the first event's LD
    dep = np.array([1 - tfit.model(np.array([0.0]), 0.0, kk, bb, aa, Ds[-1]["u"])[0] for kk, bb, aa in zip(k[::20], b[::20], aR[::20])]) * 1e3
    rho_ref = tfit.rho_star(tfit.PREF, aR)  # g/cc at PREF; scales linearly with P
    q = lambda x: [float(np.percentile(x, 50)), float(np.percentile(x, 50) - np.percentile(x, 16)), float(np.percentile(x, 84) - np.percentile(x, 50))]
    out = dict(events=events, lnp_max=float(-best.fun), acc=float(np.mean(sampler.acceptance_fraction)),
               k=q(k), b=q(b), aR_at_Pref=q(aR), T14_h=q(T14), T23_h=q(T23), T23_over_T14=q(T23 / T14),
               depth_ppt=q(dep), rho_at_Pref_gcc=q(rho_ref), sig_ppm=[float(D["sig"] * 1e6) for D in Ds],
               tc=[q(ch[:, 3 + i]) for i in range(ne)],
               b_samples_pct=[float(np.percentile(b, p)) for p in (2.5, 16, 50, 84, 97.5)])
    np.save(os.path.join(lcio.HERE, "cache", f"chain_{joint_name}.npy"), ch)
    # plot
    med = np.median(ch, 0)
    fig, ax = plt.subplots(1, ne, figsize=(5 * ne, 3.6), squeeze=False)
    for i, D in enumerate(Ds):
        ll, mb = lnlike_event(D, med[3 + i], med[0], med[1], np.exp(med[2]))
        m = tfit.model(D["t"], med[3 + i], med[0], med[1], np.exp(med[2]), D["u"], D["exp"])
        base = mb / m
        x = (D["t"] - med[3 + i]) * 24
        a = ax[0, i]
        a.plot(x, 1e3 * (D["f"] / base - 1), ".", ms=2, color="0.6")
        # bin to 20 min
        e = np.arange(x.min(), x.max() + 1 / 3., 1 / 3.)
        idx = np.digitize(x, e)
        xb = [x[idx == j].mean() for j in np.unique(idx)]; yb = [np.mean(1e3 * (D["f"] / base - 1)[idx == j]) for j in np.unique(idx)]
        a.plot(xb, yb, "o", ms=4, color="C0")
        xx = np.linspace(x.min(), x.max(), 2000)
        a.plot(xx, 1e3 * (tfit.model(med[3 + i] + xx / 24, med[3 + i], med[0], med[1], np.exp(med[2]), D["u"], D["exp"]) - 1), "r-")
        a.set_xlim(-8, 8); a.set_ylim(-3.5, 1.5)
        a.set_title(f"{events[i]}  tc={med[3+i]:.4f}", fontsize=9)
        a.set_xlabel("hours from mid-transit"); a.set_ylabel("ppt")
    plt.tight_layout(); plt.savefig(os.path.join(lcio.HERE, f"fig_fit_{joint_name}.png"), dpi=100); plt.close()
    return out


if __name__ == "__main__":
    which = sys.argv[1]
    allres = {}
    path = os.path.join(lcio.HERE, "fit_results.json")
    if os.path.exists(path):
        allres = json.load(open(path))
    if which == "single":
        for e in SETS:
            allres["single_" + e] = run([e], "single_" + e, nstep=4000, burn=1500)
            print(e, json.dumps({k: allres["single_" + e][k] for k in ["tc", "T14_h", "depth_ppt", "b", "k"]}))
    elif which == "joint":
        allres["joint_primary"] = run(["K2C5_K2SFF", "S46_SPOC2min_PDC", "S72_SPOC2min_PDC"], "joint_primary", nstep=8000, burn=3000)
        allres["joint_alt"] = run(["K2C5_EVEREST", "S46_QLP", "S72_QLP"], "joint_alt", nstep=8000, burn=3000)
        allres["joint_tess_only"] = run(["S46_SPOC2min_PDC", "S72_SPOC2min_PDC"], "joint_tess_only", nstep=8000, burn=3000)
        for k in ["joint_primary", "joint_alt", "joint_tess_only"]:
            print(k, json.dumps(allres[k]))
    elif which == "c16single":
        for e in ["K2C16_EVEREST", "K2C16_SAP"]:
            allres["single_" + e] = run([e], "single_" + e, nstep=4000, burn=1500)
            print(e, json.dumps({k: allres["single_" + e][k] for k in ["tc", "T14_h", "depth_ppt", "b", "k"]}))
    elif which == "joint4":
        allres["joint4_primary"] = run(["K2C5_K2SFF", "K2C16_EVEREST", "S46_SPOC2min_PDC", "S72_SPOC2min_PDC"], "joint4_primary", nstep=10000, burn=4000)
        allres["joint4_alt"] = run(["K2C5_EVEREST", "K2C16_SAP", "S46_QLP", "S72_QLP"], "joint4_alt", nstep=10000, burn=4000)
        for k in ["joint4_primary", "joint4_alt"]:
            print(k, json.dumps(allres[k]))
    json.dump(allres, open(path, "w"), indent=1)
