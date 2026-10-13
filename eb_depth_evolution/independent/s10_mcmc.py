"""Step 10: emcee posterior for the TESS global fit (statistical uncertainties).

Per-eclipse errors are rescaled to reduced chi2 = 1 and then multiplied by a red-noise factor
beta (time-averaging method, bins of 20-60 min).  Free: esw, r1, r2, J, i2021, i2023, i2026 and
one time shift per eclipse.  e cos w per season fixed at the LS values (set by the secondary
intervals, which the free shifts absorb anyway).
usage: s10_mcmc.py [fitjson] [nsteps]
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import emcee
import ebdata, ebfit

H = os.path.dirname(os.path.abspath(__file__))
fj = sys.argv[1] if len(sys.argv) > 1 else "fit_tess_nominal.json"
nsteps = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
pj = json.load(open(os.path.join(H, fj)))
pars = {k: v for k, v in pj.items() if not k.startswith("_")}
df = ebdata.load(pj["_src"], col=pj["_col"], seasons=[2021, 2023, 2026])
df = df[df.ecl != "S72"].reset_index(drop=True)
F0 = ebfit.Fitter(df, pars, [], **pj["_fitter_kw"])

# error rescaling per eclipse
beta_tab = {}
for g in F0.groups:
    m = F0.model_group(F0.pars, g)
    c = F0.baseline(m, g)
    res = g["y"] - m * (c[0] + c[1] * (g["t"] - g["tc"]))
    s1 = np.std(res)
    cad = np.median(np.diff(g["t"]))
    betas = []
    for binw in [20, 30, 45, 60]:
        nb = max(int(round(binw / 1440 / cad)), 2)
        nbins = len(res) // nb
        if nbins < 5:
            continue
        rb = res[:nbins * nb].reshape(nbins, nb).mean(axis=1)
        expect = s1 / np.sqrt(nb) * np.sqrt(nbins / (nbins - 1))
        betas.append(np.std(rb) / expect)
    beta = max(1.0, np.median(betas)) if betas else 1.0
    beta_tab[g["ecl"]] = (s1, beta)
    df.loc[df.ecl == g["ecl"], "err"] = s1 * beta
    print("%-4s rms=%.2f ppt beta=%.2f" % (g["ecl"], s1 * 1e3, beta))

present = sorted(df.season.unique())
free = ["esw", "r1", "r2", "J"] + ["i%d" % s for s in present] + ["dt_" + g["ecl"] for g in F0.groups]
F = ebfit.Fitter(df, pars, free, **pj["_fitter_kw"])
r = F.fit(max_nfev=3000, xtol=1e-12, ftol=1e-12)
x0 = np.array([F.pars[n] for n in free])
ndim, nw = len(free), 48
rng = np.random.default_rng(42)
scale = np.array([1e-3 if n in ("r1", "r2") else 1e-2 if n in ("esw", "J") else 0.02 if n.startswith("i")
                  else 3e-4 for n in free])
p0 = x0 + scale * 0.1 * rng.standard_normal((nw, ndim))
p0 = np.clip(p0, F.lo + 1e-9, F.hi - 1e-9)
sampler = emcee.EnsembleSampler(nw, ndim, F.lnprob,
                                moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
sampler.run_mcmc(p0, nsteps, progress=False)
chain = sampler.get_chain()
lnp = sampler.get_log_prob()
try:
    tau = sampler.get_autocorr_time(tol=0)
    print("autocorr times:", dict(zip(free, np.round(tau, 0))))
except Exception as ex:
    print("tau fail", ex)
burn = nsteps // 2
flat = chain[burn:].reshape(-1, ndim)
print("acceptance", np.mean(sampler.acceptance_fraction))
for k, n in enumerate(free):
    q = np.percentile(flat[:, k], [16, 50, 84])
    print("%-8s %.6f +%.6f -%.6f" % (n, q[1], q[2] - q[1], q[1] - q[0]))
tag = fj.replace("fit_tess_", "").replace(".json", "")
np.savez(os.path.join(H, "data", "chain_%s.npz" % tag), flat=flat[::10], names=np.array(free),
         lnp=lnp[burn:].ravel()[::10], fixed=json.dumps({k: v for k, v in F.pars.items() if k not in free}),
         errs=json.dumps({k: list(v) for k, v in beta_tab.items()}))
