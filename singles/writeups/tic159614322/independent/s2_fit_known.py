"""Joint transit fit to the two known events: shared shape, separate mid-times, local linear baselines.
Primary data: S14 TESS-SPOC (30 min, PDCSAP) and S41 SPOC (2 min, PDCSAP).
Cross-check: same fit with QLP in both sectors.
Shape is parameterised with a nominal reference period P_REF (only the transit shape matters here).
Errors from a self-written adaptive Metropolis MCMC.
"""
import json, os, sys
import numpy as np, batman
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lcio import load, OUT

P_REF = 50.0
U = [0.36, 0.23]  # approximate TESS-band quadratic LD for Teff~5800 K, logg~4.0
WIN = 1.2


def get(sector, product, t0):
    d = load(sector, product, "flux")
    m = np.abs(d.time - t0) < WIN
    t, f = d.time[m].values, d.f[m].values
    oot = np.abs(t - t0) > 0.35
    sig = 1.4826 * np.median(np.abs(f[oot] - np.median(f[oot])))
    return t, f / np.median(f[oot]), np.full(t.size, sig), d.cad.iloc[0]


def make_model(t, cad):
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp, pm.a, pm.inc, pm.ecc, pm.w = 0.0, P_REF, 0.05, 30.0, 89.0, 0.0, 90.0
    pm.u, pm.limb_dark = U, "quadratic"
    if cad > 300:
        mdl = batman.TransitModel(pm, t, supersample_factor=15, exp_time=cad / 86400.0)
    else:
        mdl = batman.TransitModel(pm, t)
    return pm, mdl


def run(prod14, prod41, label, nstep=40000, seed=1):
    rng = np.random.default_rng(seed)
    ev = []
    for (s, p, t0) in [(14, prod14, 1703.83), (41, prod41, 2423.55)]:
        t, f, e, cad = get(s, p, t0)
        ev.append(dict(t=t, f=f, e=e, cad=cad, t0=t0))
    # model setup: we shift times by guess so batman t0 is small
    for d in ev:
        d["pm"], d["mdl"] = None, None

    def model(theta, d, k):
        tc = theta[k]
        rp, lar, b = theta[2], theta[3], theta[4]
        c0, c1 = theta[5 + 2 * k], theta[6 + 2 * k]
        a = np.exp(lar)
        if not (0 <= b < 1 + rp) or a < 1.5:
            return None
        inc = np.degrees(np.arccos(b / a))
        if d["pm"] is None:
            d["pm"], d["mdl"] = make_model(d["t"] - d["t0"], d["cad"])
        pm = d["pm"]
        pm.t0, pm.rp, pm.a, pm.inc = tc - d["t0"], rp, a, inc
        lc = d["mdl"].light_curve(pm)
        return lc * (c0 + c1 * (d["t"] - d["t0"]))

    def resid(theta):
        r = []
        for k, d in enumerate(ev):
            m = model(theta, d, k)
            if m is None:
                return np.full(sum(x["t"].size for x in ev), 1e3)
            r.append((d["f"] - m) / d["e"])
        return np.concatenate(r)

    th0 = np.array([1703.83, 2423.55, 0.055, np.log(45.0), 0.3, 1, 0, 1, 0])
    best = None
    for b0 in [0.1, 0.4, 0.7]:
        th0[4] = b0
        sol = least_squares(resid, th0, x_scale=[1e-3, 1e-3, 1e-3, 0.05, 0.05, 1e-4, 1e-4, 1e-4, 1e-4])
        if best is None or sol.cost < best.cost:
            best = sol
    th = best.x
    chi2 = 2 * best.cost
    ndat = sum(x["t"].size for x in ev)
    # Self-written affine-invariant ensemble sampler (Goodman & Weare 2010 stretch move)
    lnp = lambda x: -0.5 * np.sum(resid(x) ** 2) if (0 <= x[4] < 1 + x[2] and x[2] > 0) else -np.inf
    nw, ndim = 36, th.size
    sc = np.array([5e-4, 5e-4, 1e-3, 0.03, 0.1, 5e-5, 5e-5, 5e-5, 5e-5])
    walkers = th + sc * rng.standard_normal((nw, ndim))
    walkers[:, 4] = np.abs(walkers[:, 4])
    lps = np.array([lnp(w) for w in walkers])
    chain = []
    acc = 0
    a_s = 2.0
    nsteps = nstep // nw
    for it in range(nsteps):
        for half in (0, 1):
            idx = np.arange(nw)[half::2]
            other = walkers[np.arange(nw)[1 - half::2]]
            for i in idx:
                z = ((a_s - 1) * rng.random() + 1) ** 2 / a_s
                y = other[rng.integers(len(other))]
                prop = y + z * (walkers[i] - y)
                lp = lnp(prop)
                if np.log(rng.random()) < (ndim - 1) * np.log(z) + lp - lps[i]:
                    walkers[i], lps[i] = prop, lp
                    acc += 1
        chain.append(walkers.copy())
    nstep = nsteps * nw
    chain = np.array(chain)[nsteps // 3:].reshape(-1, ndim)
    names = ["tc1", "tc2", "rp", "ln_aR", "b", "c0_1", "c1_1", "c0_2", "c1_2"]
    res = {"label": label, "products": [prod14, prod41], "chi2": chi2, "ndata": int(ndat),
           "acceptance": acc / nstep, "best": dict(zip(names, th.tolist()))}
    for n, col in zip(names, chain.T):
        q = np.percentile(col, [16, 50, 84])
        res[n] = [q[1], q[1] - q[0], q[2] - q[1]]
    dT = chain[:, 1] - chain[:, 0]
    q = np.percentile(dT, [16, 50, 84])
    res["dT"] = [q[1], q[1] - q[0], q[2] - q[1]]
    dep = chain[:, 2] ** 2
    q = np.percentile(dep, [16, 50, 84])
    res["depth_rp2"] = [q[1], q[1] - q[0], q[2] - q[1]]
    # total duration T14 for P_REF parameterisation
    a = np.exp(chain[:, 3]); b = chain[:, 4]; rp = chain[:, 2]
    arg = np.sqrt(np.clip((1 + rp) ** 2 - b ** 2, 0, None)) / a / np.sqrt(1 - (b / a) ** 2)
    T14 = P_REF / np.pi * np.arcsin(np.clip(arg, 0, 1)) * 24
    q = np.percentile(T14, [16, 50, 84])
    res["T14_hr"] = [q[1], q[1] - q[0], q[2] - q[1]]
    # implied stellar density at P_REF irrelevant; record rho* for each candidate P later via a/R* ~ P^(2/3)
    # plot
    fig, axs = plt.subplots(1, 2, figsize=(12, 4))
    for k, d in enumerate(ev):
        m = model(th, d, k)
        base = th[5 + 2 * k] + th[6 + 2 * k] * (d["t"] - d["t0"])
        axs[k].plot(d["t"] - th[k], d["f"] / base, ".", ms=2, color="0.5")
        o = np.argsort(d["t"])
        axs[k].plot(d["t"][o] - th[k], (m / base)[o], "r-")
        axs[k].set_title(f"{label}: event {k+1} Tc={th[k]:.4f}")
        axs[k].set_xlabel("t - Tc [d]")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, f"fig_fit_known_{label}.png"), dpi=110)
    plt.close()
    return res, chain


if __name__ == "__main__":
    out = {}
    r, ch = run("TESS-SPOC", "SPOC", "primary", nstep=150000)
    out["primary"] = r
    np.save(os.path.join(OUT, "chain_primary.npy"), ch[::10])
    r2, _ = run("QLP", "QLP", "qlp", nstep=100000, seed=2)
    out["qlp"] = r2
    r3, _ = run("TESS-SPOC", "TESS-SPOC", "tessspoc", nstep=100000, seed=3)
    out["tessspoc"] = r3
    with open(os.path.join(OUT, "fit_known.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for k, v in out.items():
        print(k, "tc1", v["tc1"], "tc2", v["tc2"], "dT", v["dT"], "depth", v["depth_rp2"], "T14", v["T14_hr"], "b", v["b"],
              "chi2/N", v["chi2"], v["ndata"], "acc", v["acceptance"])
