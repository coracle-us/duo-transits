"""
Eclipse-shape fit for TIC 147035057 (second-reviewer item 1, 2026-10-11).

Question: is the 37.315-d eclipse grazing?
Check: with the FLAME star (R = 2.00, M = 1.44), a circular orbit has a/R* = 26.5 and
T14(b=0) = 11.4 h, against an observed ~4.5 h.

Data: the three eclipses (S7 TESS-SPOC 30-min, S33 SPOC 2-min, S108 own TESScut 200-s),
each in a +/-0.6 d window. Each window has its own linear baseline, solved by linear least
squares at every likelihood call (profiled out).

Model:
- batman, circular orbit, P fixed (37.31497 d, error 3e-5 d), T0 free;
- k = Rp/R*, b, rho_circ (stellar density for a circular orbit), quadratic LD via Kipping
  (q1, q2) with uniform priors, one white-noise scale per window.

Fits:
  A  MCMC with a Gaussian prior on rho* from Gaia FLAME (circular orbit)
  B  MCMC with a log-uniform prior on rho_circ (no stellar information)
  A2/B2 as A/B, plus a Gaussian prior on quadratic limb darkening (u1 = 0.24 +/- 0.10,
     u2 = 0.30 +/- 0.10; approximate Claret 2017 TESS values for Teff 6400 K, log g 4.0,
     [M/H] -0.5). Without it, A pushes LD to unphysical extremes to fake a V shape.
  C  maximum likelihood, free rho, constrained to (i) full transit b <= 1 - k
     and (ii) grazing b > 1 - k; delta chi2 / delta BIC
  D  eccentricity needed for a non-grazing solution with the FLAME density:
     (rho_circ / rho*)^(1/3) = (1 + e sin w) / sqrt(1 - e^2)

Output: shape_fit.json, fig_shape.png
Run:    python singles/writeups/tic147035057/shape_fit.py
"""
import json
import sys
from pathlib import Path

import batman
import emcee
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import vet_candidate as vc  # noqa: E402

TIC = 147035057
VET = ROOT / "singles" / "results" / "vetting"
EPH = json.load(open(VET / "147035057_ephemeris.json"))
P = EPH["P"]
G = 6.674e-8
# Gaia DR3 FLAME (astrophysical_parameters): R 1.998 (1.958-2.039), M 1.440 (1.400-1.480)
R_FL, R_FL_E, M_FL, M_FL_E = 1.998, 0.041, 1.440, 0.040
RHO_FL = M_FL * 1.989e33 / (4 / 3 * np.pi * (R_FL * 6.957e10) ** 3)
RHO_FL_E = RHO_FL * np.hypot(M_FL_E / M_FL, 3 * R_FL_E / R_FL)


def load():
    prods = vc.all_products(TIC)
    out = []
    for Tm, sec, pref in zip(EPH["epochs"], (7, 33, 108), ("TESS-SPOC", "SPOC_120", "OWN")):
        p = next(q for q in prods if q["sector"] == sec and q["product"] == pref)
        t = p["time"]
        m = (np.abs(t - Tm) < 0.6) & np.isfinite(p["flux"])
        f = p["flux"][m] / np.nanmedian(p["flux"][m])
        out.append(dict(t=t[m], f=f, dt=t[m] - Tm, ex=float(p["exptime"]), sector=sec, Tm=Tm))
    return out


DATA = load()
MODELS = []
for d in DATA:
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp, pm.a, pm.inc, pm.ecc, pm.w = 0.0, P, 0.05, 26.0, 89.0, 0.0, 90.0
    pm.u, pm.limb_dark = [0.3, 0.2], "quadratic"
    ss = max(1, int(round(d["ex"] / 120)))
    # model time relative to each epoch's own predicted mid-time; a common shift dT0 is fitted
    MODELS.append((pm, d))


def aR_from_rho(rho):
    return (G * rho * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)


def light(theta, d, pm):
    dT0, k, b, lrho, q1, q2 = theta[:6]
    aR = aR_from_rho(10 ** lrho)
    pm.rp, pm.a = k, aR
    pm.inc = np.degrees(np.arccos(np.clip(b / aR, 0, 1)))
    sq = np.sqrt(q1)
    pm.u = [2 * sq * q2, sq * (1 - 2 * q2)]
    ss = max(1, int(round(d["ex"] / 120)))
    m = batman.TransitModel(pm, d["dt"] - dT0, supersample_factor=ss, exp_time=d["ex"] / 86400.0)
    return m.light_curve(pm)


def chi2_epochs(theta, return_models=False):
    tot, mods = 0.0, []
    for i, (pm, d) in enumerate(MODELS):
        lm = light(theta, d, pm)
        A = np.column_stack([lm, lm * d["dt"]])
        c, *_ = np.linalg.lstsq(A, d["f"], rcond=None)
        mod = A @ c
        s = theta[6 + i]
        tot += np.sum(((d["f"] - mod) / s) ** 2) + 2 * len(d["f"]) * np.log(s)
        mods.append(mod)
    return (tot, mods) if return_models else tot


def in_bounds(theta):
    dT0, k, b, lrho, q1, q2 = theta[:6]
    return (abs(dT0) < 0.05 and 0.01 < k < 0.8 and 0 <= b < 1 + k and -3 < lrho < 2
            and 0 < q1 < 1 and 0 < q2 < 1 and all(1e-5 < s < 0.02 for s in theta[6:]))


def lnprob(theta, prior):
    if not in_bounds(theta):
        return -np.inf
    lp = 0.0
    if prior.startswith("flame"):
        lp = -0.5 * ((10 ** theta[3] - RHO_FL) / RHO_FL_E) ** 2
    if prior.endswith("_ld"):
        sq = np.sqrt(theta[4])
        u1, u2 = 2 * sq * theta[5], sq * (1 - 2 * theta[5])
        lp += -0.5 * (((u1 - 0.24) / 0.10) ** 2 + ((u2 - 0.30) / 0.10) ** 2)
    return lp - 0.5 * chi2_epochs(theta)


def sig0():
    out = []
    for d in DATA:
        oot = np.abs(d["dt"]) > 0.15
        out.append(np.std(d["f"][oot] - np.polyval(np.polyfit(d["dt"][oot], d["f"][oot], 1), d["dt"][oot])))
    return out


def run_mcmc(prior, start, nwalk=40, nstep=4000, seed=1):
    rng = np.random.default_rng(seed)
    ndim = len(start)
    p0 = []
    while len(p0) < nwalk:
        trial = np.array(start) + rng.normal(0, 1, ndim) * np.r_[1e-3, 5e-3, 0.02, 0.05, 0.05, 0.05, [2e-5] * 3]
        if np.isfinite(lnprob(trial, prior)):
            p0.append(trial)
    sm = emcee.EnsembleSampler(nwalk, ndim, lnprob, args=(prior,))
    sm.run_mcmc(p0, nstep, progress=False)
    ch = sm.get_chain(discard=nstep // 2, thin=10, flat=True)
    return ch, float(np.mean(sm.acceptance_fraction))


def q(x):
    lo, me, hi = np.percentile(x, [16, 50, 84])
    return dict(median=float(me), minus=float(me - lo), plus=float(hi - me), p2_5=float(np.percentile(x, 2.5)),
                p97_5=float(np.percentile(x, 97.5)))


def summarize(ch):
    k, b, rho = ch[:, 1], ch[:, 2], 10 ** ch[:, 3]
    graz = b > 1 - k
    return dict(k=q(k), b=q(b), rho_circ=q(rho), Rp_Rjup_FLAME=q(k * R_FL / 0.10045),
                frac_grazing=float(graz.mean()), q1=q(ch[:, 4]), q2=q(ch[:, 5]))


def ml_fit(constraint, starts):
    best = None
    for st in starts:
        def f(th):
            if not in_bounds(th):
                return 1e30
            k, b = th[1], th[2]
            if constraint == "full" and b > 1 - k:
                return 1e30
            if constraint == "grazing" and b <= 1 - k:
                return 1e30
            return chi2_epochs(th)
        r = minimize(f, st, method="Nelder-Mead", options=dict(maxiter=40000, maxfev=40000, xatol=1e-7, fatol=1e-4))
        r = minimize(f, r.x, method="Nelder-Mead", options=dict(maxiter=40000, maxfev=40000, xatol=1e-8, fatol=1e-5))
        if best is None or r.fun < best.fun:
            best = r
    th = best.x
    chi2_plain = sum(np.sum(((d["f"] - m) / th[6 + i]) ** 2) for i, (d, m) in enumerate(zip(DATA, chi2_epochs(th, True)[1])))
    return dict(theta=th.tolist(), m2lnL=float(best.fun), k=float(th[1]), b=float(th[2]),
                rho_circ=float(10 ** th[3]), chi2=float(chi2_plain))


def main():
    s0 = sig0()
    res = dict(P=P, rho_FLAME=RHO_FL, rho_FLAME_err=RHO_FL_E,
               aR_FLAME=float(aR_from_rho(RHO_FL)), n_points=[len(d["f"]) for d in DATA])
    # C: maximum-likelihood full vs grazing (no density prior)
    full_starts = [[0, 0.056, b, lr, 0.3, 0.3, *s0] for b in (0.0, 0.5, 0.85) for lr in (0.0, 0.7, 1.2)]
    graz_starts = [[0, k, b, lr, 0.3, 0.3, *s0] for k, b in ((0.08, 1.0), (0.15, 1.05), (0.3, 1.2)) for lr in (-0.6, 0.0, 0.5)]
    full = ml_fit("full", full_starts)
    graz = ml_fit("grazing", graz_starts)
    npts = sum(res["n_points"])
    res["ml_full"] = full
    res["ml_grazing"] = graz
    res["delta_m2lnL_full_minus_grazing"] = full["m2lnL"] - graz["m2lnL"]
    res["delta_BIC_full_minus_grazing"] = full["m2lnL"] - graz["m2lnL"]   # same number of parameters
    # D: eccentricity needed for the best non-grazing (full) solution with the FLAME density
    g = (full["rho_circ"] / RHO_FL) ** (1 / 3)
    res["full_solution_ecc"] = dict(g=float(g), esinw_min_approx=float(g - 1) if g > 1 else 0.0,
                                    e_min_omega90=float((g ** 2 - 1) / (g ** 2 + 1)) if g > 1 else 0.0)
    # B and A: MCMC
    chB, accB = run_mcmc("free", graz["theta"] if graz["m2lnL"] < full["m2lnL"] else full["theta"])
    res["mcmc_free_rho"] = summarize(chB) | dict(acceptance=accB)
    # e needed, over non-grazing samples of B
    ng = chB[chB[:, 2] <= 1 - chB[:, 1]]
    if len(ng):
        gg = (10 ** ng[:, 3] / RHO_FL) ** (1 / 3)
        res["mcmc_free_rho"]["nongrazing_e_min_omega90"] = q(np.where(gg > 1, (gg ** 2 - 1) / (gg ** 2 + 1), 0))
        res["mcmc_free_rho"]["n_nongrazing_samples"] = int(len(ng))
    stA = graz["theta"][:]
    stA[3] = np.log10(RHO_FL)
    chA, accA = run_mcmc("flame", stA, seed=2)
    res["mcmc_flame_prior"] = summarize(chA) | dict(acceptance=accA)
    chB2, accB2 = run_mcmc("free_ld", graz["theta"], seed=3)
    res["mcmc_free_rho_ldprior"] = summarize(chB2) | dict(acceptance=accB2)
    chA2, accA2 = run_mcmc("flame_ld", stA, seed=4)
    res["mcmc_flame_prior_ldprior"] = summarize(chA2) | dict(acceptance=accA2)
    # best fit with FLAME density + LD prior (circular): how bad is it vs the free grazing fit?
    def negpost(th, prior):
        v = lnprob(th, prior)
        return 1e30 if not np.isfinite(v) else -2 * v
    bestA2 = minimize(negpost, chA2[np.argmax([lnprob(x, 'flame_ld') for x in chA2[::50]]) * 50], args=("flame_ld",),
                      method="Nelder-Mead", options=dict(maxiter=20000, maxfev=20000))
    bestB2 = minimize(negpost, chB2[np.argmax([lnprob(x, 'free_ld') for x in chB2[::50]]) * 50], args=("free_ld",),
                      method="Nelder-Mead", options=dict(maxiter=20000, maxfev=20000))
    res["map_flame_ld_vs_free_ld_delta_m2lnpost"] = float(bestA2.fun - bestB2.fun)
    # D2: eccentricity needed for a NON-grazing eclipse with the FLAME star, as a function of b.
    # Observed first-to-last-contact duration from the free grazing ML model.
    th = np.array(graz["theta"])
    tt = np.linspace(-0.4, 0.4, 20001)
    dd = dict(dt=tt, ex=120.0)
    pm = MODELS[1][0]
    lm = light(th, dd, pm)
    inn = tt[lm < 1 - 1e-7]
    T_obs = float((inn.max() - inn.min()) * 24)
    k_ = 0.06
    tab = []
    for b in (0.0, 0.5, 0.7, 0.8, 0.9):
        Tc = P * 24 / np.pi * np.arcsin(np.sqrt((1 + k_) ** 2 - b ** 2) / res["aR_FLAME"])
        gg = Tc / T_obs        # (1 + e sin w) / sqrt(1 - e^2) needed (transit near periastron)
        tab.append(dict(b=b, T14_circ_h=float(Tc), factor=float(gg), e_min_omega90=float((gg ** 2 - 1) / (gg ** 2 + 1)),
                        esinw_if_e_small=float(gg - 1)))
    res["T14_observed_h_grazing_model"] = T_obs
    res["eccentricity_needed_nongrazing"] = tab
    (HERE / "shape_fit.json").write_text(json.dumps(res, indent=1))
    np.savez_compressed(HERE / "shape_fit_chains.npz", free=chB, flame=chA, free_ld=chB2, flame_ld=chA2)
    print(json.dumps({k: v for k, v in res.items() if k not in ("ml_full", "ml_grazing")}, indent=1))
    print("ML full:", {k: v for k, v in full.items() if k != "theta"})
    print("ML grazing:", {k: v for k, v in graz.items() if k != "theta"})

    # figure: data + best full and best grazing models; posterior b vs k
    fig = plt.figure(figsize=(10, 6))
    _, mf = chi2_epochs(np.array(full["theta"]), True)
    _, mg = chi2_epochs(np.array(graz["theta"]), True)
    for i, d in enumerate(DATA):
        ax = fig.add_subplot(2, 3, i + 1)
        o = np.argsort(d["dt"])
        x = d["dt"][o] * 24
        if d["ex"] < 600:
            bins = np.arange(-14.4, 14.41, 0.5)
            idx = np.digitize(x, bins)
            bx = [x[idx == j].mean() for j in np.unique(idx)]
            by = [d["f"][o][idx == j].mean() for j in np.unique(idx)]
            ax.plot(x, (d["f"][o] - 1) * 1e3, ".", ms=1, color="#bcc6d0")
            ax.plot(bx, (np.array(by) - 1) * 1e3, "o", ms=3, color="#1d2b3a")
        else:
            ax.plot(x, (d["f"][o] - 1) * 1e3, "o", ms=3, color="#1d2b3a")
        ax.plot(x, (mf[i][o] - 1) * 1e3, "-", color="#2b8a5e", lw=1.3, label=f"full (b={full['b']:.2f}, k={full['k']:.3f})")
        ax.plot(x, (mg[i][o] - 1) * 1e3, "--", color="#d9682b", lw=1.3, label=f"grazing (b={graz['b']:.2f}, k={graz['k']:.3f})")
        ax.set_xlim(-6, 6)
        ax.set_title(f"S{d['sector']}", fontsize=9)
        ax.set_xlabel("hours from mid-eclipse", fontsize=8)
        if i == 0:
            ax.set_ylabel("relative flux (ppt)", fontsize=8)
            ax.legend(fontsize=6.5, loc="lower left")
    ax = fig.add_subplot(2, 3, 4)
    ax.plot(chB2[:, 1], chB2[:, 2], ",", color="#5d6b7c", alpha=0.3, label="free density (LD prior)")
    ax.plot(chA2[:, 1], chA2[:, 2], ",", color="#d9682b", alpha=0.3, label="FLAME density + LD prior")
    kk = np.linspace(0.01, 0.8, 100)
    ax.plot(kk, 1 - kk, "k--", lw=0.8)
    ax.text(0.35, 0.55, "b = 1 - k\n(grazing above)", fontsize=7)
    ax.set_xlabel("k = R2/R1", fontsize=8)
    ax.set_ylabel("b", fontsize=8)
    ax.set_xlim(0, 0.8)
    ax.set_ylim(0, 1.6)
    ax.legend(fontsize=6.5, markerscale=20)
    ax = fig.add_subplot(2, 3, 5)
    ax.hist(chB2[:, 1] * R_FL / 0.10045, bins=60, color="#5d6b7c", alpha=0.6, density=True, label="free density (LD prior)")
    ax.hist(chA2[:, 1] * R_FL / 0.10045, bins=60, color="#d9682b", alpha=0.6, density=True, label="FLAME + LD prior")
    ax.set_xlabel("R2 [R_Jup] for R1 = 2.0 R_sun", fontsize=8)
    ax.legend(fontsize=6.5)
    ax = fig.add_subplot(2, 3, 6)
    ax.hist(np.log10(10 ** chB2[:, 3]), bins=60, color="#5d6b7c", alpha=0.6, density=True)
    ax.axvline(np.log10(RHO_FL), color="#d9682b")
    ax.set_xlabel("log10 rho_circ [g/cc] (free fit); orange: FLAME", fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "fig_shape.png", dpi=150)


if __name__ == "__main__":
    main()
