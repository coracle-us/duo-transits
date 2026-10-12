"""(a) Aliases n = 4m (m = 4..12) that would also put a transit at the S21 epoch: test their other epochs.
(b) Transit-duration / stellar-density consistency for each alias (circular orbits, rho* from M=1.1, R=1.65).
(c) Linear ephemerides for the surviving aliases using Tc1, the S21 mid-time and Tc2.
"""
import json, os
import numpy as np, batman
from scipy.optimize import least_squares
from model import shape
from lcio import load, OUT

E = json.load(open(os.path.join(OUT, "epochs.json")))
FKALL = json.load(open(os.path.join(OUT, "fit_known.json")))["primary"]
TC1, TC2, DT, T14 = E["TC1"], E["TC2"], E["DT"], E["T14_d"]
S1, S2 = max(FKALL["tc1"][1:]), max(FKALL["tc2"][1:])
SECTORS = [14, 15, 16, 21, 22, 23, 41, 48, 49, 75, 76]
res = {}

# ---------------- (a) high-n multiples of 4 ----------------
data = {}
for s in SECTORS:
    for p in ["SPOC", "TESS-SPOC", "QLP"]:
        d = load(s, p, "flux", strict_quality=True)
        if d is not None:
            data[(s, p)] = (d.time.values, d.f.values, d.cad.iloc[0])


def test_epoch(tc):
    out = []
    for (s, p), (t, f, cad) in data.items():
        m = np.abs(t - tc) < 1.0
        if m.sum() < 20:
            continue
        tt, ff = t[m], f[m] / np.median(f[m])
        S = shape(tt, tc, cad)
        n_exp = T14 * 86400 / cad
        nin = np.sum(S > 0.05)
        if nin < 0.8 * n_exp or np.sum(S == 0) < 10:
            continue
        x = tt - tc
        X = np.vstack([np.ones_like(x), x, x ** 2, -S]).T
        b, *_ = np.linalg.lstsq(X, ff, rcond=None)
        r = ff - X @ b
        rms = np.std(r, ddof=4)
        sig = np.sqrt(np.linalg.inv(X.T @ X)[-1, -1]) * rms
        out.append(dict(sector=s, product=p, depth=float(b[-1]), sig=float(sig)))
    return out


res["high_n"] = {}
print("(a) n = 4m aliases beyond 12:")
for n in range(16, 52, 4):
    P = DT / n
    verdict, tested = "allowed", []
    for k in range(-3, int((3400 - TC1) / P) + 2):
        tc = TC1 + k * P
        if k in (0, n, n // 4):
            continue
        tr = test_epoch(tc)
        if tr:
            tested.append(dict(k=k, tc=tc, results=tr))
            # exclude if a fully-covered product shows depth < 1 ppt and is >5 sigma below 3 ppt
            if any(r["depth"] < 0.001 and (0.003 - r["depth"]) / r["sig"] > 5 for r in tr):
                verdict = "excluded"
    n_excl = sum(any(r["depth"] < 0.001 and (0.003 - r["depth"]) / r["sig"] > 5 for r in e["results"]) for e in tested)
    res["high_n"][n] = dict(P=P, verdict=verdict, n_epochs_full_coverage=len(tested), n_epochs_flat=int(n_excl),
                            epochs=tested)
    print(f"   n={n:2d} P={P:7.3f} d: {len(tested)} fully covered epochs, {n_excl} flat -> {verdict}")

# ---------------- (b) duration / density consistency ----------------
G = 6.674e-8
Msun, Rsun = 1.989e33, 6.957e10
M, R = 1.1, 1.65
rho = M * Msun / (4 / 3 * np.pi * (R * Rsun) ** 3)
ev = []
for s, p, t0 in [(14, "TESS-SPOC", TC1), (41, "SPOC", TC2)]:
    d = load(s, p, "flux")
    m = np.abs(d.time - t0) < 1.2
    t, f = d.time[m].values, d.f[m].values
    oot = np.abs(t - t0) > 0.35
    f = f / np.median(f[oot])
    e = 1.4826 * np.median(np.abs(f[oot] - 1))
    ev.append((t, f, e, d.cad.iloc[0], t0))


def chi2_for(P, aR_fixed=None):
    def resid(th):
        rp, b = th[0], th[1]
        aR = aR_fixed if aR_fixed is not None else np.exp(th[2])
        rr = []
        for k, (t, f, e, cad, t0) in enumerate(ev):
            pm = batman.TransitParams()
            pm.t0, pm.per, pm.rp, pm.a = t0 + th[3 + k], P, rp, aR
            pm.inc = np.degrees(np.arccos(min(b / aR, 1)))
            pm.ecc, pm.w, pm.u, pm.limb_dark = 0, 90, [0.36, 0.23], "quadratic"
            mdl = batman.TransitModel(pm, t, supersample_factor=15, exp_time=cad / 86400) if cad > 300 else batman.TransitModel(pm, t)
            c0, c1 = th[5 + 2 * k], th[6 + 2 * k]
            rr.append((f - mdl.light_curve(pm) * (c0 + c1 * (t - t0))) / e)
        return np.concatenate(rr)
    best = None
    for b0 in [0.1, 0.5, 0.8, 0.9, 0.95]:
        th0 = [0.055, b0, np.log(40.0), 0, 0, 1, 0, 1, 0]
        lb = [0.01, 0, np.log(2), -0.05, -0.05, 0.9, -1, 0.9, -1]
        ub = [0.3, 1.3, np.log(2000), 0.05, 0.05, 1.1, 1, 1.1, 1]
        sol = least_squares(resid, th0, bounds=(lb, ub), x_scale=[1e-3, 0.05, 0.05, 1e-3, 1e-3, 1e-4, 1e-4, 1e-4, 1e-4])
        if best is None or sol.cost < best.cost:
            best = sol
    return 2 * best.cost, best.x


chi_free, x_free = chi2_for(50.0, None)
res["density"] = dict(rho_star_cgs=rho, chi2_free_shape=chi_free, aR_free_at_P50=float(np.exp(x_free[2])), per_alias={})
print(f"\n(b) rho* = {rho:.3f} g/cc; free-shape chi2 = {chi_free:.1f}")
for n in range(1, 13):
    P = DT / n
    aR = (G * rho * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    c, x = chi2_for(P, aR)
    res["density"]["per_alias"][n] = dict(P=P, aR_circ=aR, chi2=c, dchi2=c - chi_free, b=float(x[1]), rp=float(x[0]))
    print(f"   n={n:2d} P={P:8.3f}  a/R*={aR:6.1f}  best b={x[1]:.3f} rp={x[0]:.4f}  Δχ²={c - chi_free:6.1f}")

# ---------------- (c) ephemerides ----------------
ing = json.load(open(os.path.join(OUT, "ingress_fap.json")))["ingress_timing"]
vals = np.array([v["tc"] for k, v in ing.items() if "free depth" in k])
tc3 = float(np.median(vals))
formal = np.median([(v["hi"] - v["lo"]) / 2 for k, v in ing.items() if "free depth" in k])
spread = 1.4826 * np.median(np.abs(vals - tc3))
s3 = float(np.hypot(formal, spread))
res["tc_S21"] = dict(value=tc3, sigma=s3, formal=float(formal), product_spread=float(spread), all=vals.tolist())
pred = TC1 + DT / 4
print(f"\n(c) S21 mid-time (median over products, ingress-side fits) = {tc3:.4f} ± {s3:.4f} "
      f"(formal {formal:.4f}, spread {spread:.4f}); predicted {pred:.4f}; O-C = {(tc3-pred)*24:+.2f} h = {(tc3-pred)/np.hypot(s3, 0.003):.1f} sigma")
res["ephemeris"] = {}
for n in [4, 8]:
    E_ = np.array([0, n // 4, n])
    T = np.array([TC1, tc3, TC2])
    W = 1 / np.array([S1, s3, S2]) ** 2
    A = np.vstack([np.ones(3), E_]).T
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    beta = cov @ A.T @ (W * T)
    chi = np.sum(W * (T - A @ beta) ** 2)
    # reference epoch near weighted centre to decorrelate
    res["ephemeris"][n] = dict(P=float(beta[1]), sig_P=float(np.sqrt(cov[1, 1])), T0=float(beta[0]), sig_T0=float(np.sqrt(cov[0, 0])),
                               chi2=float(chi), epochs=E_.tolist())
    print(f"   n={n}: P = {beta[1]:.5f} ± {np.sqrt(cov[1,1]):.5f} d, T0 = {beta[0]:.4f} ± {np.sqrt(cov[0,0]):.4f} BTJD, chi2={chi:.2f}")
json.dump(res, open(os.path.join(OUT, "highn_density_ephem.json"), "w"), indent=1, default=float)
