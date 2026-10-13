"""Step 11: K2 Campaign 0 (2014) -- inclination from the single primary eclipse, and the
secondary-eclipse limit.

Shape parameters (esw, r1, r2, J_TESS) are taken from the TESS posterior; J is converted to the
Kepler band with blackbody surface-brightness ratios (top-hat bandpasses), T1 = 7754 K (TIC).
Free per K2 fit: i2014, primary time shift, linear baseline.  e cos w for 2014 = linear
extrapolation of the 2021->2023 values (it only enters through e; tested +-0.03).
Third light in the K2 apertures: CROWDSAP = 0.9989 (PDCSAP already corrected) -> l3 = 0.
Outputs: k2_results.json, fig_k2.png
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ebmodel as M

H = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(H, "data")
tag = sys.argv[1] if len(sys.argv) > 1 else "nominal"
pj = json.load(open(os.path.join(H, "fit_tess_%s.json" % tag)))
chp = os.path.join(D, "chain_%s.npz" % tag)
HAVE_CHAIN = os.path.exists(chp)
if HAVE_CHAIN:
    ch = np.load(chp)
    names = list(ch["names"]); flat = ch["flat"]
UK = (0.28, 0.30)          # Kepler-band quadratic LD, ~7500-7800 K, logg 3.9 (approximate)
T1 = 7754.0


def bb_band(T, lo, hi):
    lam = np.linspace(lo, hi, 400) * 1e-9
    x = 1.4388e-2 / (lam * T)
    return np.trapezoid(lam ** -4 / np.expm1(x), lam)   # photon-weighted (lambda^-5 * lambda)


_JK = {}
def j_kepler(J_tess):
    key = round(float(J_tess), 6)
    if key not in _JK:
        _JK[key] = _j_kepler(J_tess)
    return _JK[key]


def _j_kepler(J_tess):
    f = lambda T2: bb_band(T2, 600, 1000) / bb_band(T1, 600, 1000) - J_tess
    T2 = brentq(f, 2500, 30000)
    return bb_band(T2, 420, 900) / bb_band(T1, 420, 900), T2


W = pd.read_csv(os.path.join(D, "windows.csv"))
W = W[W.src.isin(["k2pdc", "k2sff", "everest"]) & (W.season == 2014) & (W.epoch.isin([-124, -125]))]
ecw14 = pj["ecw2021"] - (pj["ecw2023"] - pj["ecw2021"]) * (2510.0 + 211.6) / (3250.0 - 2510.0)
P, T0 = pj["P"], pj["T0"]


def prep(src, epoch, typ, half):
    d = W[(W.src == src) & (W.epoch == epoch) & (W.type == typ)].sort_values("time")
    tc = T0 + epoch * P + (pj["_dsec2021"] if typ == "S" else 0)
    d = d[np.abs(d.time - tc) < half]
    t, y = d.time.values, d.flux.values
    oot = np.abs(t - tc) > 0.3
    y = y / np.median(y[oot])
    # 5-sigma clip of out-of-eclipse outliers only
    s = 1.4826 * np.median(np.abs(y[oot] - np.median(y[oot])))
    keep = ~(oot & (np.abs(y - 1) > 5 * s))
    return t[keep], y[keep], tc, s


def model(t, inc, dt, sh, ecw):
    esw, r1, r2, J = sh
    e = np.hypot(ecw, esw); w = np.arctan2(esw, ecw)
    Jk, _ = j_kepler(J)
    return M.lc(t - dt, P, T0, e, w, inc, r1, r2, Jk, UK[0], UK[1], UK[0], UK[1], l3=0.0,
                texp=1765.5 / 86400, nsub=15)


def fit_primary(t, y, tc, sig, sh, ecw):
    def res(x):
        m = model(t, x[0], x[1], sh, ecw)
        A = np.vstack([m, m * (t - tc)]).T
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        return (y - A @ c) / sig
    best = None
    for i0 in np.arange(84.0, 88.01, 0.5):
        for dt0 in (-0.03, 0.0, 0.03):
            try:
                r = least_squares(res, [i0, dt0], bounds=([75, -0.2], [90, 0.2]), x_scale=[0.1, 0.01])
            except Exception:
                continue
            if best is None or r.cost < best.cost:
                best = r
    return best


# need dsec for 2021 (conjunction offset) for the secondary window centre
sh0 = (pj["esw"], pj["r1"], pj["r2"], pj["J"])
e21 = np.hypot(pj["ecw2021"], pj["esw"]); w21 = np.arctan2(pj["esw"], pj["ecw2021"])
pj["_dsec2021"] = M.secondary_offset(P, e21, w21)

out = {"ecw2014_assumed": ecw14, "UK": UK}
print("ecw2014 (extrapolated) = %.4f" % ecw14)
Jk0, T2 = j_kepler(pj["J"])
print("J_TESS=%.3f -> T2=%.0f K -> J_Kp=%.3f" % (pj["J"], T2, Jk0))
out["T2_K"] = T2; out["J_Kp_nominal"] = Jk0

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
res_src = {}
for src in ["k2pdc", "k2sff"]:   # EVEREST quality flags remove the in-eclipse cadences
    t, y, tc, sig = prep(src, -124, "P", 1.0)
    r = fit_primary(t, y, tc, sig, sh0, ecw14)
    m = model(t, r.x[0], r.x[1], sh0, ecw14)
    A = np.vstack([m, m * (t - tc)]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    chi2 = np.sum(((y - A @ c) / sig) ** 2)
    tt = np.linspace(t.min(), t.max(), 3000)
    mm = model(tt, r.x[0], r.x[1], sh0, ecw14)
    depth = 1 - mm.min()
    J_ = np.linalg.pinv(r.jac.T @ r.jac) * max(1, chi2 / (len(t) - 4))
    res_src[src] = dict(i2014=r.x[0], i_err=np.sqrt(J_[0, 0]), dt=r.x[1], dt_err=np.sqrt(J_[1, 1]),
                        tmin=tt[np.argmin(mm)], depth_model=depth, depth_obs_min=1 - np.min(y / (A @ c / m)),
                        n=len(t), n_in=int(np.sum(m < 0.999)), chi2=chi2, oot_rms=sig)
    print(src, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in res_src[src].items()})
    axes[0].plot(t - tc, y / (A @ c / m), ".", label=src)
    axes[0].plot(tt - tc, mm, "-", lw=0.8)
axes[0].set_xlim(-0.6, 0.6); axes[0].legend(); axes[0].set_title("K2 primary (epoch -124)")
out["primary_by_pipeline"] = res_src

# posterior propagation (k2sff = cleanest; also pdc)
rng = np.random.default_rng(1)
if HAVE_CHAIN:
    idx = rng.choice(len(flat), 300, replace=False)
    k = {n: names.index(n) for n in ["esw", "r1", "r2", "J"]}
for src in (["k2sff", "k2pdc"] if HAVE_CHAIN else []):
    t, y, tc, sig = prep(src, -124, "P", 1.0)
    iis = []
    for j in idx:
        sh = tuple(flat[j, k[n]] for n in ["esw", "r1", "r2", "J"])
        ecw = ecw14 + rng.normal(0, 0.02)
        r = fit_primary(t, y, tc, sig, sh, ecw)
        jac = r.jac
        cov = np.linalg.pinv(jac.T @ jac) * max(1, 2 * r.cost / (len(t) - 4))
        iis.append(r.x[0] + rng.normal(0, np.sqrt(cov[0, 0])))
    iis = np.array(iis)
    q = np.percentile(iis, [16, 50, 84])
    out["i2014_post_" + src] = dict(p16=q[0], p50=q[1], p84=q[2])
    print("i2014 posterior (%s): %.3f +%.3f -%.3f" % (src, q[1], q[2] - q[1], q[1] - q[0]))
    np.save(os.path.join(D, "i2014_samples_%s_%s.npy" % (src, tag)), iis)

# LD and ecw systematics for k2sff
t, y, tc, sig = prep("k2sff", -124, "P", 1.0)
syst = {}
for lab, uk, ecw in [("UK-0.1", (UK[0] - 0.1, UK[1]), ecw14), ("UK+0.1", (UK[0] + 0.1, UK[1]), ecw14),
                     ("u2-0.1", (UK[0], UK[1] - 0.1), ecw14), ("u2+0.1", (UK[0], UK[1] + 0.1), ecw14),
                     ("ecw-0.03", UK, ecw14 - 0.03), ("ecw+0.03", UK, ecw14 + 0.03),
                     ("ecw=2021", UK, pj["ecw2021"]), ("apsidal-extrap", UK, "aps")]:
    UK_save = UK
    UK = uk
    sh_use = sh0
    if ecw == "aps":
        # w advancing at the 2021->2023 rate of the esw-per-season fit, extrapolated back to 2014
        pe = json.load(open(os.path.join(H, "fit_tess_esw_season.json")))
        w21 = np.arctan2(pe["esw2021"], pe["ecw2021"]); w23 = np.arctan2(pe["esw2023"], pe["ecw2023"])
        e21 = np.hypot(pe["esw2021"], pe["ecw2021"])
        w14 = w21 - (w23 - w21) * (2510.0 + 211.6) / (3258.0 - 2510.0)
        ecw = e21 * np.cos(w14); sh_use = (e21 * np.sin(w14),) + tuple(sh0[1:])
        print("apsidal extrapolation: w2021=%.2f w2023=%.2f -> w2014=%.2f deg, ecw=%.4f esw=%.4f"
              % (np.degrees(w21), np.degrees(w23), np.degrees(w14), ecw, sh_use[0]))
    r = fit_primary(t, y, tc, sig, sh_use, ecw)
    UK = UK_save
    syst[lab] = r.x[0]
    print("syst %-9s i2014=%.3f" % (lab, r.x[0]))
out["i2014_syst"] = syst

# secondary-eclipse limit: inject the model secondary for a grid of i and test
sec_res = []
for src in ["k2sff", "k2pdc"]:
    for ep in (-125, -124):
        ts, ys, tcs, sigs = prep(src, ep, "S", 1.0)
        sec_res.append((src, ep, ts, ys, tcs, sigs))
        axes[1].plot(ts - tcs, ys, ".", ms=3, label="%s %d" % (src, ep))
i_grid = np.arange(84.0, 88.01, 0.05)
chi_grid = []
for inc in i_grid:
    c2 = 0
    for (src, ep, ts, ys, tcs, sigs) in sec_res:
        if src != "k2sff":
            continue
        best = np.inf
        for dt in res_src["k2sff"]["dt"] + np.arange(-0.10, 0.101, 0.01):     # allow timing freedom
            m = model(ts, inc, dt, sh0, ecw14)
            A = np.vstack([m, m * (ts - tcs)]).T
            c, *_ = np.linalg.lstsq(A, ys, rcond=None)
            best = min(best, np.sum(((ys - A @ c) / sigs) ** 2))
        c2 += best
    chi_grid.append(c2)
chi_grid = np.array(chi_grid)
dchi = chi_grid - chi_grid.min()
i_sec_max = i_grid[np.where(dchi < 9)[0].max()] if np.any(dchi < 9) else np.nan
# predicted secondary depth vs i
dep = []
for inc in i_grid:
    tt = tcs + np.linspace(-0.6, 0.6, 2000)
    dep.append(1 - model(tt, inc, res_src['k2sff']['dt'], sh0, ecw14).min())
dep = np.array(dep)
first_sec = i_grid[np.argmax(dep > 1e-5)]
print("secondary: Delta chi2 < 9 for i <= %.2f ; model secondary appears for i > %.2f" % (i_sec_max, first_sec))
# empirical depth limit: box of 0.3 d width at any time within +-0.15 d of predicted centre
lims = []
for (src, ep, ts, ys, tcs, sigs) in sec_res:
    if src != "k2sff":
        continue
    oot = np.abs(ts - tcs) > 0.45
    for c0 in np.arange(-0.15, 0.151, 0.02):
        inn = np.abs(ts - tcs - c0) < 0.1
        lims.append((np.median(ys[oot]) - np.mean(ys[inn])) / (sigs / np.sqrt(inn.sum())))
print("box-depth S/N in secondary windows: min %.2f max %.2f" % (np.min(lims), np.max(lims)))
sec_lim = []
for (src, ep, ts, ys, tcs, sigs) in sec_res:
    if src != "k2sff":
        continue
    inn = np.abs(ts - tcs) < 0.15
    sec_lim.append(3 * sigs / np.sqrt(inn.sum()))
out["secondary"] = dict(i_max_dchi9=float(i_sec_max), i_first_secondary=float(first_sec),
                        depth_3sigma_limit_per_window=[float(x) for x in sec_lim],
                        depth_mean_signal_over_noise_range=[float(np.min(lims)), float(np.max(lims))])
axes[1].set_title("K2 secondary windows"); axes[1].legend(fontsize=7)
axes[2].plot(i_grid, dchi, "k-"); axes[2].axhline(9, color="r", ls=":")
axes[2].plot(i_grid, dep * 1e3 * 10, "b-", label="model sec depth x10 (ppt)")
axes[2].set_xlabel("i2014 (deg)"); axes[2].set_ylabel("Delta chi2 (secondary windows)"); axes[2].legend()
axes[2].set_ylim(0, 100)
plt.tight_layout(); plt.savefig(os.path.join(H, "fig_k2_%s.png" % tag), dpi=80)
json.dump(out, open(os.path.join(H, "k2_results_%s.json" % tag), "w"), indent=1, default=float)
