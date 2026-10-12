"""EPIC 211892898 / TIC 184895645: independent eclipse timing from the cached K2 light curves,
period fit, alias check against K2+TESS coverage, and comparison with the Gaia DR3 SB1 orbit.

Run from the repo root with the venv python. Output: ecl_211892898.txt, ecl_211892898.png."""
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EPIC = int(sys.argv[1]) if len(sys.argv) > 1 else 211892898
TIC = {211892898: 184895645, 211390677: 366568535}[EPIC]
# catalogue event + pipeline repeats (grid times, step D/4) to seed the fits
SEEDS = {211892898: [189.3519, 1132.1486, 1276.8096], 211390677: [151.7983, 1257.5695]}[EPIC]
D_CAT = {211892898: 62.2756 / 24, 211390677: 39.7195 / 24}[EPIC]
out = open(HERE / f"ecl_{EPIC}.txt", "w")


def say(*a):
    print(*a)
    print(*a, file=out)


def load_k2():
    parts = []
    for fn in sorted(glob.glob(str(ROOT / "k2singles" / "cache" / str(EPIC) / "k2_c*.npz"))):
        z = np.load(fn)
        c = int(Path(fn).stem.split("_c")[1])
        parts.append(pd.DataFrame(dict(t=z["t"], f=z["f"], camp=c)))
    return pd.concat(parts, ignore_index=True).sort_values("t").reset_index(drop=True)


def trap(t, t0, depth, T14, Tin):
    """Trapezoid: total duration T14, ingress/egress Tin (<= T14/2)."""
    Tin = np.clip(Tin, 1e-3, T14 / 2)
    x = np.abs(t - t0)
    m = np.zeros_like(t)
    full = x <= T14 / 2 - Tin
    ramp = (x > T14 / 2 - Tin) & (x < T14 / 2)
    m[full] = depth
    m[ramp] = depth * (T14 / 2 - x[ramp]) / Tin
    return 1 - m


def fit_one(t, f, t_seed, D, win=4.0):
    """Fit trapezoid x quadratic baseline to data within +-(D/2+win) of the seed. Returns dict."""
    sel = np.abs(t - t_seed) < D / 2 + win
    tt, ff = t[sel], f[sel]
    tm = tt.mean()

    def model(p, tt):
        t0, depth, T14, Tin, c0, c1, c2 = p
        return trap(tt, t0, depth, T14, Tin) * (c0 + c1 * (tt - tm) + c2 * (tt - tm) ** 2)

    def resid(p):
        return ff - model(p, tt)

    p0 = [t_seed, 0.008, D, 0.2 * D, 1.0, 0.0, 0.0]
    lo = [t_seed - 1.5, 0, 0.3 * D, 0.01, 0.9, -1, -1]
    hi = [t_seed + 1.5, 0.1, 2 * D, D, 1.1, 1, 1]
    r = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac")
    # refit with a sigma-clipped set
    res = resid(r.x)
    sig = 1.4826 * np.median(np.abs(res - np.median(res)))
    keep = np.abs(res) < 4 * sig
    r = least_squares(lambda p: (ff - model(p, tt))[keep], r.x, bounds=(lo, hi), x_scale="jac")
    res = (ff - model(r.x, tt))[keep]
    sig = np.std(res, ddof=len(p0))
    # residual bootstrap for the t0 error (cyclic shift keeps the red-noise structure)
    rng = np.random.default_rng(1)
    boots = []
    mod = model(r.x, tt)
    for _ in range(300):
        shift = rng.integers(len(res))
        fb = mod[keep] + np.roll(res, shift)
        rb = least_squares(lambda p: fb - model(p, tt[keep]), r.x, bounds=(lo, hi), x_scale="jac")
        boots.append(rb.x[:4])
    boots = np.array(boots)
    return dict(t0=r.x[0], t0_err=boots[:, 0].std(), depth=r.x[1], depth_err=boots[:, 1].std(), T14=r.x[2],
                T14_err=boots[:, 2].std(), Tin=r.x[3], sig=sig, n=int(keep.sum()), tt=tt, ff=ff, mod=mod, keep=keep)


k2 = load_k2()
say(f"EPIC {EPIC}: K2 campaigns {sorted(k2.camp.unique())}")
for c in sorted(k2.camp.unique()):
    s = k2[k2.camp == c]
    say(f"  C{c}: BTJD {s.t.min():.3f}-{s.t.max():.3f}, {len(s)} cadences")

fits = []
for T in SEEDS:
    r = fit_one(k2.t.to_numpy(), k2.f.to_numpy(), T, D_CAT)
    fits.append(r)
    say(f"seed {T:.3f}: t0 = {r['t0']:.4f} +- {r['t0_err']:.4f}  depth = {r['depth']*1e3:.2f} +- {r['depth_err']*1e3:.2f} ppt"
        f"  T14 = {r['T14']*24:.1f} +- {r['T14_err']*24:.1f} h  Tin = {r['Tin']*24:.1f} h  rms {r['sig']*1e3:.2f} ppt  n={r['n']}")

# ---- plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, len(fits), figsize=(5 * len(fits), 3.5))
for ax, r in zip(np.atleast_1d(axs), fits):
    ax.plot(r["tt"] - r["t0"], r["ff"], ".", ms=2, color="0.5")
    ax.plot(r["tt"] - r["t0"], r["mod"], "r-", lw=1)
    ax.set_title(f"t0={r['t0']:.3f}+-{r['t0_err']:.3f}")
    ax.set_xlabel("t - t0 (d)")
fig.tight_layout()
fig.savefig(HERE / f"ecl_{EPIC}.png", dpi=110)

t0s = np.array([r["t0"] for r in fits])
e0s = np.array([r["t0_err"] for r in fits])
say("\nEclipse mid-times (BTJD):", np.round(t0s, 4), "+-", np.round(e0s, 4))

# ---- period fit over all integer cycle assignments consistent with pairwise separations
say("\nLinear-ephemeris fits P = dT/n (anchored on the first two times), chi2 of the third:")
dT = t0s[1] - t0s[0]
rows = []
for n in range(1, int(dT / 5) + 1):
    P = dT / n
    cyc = np.round((t0s - t0s[0]) / P)
    if len(set(cyc)) < len(cyc):
        continue
    A = np.vstack([np.ones(len(t0s)), cyc]).T
    W = 1 / e0s
    coef, *_ = np.linalg.lstsq(A * W[:, None], t0s * W, rcond=None)
    T0, Pf = coef
    resid = t0s - (T0 + cyc * Pf)
    chi2 = np.sum((resid / e0s) ** 2)
    cov = np.linalg.inv((A * W[:, None]).T @ (A * W[:, None]))
    rows.append(dict(n=n, P=Pf, P_err=np.sqrt(cov[1, 1]), T0=T0, chi2=chi2, max_oc_d=np.abs(resid).max()))
pf = pd.DataFrame(rows)
say(pf[pf.chi2 < 50].to_string(index=False))
pf.to_csv(HERE / f"ecl_{EPIC}_periods.csv", index=False)

# ---- Coverage test of each surviving period against ALL K2 + TESS data (by eye-proof box test)
# For each candidate P (chi2 < 25) count predicted eclipses that fall fully inside K2/TESS data with
# a flat light curve: in-eclipse mean within 3 sigma of zero and depth < 0.5*observed.
sys.path.insert(0, str(ROOT / "scripts"))
import alias_pipeline as ap  # noqa: E402

tess = ap.load_lightcurves(TIC, products=ap.PRODUCT_RANK + ["OWN"])
say(f"\nTESS sectors loaded for TIC {TIC}: {sorted(tess.sector.unique())}")
k2lc = pd.DataFrame(dict(time=k2.t, flux=k2.f, sector=1000 + k2.camp, segment=1000 + k2.camp))
allc = pd.concat([tess[["time", "flux", "sector", "segment"]], k2lc]).sort_values("time")
t_all, f_all = allc.time.to_numpy(), allc.flux.to_numpy()
depth_obs = np.median([r["depth"] for r in fits])
D_fit = np.median([r["T14"] for r in fits])


def window_test(Tc, D, delta):
    """Simple: data within +-0.35 D of Tc (flat part) and baseline in [0.5D, 0.5D+1.5d] each side."""
    inn = np.abs(t_all - Tc) < 0.35 * D
    oot = (np.abs(t_all - Tc) > 0.5 * D + 0.1) & (np.abs(t_all - Tc) < 0.5 * D + 1.5)
    if inn.sum() < 10 or (oot & (t_all < Tc)).sum() < 10 or (oot & (t_all > Tc)).sum() < 10:
        return "nodata", np.nan, np.nan
    c = np.polyfit(t_all[oot] - Tc, f_all[oot], 1)
    r = f_all / np.polyval(c, t_all - Tc) - 1
    sig = 1.4826 * np.median(np.abs(r[oot] - np.median(r[oot])))
    d = -np.mean(r[inn])
    s = sig / np.sqrt(inn.sum()) * 2.0  # x2 red-noise allowance
    if d > 0.5 * delta and d / s > 5:
        return "DIP", d, s
    if d + 3 * s < 0.5 * delta:
        return "flat", d, s
    return "ambig", d, s


say(f"\nAlias check (observed depth {depth_obs*1e3:.2f} ppt, T14 {D_fit*24:.1f} h) against K2 C5/C16/C18 + TESS:")
for _, row in pf[pf.chi2 < 25].iterrows():
    P, T0 = row.P, row.T0
    ks = np.arange(np.floor((t_all.min() - T0) / P), np.ceil((t_all.max() - T0) / P) + 1)
    verd = []
    for k in ks:
        Tc = T0 + k * P
        if Tc < t_all.min() - 1 or Tc > t_all.max() + 1:
            continue
        v, d, s = window_test(Tc, D_fit, depth_obs)
        if v != "nodata":
            sec = allc.sector.to_numpy()[np.searchsorted(t_all, Tc) - 1]
            verd.append((int(k), round(Tc, 2), int(sec), v, None if np.isnan(d) else round(d * 1e3, 2)))
    nflat = sum(1 for v in verd if v[3] == "flat")
    ndip = sum(1 for v in verd if v[3] == "DIP")
    say(f"  n={int(row.n):3d} P={P:9.4f}: {ndip} dips, {nflat} FLAT windows, {len(verd)-ndip-nflat} ambiguous -> "
        f"{'EXCLUDED' if nflat else 'allowed'}")
    if nflat and nflat <= 12:
        say("        flat at:", [v for v in verd if v[3] == "flat"])
    say("        dips at:", [v for v in verd if v[3] == "DIP"])

# ---- Gaia DR3 SB1 orbit -> predicted conjunction times
nss = pd.read_csv(HERE / "gaia_nss.csv")
nss = nss[nss.source_id == {211892898: 660920219027943552}.get(EPIC, -1)]
if len(nss):
    g = nss.iloc[0]
    say("\nGaia DR3 NSS SB1:", f"P={g.period:.4f}+-{g.period_error:.4f} d, Tp={g.t_periastron:.3f}+-{g.t_periastron_error:.3f} (BTJD {389.0+g.t_periastron:.3f}),"
        f" e={g.eccentricity:.3f}+-{g.eccentricity_error:.3f}, omega={g.arg_periastron:.2f}+-{g.arg_periastron_error:.2f} deg, K={g.semi_amplitude_primary:.2f} km/s")

    def conj_times(P, Tp, e, w_deg, which):
        """Time (within the first orbit after Tp) when nu + omega = 90 deg (which=0: primary at
        superior conjunction = PRIMARY ECLIPSED, Gaia/Hilditch RV convention v=K[cos(nu+w)+e cos w])
        or 270 deg (which=1: secondary eclipsed)."""
        w = np.deg2rad(w_deg)
        nu = (np.pi / 2 if which == 0 else 3 * np.pi / 2) - w
        E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu / 2))
        M = E - e * np.sin(E)
        return Tp + (M / (2 * np.pi) % 1.0) * P

    rng = np.random.default_rng(2)
    N = 20000
    Ps = rng.normal(g.period, g.period_error, N)
    Tps = rng.normal(389.0 + g.t_periastron, g.t_periastron_error, N)
    es = np.clip(rng.normal(g.eccentricity, g.eccentricity_error, N), 0, 0.95)
    ws = rng.normal(g.arg_periastron, g.arg_periastron_error, N)
    for which, name in [(0, "primary eclipsed (nu+w=90)"), (1, "secondary eclipsed (nu+w=270)")]:
        tc = conj_times(Ps, Tps, es, ws, which)
        say(f"\n  Conjunction type: {name}")
        for t_obs, e_obs in zip(t0s, e0s):
            k = np.round((t_obs - np.median(tc)) / g.period)
            pred = tc + k * Ps
            say(f"    K2 eclipse at {t_obs:.3f}+-{e_obs:.3f}: predicted {np.median(pred):.3f} +- {np.std(pred):.3f} (cycle {int(np.median(k))}),"
                f" O-C = {t_obs-np.median(pred):+.3f} d = {(t_obs-np.median(pred))/np.hypot(np.std(pred), e_obs):+.1f} sigma")
        # where would this conjunction fall in TESS S72 (3260.5-3285.6)?
        for k in range(30, 50):
            tt = np.median(tc) + k * g.period
            if 3255 < tt < 3290:
                say(f"    -> this conjunction predicted in TESS S72 window at BTJD {tt:.2f} +- {np.std(tc + k*Ps):.2f}")
    # eclipse duration ratio / separation sanity
    e, w = g.eccentricity, np.deg2rad(g.arg_periastron)
    say(f"\n  Expected duration ratio secondary/primary ~ (1+e sin w)/(1-e sin w) = {(1+e*np.sin(w))/(1-e*np.sin(w)):.2f}"
        f" (convention: primary eclipse at nu+w=90)")
out.close()
