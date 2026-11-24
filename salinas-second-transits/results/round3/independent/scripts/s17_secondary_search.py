"""All-phase secondary-eclipse search, and what it says about the companion.

A self-luminous companion (an M dwarf) passing behind the star gives a secondary eclipse of depth
(transit depth) x S_c/S_* (S = surface brightness in the TESS band). A planet or brown dwarf at
Teq ~ 1100 K gives a few ppm at most. On a circular orbit the eclipse is at phase 0.5 with the
transit's duration; on an eccentric one it moves in phase and duration. This script:
 1. uses TESS-SPOC in every sector it exists and QLP only in sectors without it; bins to 10 min;
    drops the primary transits;
 2. fits a trapezoid eclipse depth at a trial (phase, duration D) jointly with a quadratic baseline
    per orbit over +-2D (no prior detrending, so an eclipse is not partly fitted away); does this at
    every phase (step D/4) for a grid of D; the spread of depth/err over phase sets an empirical
    noise scale per D, by which the errors are inflated;
 3. measures the depth along the (e, omega) locus allowed by the transit density ratio
    g = (rho_circ/rho_star)^(1/3) = (1 + e sin w)/sqrt(1 - e^2), with the occultation phase,
    duration and impact parameter from the orbit, for transit b = 0 and b at its 2-sigma limit (s18);
 4. injects the expected M-dwarf eclipse at a few places to check it is recovered.
Outputs out/secondary_<tic>.txt and .png.
Usage: s17_secondary_search.py TIC"""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import brentq
from common import OUT, load_lc_files, load_tesscut_phot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

tic = int(sys.argv[1])
L = []
pr = L.append
eph = json.load(open(os.path.join(OUT, f"ephem_{tic}.json")))
sh = json.load(open(os.path.join(OUT, f"shape_{tic}.json")))
P, T0 = eph["P"], eph["T0_E0"]
k, a_circ, b_fit = sh["rp"], sh["a"], sh["b"]
DEPTH = sh["depth"]
Dpri = P / np.pi * np.arcsin(np.sqrt((1 + k) ** 2 - b_fit**2) / a_circ)
rtxt = open(os.path.join(OUT, f"rhofixed_{tic}.txt")).read()
a_star = float(re.search(r"-> a/R\* = ([\d.]+)", rtxt).group(1))  # a/R* from the Gaia density
radtxt = open(os.path.join(OUT, f"radius_{tic}.txt")).read()
teff = float(re.search(r"Teff=([\d.]+)", radtxt).group(1))
rstar = float(re.search(r"Gaia-based: L = .*?R = ([\d.]+)", radtxt).group(1))
bp = os.path.join(OUT, f"bprofile_{tic}.txt")
btxt = open(bp).read() if os.path.exists(bp) else ""
m = re.search(r"b < ([\d.]+) at 2 sigma", btxt)
b_up = float(m.group(1)) if m else None
if btxt and not m:  # no 2-sigma limit inside the grid: use the largest b still within dchi2 < 4
    ok = [float(b) for b, d in re.findall(r"b = ([\d.]+):.*scaled dchi2 = +([\d.]+)", btxt) if float(d) < 4]
    b_up = max(ok) if ok else None
bprof = np.array([[float(x) for x in r] for r in
                  re.findall(r"b = ([\d.]+): Rp/R\* = [\d.]+  a/R\* = +([\d.]+)", btxt)]) if btxt else None


def g_of(b):
    """Density ratio (1+e sin w)/sqrt(1-e^2) for transit impact parameter b (a/R* from s18 at that b)."""
    return (np.interp(b, bprof[:, 0], bprof[:, 1]) if bprof is not None else a_circ) / a_star


pr(f"TIC {tic}: P = {P:.6f} d, T14 = {Dpri*24:.2f} h, Rp/R* = {k:.4f}, transit depth {DEPTH*1e6:.0f} ppm, "
   f"Teff = {teff:.0f} K")
pr(f"a/R* from the Gaia density = {a_star:.2f}; largest transit b allowed at 2 sigma (s18): {b_up}")
# main-sequence M dwarf Teff for a radius: coarse dwarf-sequence table, good to ~150 K
MS_R = np.array([0.10, 0.12, 0.14, 0.17, 0.20, 0.27, 0.36, 0.45])
MS_T = np.array([2400, 2650, 2850, 2990, 3060, 3200, 3400, 3550])
Rc = k * rstar
T_ms = float(np.interp(Rc, MS_R, MS_T)) if Rc >= 0.09 else None
pr(f"transiting object: Rp/R* x R* (Gaia {rstar:.2f} Rsun) = {Rc:.3f} Rsun = {Rc*109.1:.1f} Rearth; "
   + (f"a main-sequence star that size would have Teff ~ {T_ms:.0f} K" if T_ms else
      "smaller than any hydrogen-burning star"))

# ---------- data ----------
BIN = 10 / 1440
lcs = load_lc_files(tic)
spoc = {lc["sector"] for lc in lcs if lc["source"] == "TESS-SPOC"}
lcs = [lc for lc in lcs if len(lc["t"]) > 100 and (lc["source"] == "TESS-SPOC" or lc["sector"] not in spoc)]
have = {lc["sector"] for lc in lcs}
lcs += [lc for lc in load_tesscut_phot(tic) if len(lc["t"]) > 100 and lc["sector"] not in have]


def binlc(t, f, e):
    if np.nanmedian(np.diff(t)) >= 0.9 * BIN:
        return t, f, e
    i = np.floor((t - t[0]) / BIN).astype(int)
    w = 1 / e**2
    sw = np.bincount(i, w)
    ok = sw > 0
    return (np.bincount(i, w * t)[ok] / sw[ok], np.bincount(i, w * f)[ok] / sw[ok], 1 / np.sqrt(sw[ok]))


T, F, E, S = [], [], [], []
for lc in lcs:
    t, f, e = binlc(lc["t"], lc["f"], lc["e"])
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    keep = np.abs(ph * P) > 0.75 * Dpri  # drop the primary transits
    T.append(t[keep]); F.append(f[keep]); E.append(e[keep]); S.append(np.full(keep.sum(), lc["sector"]))
T, F, E, S = map(np.concatenate, (T, F, E, S))
pr(f"data: {len(lcs)} sectors ({sum(1 for l in lcs if l['source']=='TESS-SPOC')} TESS-SPOC, "
   f"{sum(1 for l in lcs if l['source']=='QLP')} QLP-only, {sum(1 for l in lcs if l['source']=='TESScut')} "
   f"TESScut-only), {len(T)} points after 10-min binning")


def trap(dt, D):
    """Eclipse shape in [0,1]; ingress/egress length D*k/(1+k)."""
    return np.clip((D / 2 - np.abs(dt)) / (D * k / (1 + k)), 0, 1)


def measure(phi, D, inj=None, sel=None, nmin=3):
    """Depth at (phase, duration) with a quadratic baseline per orbit over +-2D.
    inj = (phase, duration, depth) multiplied into the flux first. Returns (depth, err, n_orbits)."""
    x = (T - T0) / P - phi
    cyc = np.round(x)
    dt = (x - cyc) * P
    w = np.abs(dt) < 2 * D
    if sel is not None:
        w &= sel
    t, dt, f, e, key = T[w], dt[w], F[w].copy(), E[w], (cyc[w] * 1000 + S[w])
    if inj:
        dti = (((t - T0) / P - inj[0] + 0.5) % 1 - 0.5) * P
        f *= 1 - inj[2] * trap(dti, inj[1])
    s = trap(dt, D)
    FR, SR = [], []
    for kk in np.unique(key):
        mm = key == kk
        if (s[mm] > 0.5).sum() < 5 or (dt[mm] < -D / 2).sum() < 5 or (dt[mm] > D / 2).sum() < 5:
            continue
        u = dt[mm] / D
        W = 1 / e[mm]
        X = np.vstack([np.ones_like(u), u, u**2]).T * W[:, None]
        Q, _ = np.linalg.qr(X)
        fw, sw = f[mm] * W, s[mm] * W
        FR.append(fw - Q @ (Q.T @ fw)); SR.append(-(sw - Q @ (Q.T @ sw)))
    if len(FR) < nmin:
        return None
    fr, sr = np.concatenate(FR), np.concatenate(SR)
    den = sr @ sr
    d = (sr @ fr) / den
    chi2 = np.sum((fr - d * sr) ** 2) / max(1, len(fr) - 3 * len(FR) - 1)
    return d, np.sqrt(max(1.0, chi2) / den), len(FR)


# ---------- 1. all-phase search ----------
DG = Dpri * np.array([0.35, 0.5, 0.7, 1.0, 1.4, 2.0])
scans, sig_eff = {}, {}
for D in DG:
    edge = (D + Dpri) / 2 / P + 0.01
    rows = [(phi, *r) for phi in np.arange(edge, 1 - edge, D / 4 / P) if (r := measure(phi, D))]
    a = np.array(rows)
    snr = a[:, 1] / a[:, 2]
    se = 1.4826 * np.median(np.abs(snr - np.median(snr)))
    sig_eff[D] = max(1.0, se)
    a[:, 2] *= sig_eff[D]
    scans[D] = a
    z = a[:, 1] / a[:, 2]
    j = np.argmax(z)
    pr(f"D = {D*24:5.2f} h ({D/Dpri:.2f} T14): {len(a)} phases, noise scale {se:.2f}; "
       f"max depth/err {z[j]:+.2f} at phase {a[j,0]:.4f} ({a[j,1]*1e6:.0f} ppm, {a[j,3]:.0f} orbits); "
       f"min {z.min():+.2f}; median err {np.median(a[:,2])*1e6:.0f} ppm")


def sig_for(D):
    return sig_eff[DG[np.argmin(np.abs(np.log(DG / D)))]]


def meas(phi, D, inj=None, sel=None, nmin=3):
    r = measure(phi, D, inj, sel, nmin)
    return None if r is None else (r[0], r[1] * sig_for(D), r[2])


# ---------- expected depths ----------
lam = np.linspace(600e-9, 1000e-9, 400)


def sratio(Tc):
    bb = lambda T_: 1 / lam**5 / np.expm1(1.4388e-2 / (lam * T_))
    return np.trapezoid(bb(Tc), lam) / np.trapezoid(bb(teff), lam)


pr("expected secondary depth = transit depth x S_c/S_* (TESS band as a 600-1000 nm top hat, blackbodies);")
pr("the ratio to the transit depth is the same whether or not another star dilutes both:")
EXP = {}
for Tc in sorted({2700, 2900, 3100, 3300} | ({int(round(T_ms, -1))} if T_ms else set())):
    EXP[Tc] = DEPTH * sratio(Tc)
    pr(f"  companion Teff {Tc} K: {EXP[Tc]*1e6:.0f} ppm" + ("  <- main-sequence star of this radius"
                                                            if T_ms and Tc == int(round(T_ms, -1)) else ""))
TREF = int(round(T_ms, -1)) if T_ms else 2700
pr(f"  planet or brown dwarf at Teq ~1100 K: thermal {DEPTH*sratio(1100)*1e6:.2f} ppm, reflected (Ag = 0.3) "
   f"{0.3*(k/a_star)**2*1e6:.0f} ppm")


# ---------- 2. eccentric locus ----------
def occ(e, w, b_tr):
    """Occultation phase after transit, duration, impact parameter."""
    def M(f):
        Ea = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(f / 2))
        return Ea - e * np.sin(Ea)
    dphi = ((M(1.5 * np.pi - w) - M(0.5 * np.pi - w)) / (2 * np.pi)) % 1
    r = (1 + e * np.sin(w)) / (1 - e * np.sin(w))
    b_oc = b_tr * r
    ch = (1 + k) ** 2 - b_oc**2
    Docc = Dpri * r * np.sqrt(ch / ((1 + k) ** 2 - b_tr**2)) if ch > 0 else 0.0
    return dphi, Docc, b_oc


e_max = 1 - 2.0 / a_star  # periastron outside 2 stellar radii
pr(f"eccentric orbits with (1 + e sin w)/sqrt(1-e^2) = g: g = {g_of(0):.3f} at transit b = 0"
   + (f", {g_of(b_up):.3f} at b = {b_up:.2f}" if b_up else "") + f"; e < {e_max:.2f} (periastron > 2 R*)")
loc = []
for wd in np.arange(-90, 271, 2.0):
    w = np.radians(wd)
    for b_tr in [0.0] + ([b_up] if b_up else []):
        gb = g_of(b_tr)
        fn = lambda e: (1 + e * np.sin(w)) / np.sqrt(1 - e**2) - gb
        eg = np.linspace(0, 0.995, 996)
        fv = fn(eg)
        roots = [0.0] if abs(fv[0]) < 1e-9 else []
        roots += [brentq(fn, eg[i], eg[i + 1]) for i in np.where(np.sign(fv[:-1]) * np.sign(fv[1:]) < 0)[0]]
        for e in roots:
            if e > e_max:
                continue
            dphi, Docc, b_oc = occ(e, w, b_tr)
            r = meas(dphi, Docc) if Docc > 0 else None
            loc.append((wd, e, b_tr, dphi, Docc, b_oc) + ((r[0], r[1]) if r else (np.nan, np.nan)))
loc = np.array(loc)
for b_tr in sorted(set(loc[:, 2])):
    l = loc[loc[:, 2] == b_tr]
    noecl = l[:, 4] == 0
    good = np.isfinite(l[:, 6])
    pr(f"  transit b = {b_tr:.2f}: omega {l[0,0]:.0f}..{l[-1,0]:.0f} deg, e {l[:,1].min():.2f}..{l[:,1].max():.2f}; "
       f"no occultation for {noecl.sum()}/{len(l)}; not covered {(~noecl & ~good).sum()}; measured {good.sum()}")
    lg = l[good]
    z = lg[:, 6] / lg[:, 7]
    pr(f"    depth/err along the locus {z.min():+.2f}..{z.max():+.2f}; 3-sigma upper limits "
       f"{(lg[:,6]+3*lg[:,7]).min()*1e6:.0f}..{(lg[:,6]+3*lg[:,7]).max()*1e6:.0f} ppm")
    for Tc, ex in EXP.items():
        rej = (ex - lg[:, 6]) / lg[:, 7]
        pr(f"    Teff {Tc} K ({ex*1e6:.0f} ppm): excluded at > 3 sigma for {np.sum(rej > 3)}/{len(lg)} omegas, "
           f"> 2 sigma for {np.sum(rej > 2)}; weakest {rej.min():.1f} sigma at omega {lg[np.argmin(rej),0]:.0f}")
    for x in l[np.linspace(0, len(l) - 1, 9).astype(int)]:
        dtxt = (f"{x[6]*1e6:+5.0f} +- {x[7]*1e6:3.0f} ppm" if np.isfinite(x[6])
                else ("no occultation" if x[4] == 0 else "not covered"))
        pr(f"    omega {x[0]:4.0f}: e = {x[1]:.2f}, phase {x[3]:.3f}, {x[4]*24:5.2f} h, b_occ {x[5]:.2f}: {dtxt}")

# ---------- 3. circular orbit (dilution case) ----------
r5 = meas(0.5, Dpri)
if r5 is None:
    pr(f"circular orbit: phase 0.5 has fewer than 3 orbits with in- and out-of-eclipse data")
else:
    pr(f"circular orbit: phase 0.5, {Dpri*24:.2f} h: {r5[0]*1e6:+.0f} +- {r5[1]*1e6:.0f} ppm over {r5[2]} orbits; "
       f"3-sigma upper limit {(r5[0]+3*r5[1])*1e6:.0f} ppm; excluded: "
       + ", ".join(f"{Tc} K at {(ex-r5[0])/r5[1]:.1f} sigma" for Tc, ex in EXP.items()))

# ---------- 4. injection ----------
dep = EXP[TREF]
l0 = loc[(loc[:, 2] == 0) & np.isfinite(loc[:, 6])]
tests = [(0.5, Dpri)] + [(x[3], x[4]) for x in l0[np.linspace(0, len(l0) - 1, 5).astype(int)[1:4]]]
pr(f"injection of a {dep*1e6:.0f} ppm eclipse into the flux, recovered with the same fit:")
for phi, D in tests:
    r0, r = meas(phi, D), meas(phi, D, inj=(phi, D, dep))
    if r0 is None:
        pr(f"  phase {phi:.3f}, {D*24:5.2f} h: not covered by 3 orbits")
        continue
    pr(f"  phase {phi:.3f}, {D*24:5.2f} h: {r0[0]*1e6:+.0f} -> {r[0]*1e6:+.0f} +- {r[1]*1e6:.0f} ppm "
       f"(recovered {(r[0]-r0[0])/dep:.2f} of injected)")

# ---------- 5. strongest peak: does it hold in independent halves? ----------
best = max(((D, row) for D in DG for row in scans[D]), key=lambda x: x[1][1] / x[1][2])
D, row = best
pr(f"strongest peak of the search: phase {row[0]:.4f}, {D*24:.2f} h, {row[1]*1e6:.0f} +- {row[2]*1e6:.0f} ppm "
   f"({row[1]/row[2]:+.2f}); most negative over all durations "
   f"{min(np.min(scans[d][:,1]/scans[d][:,2]) for d in DG):+.2f}")
tmid = np.median(T)
for name, sel in [(f"before BTJD {tmid:.0f}", T < tmid), (f"after BTJD {tmid:.0f}", T >= tmid),
                  ("sectors < 40 (30-min)", S < 40), ("sectors 40-60", (S >= 40) & (S <= 60)), ("sectors > 60", S > 60)]:
    r = meas(row[0], D, sel=sel)
    if r is None:
        continue
    pr(f"  {name:22s}: {r[0]*1e6:+.0f} +- {r[1]*1e6:.0f} ppm ({r[2]} orbits)")
cy = np.round((T - T0) / P - row[0])
pr("  each orbit alone (white-noise error x the noise scale):")
for c in np.unique(cy):
    r = meas(row[0], D, sel=(cy == c), nmin=1)
    if r is not None:
        pr(f"    E{int(c):4d} (S{int(np.median(S[cy == c]))}, BTJD {T0 + (c + row[0]) * P:.2f}): "
           f"{r[0]*1e6:+5.0f} +- {r[1]*1e6:4.0f} ppm")

# ---------- plot ----------
fig, ax = plt.subplots(3, 1, figsize=(10, 11))
for i, D in enumerate(DG):
    a = scans[D]
    ax[0].plot(a[:, 0], a[:, 1] / a[:, 2] + 6 * i, lw=0.7)
    ax[0].text(1.0, 6 * i, f"{D*24:.1f} h", fontsize=8, va="center")
ax[0].set_xlabel("phase after transit"); ax[0].set_ylabel("depth / err (offset 6 per duration)")
ax[0].set_title(f"TIC {tic}: secondary-eclipse search, quadratic baseline per orbit")
for b_tr, c in zip(sorted(set(loc[:, 2])), ["C0", "C3"]):
    l = loc[(loc[:, 2] == b_tr) & np.isfinite(loc[:, 6])]
    ax[1].errorbar(l[:, 0] + (0.6 if c == "C3" else 0), l[:, 6] * 1e6, l[:, 7] * 1e6, fmt=".", color=c, ms=3,
                   lw=0.5, label=f"transit b = {b_tr:.2f}")
for Tc, ex in EXP.items():
    ax[1].axhline(ex * 1e6, ls="--", lw=0.8, color="k")
    ax[1].text(274, ex * 1e6, f"{Tc} K", fontsize=8, va="center")
ax[1].axhline(0, color="gray", lw=0.5)
ax[1].set_xlim(np.nanmin(loc[:, 0]) - 8, 290)
ax[1].set_xlabel("omega (deg)"); ax[1].set_ylabel("occultation depth on the e-omega locus (ppm)")
ax[1].legend(fontsize=8)
# phase-0.5 fold, quadratic per orbit removed out of eclipse
x = (T - T0) / P - 0.5
cyc = np.round(x)
dt = (x - cyc) * P
w = np.abs(dt) < 2 * Dpri
xs, ys, ws = [], [], []
for kk in np.unique(cyc[w] * 1000 + S[w]):
    mm = w & (cyc * 1000 + S == kk)
    oot = mm & (np.abs(dt) > Dpri / 2)
    if oot.sum() < 10:
        continue
    c = np.polyfit(dt[oot], F[oot], 2, w=1 / E[oot])
    xs.append(dt[mm] * 24); ys.append(F[mm] / np.polyval(c, dt[mm]) - 1); ws.append(1 / E[mm] ** 2)
xs, ys, ws = np.concatenate(xs), np.concatenate(ys), np.concatenate(ws)
bins = np.arange(-2 * Dpri * 24, 2 * Dpri * 24 + 0.01, 0.75)
idx = np.digitize(xs, bins)
B = []
for j in range(1, len(bins)):
    m_ = idx == j
    if m_.sum() > 5:
        mu = np.sum(ws[m_] * ys[m_]) / ws[m_].sum()
        sd = np.sqrt(np.sum(ws[m_] * (ys[m_] - mu) ** 2) / ws[m_].sum())
        B.append((xs[m_].mean(), mu * 1e6, sd / np.sqrt(m_.sum()) * 1e6))
B = np.array(B)
ax[2].errorbar(B[:, 0], B[:, 1], B[:, 2], fmt="o", ms=3,
               label="data, 45-min weighted bins (quadratic per orbit fitted out of eclipse)")
xx = np.linspace(bins[0], bins[-1], 500)
ax[2].plot(xx, -EXP[TREF] * 1e6 * trap(xx / 24, Dpri), "r--",
           label=f"{TREF} K companion, circular orbit ({EXP[TREF]*1e6:.0f} ppm)")
ax[2].plot(xx, -EXP[3300] * 1e6 * trap(xx / 24, Dpri), "m:", label=f"3300 K ({EXP[3300]*1e6:.0f} ppm)")
ax[2].axhline(0, color="gray", lw=0.5)
ax[2].set_xlabel("hours from phase 0.5"); ax[2].set_ylabel("flux - 1 (ppm)"); ax[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, f"secondary_{tic}.png"), dpi=120)
txt = "\n".join(L)
open(os.path.join(OUT, f"secondary_{tic}.txt"), "w").write(txt + "\n")
print(txt)
