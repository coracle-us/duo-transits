"""Step 14: assemble final numbers -> results.csv, timing.csv, summary.json, fig_summary.png.

Reads the nominal fit/posterior and every systematic variant.  All numbers in REPORT.txt come
from the printed output of this script (s14_summary.out.txt).
"""
import os, json, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ebmodel as M

H = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(H, "data")
yr = lambda btjd: 2000.0 + (np.asarray(btjd) + 2457000.0 - 2451545.0) / 365.25
nom = json.load(open(os.path.join(H, "fit_tess_nominal.json")))
sd = json.load(open(os.path.join(H, "season_depths_nominal.json")))
E = pd.read_csv(os.path.join(H, "per_eclipse_nominal.csv"))
eph = json.load(open(os.path.join(H, "ephemeris_nominal.json")))
k2 = json.load(open(os.path.join(H, "k2_results_nominal.json")))
P = nom["P"]
q = lambda key: (sd[key][2], sd[key][3] - sd[key][2], sd[key][2] - sd[key][1])


def dmin(inc, e, w, which):
    fc = (np.pi / 2 if which == "P" else 3 * np.pi / 2) - w
    f = fc + np.linspace(-0.3, 0.3, 4001)
    r = (1 - e ** 2) / (1 + e * np.cos(f))
    s = np.sin(w + f)
    return np.min(r * np.sqrt(1 - np.sin(np.radians(inc)) ** 2 * s ** 2))


def thresholds(p, s=2023):
    esw = p.get("esw%d" % s, p["esw"]); ecw = p["ecw%d" % s]
    e, w = np.hypot(ecw, esw), np.arctan2(esw, ecw)
    g = lambda inc, which, target: dmin(inc, e, w, which) - target
    out = {}
    for which in ("P", "S"):
        out["ion_" + which] = brentq(g, 60, 89.9999, args=(which, p["r1"] + p["r2"]))
        dd = abs(p["r2"] - p["r1"])
        out["itot_" + which] = brentq(g, 60, 90, args=(which, dd)) if g(90, which, dd) < 0 else 90.0
    return out


def depths(p, s, l3, u=None):
    u = u or p["_fitter_kw"]
    esw = p.get("esw%d" % s, p["esw"]); ecw = p["ecw%d" % s]
    e, w = np.hypot(ecw, esw), np.arctan2(esw, ecw)
    out = {}
    for which in ("P", "S"):
        tc = 0.0 if which == "P" else M.secondary_offset(P, e, w)
        t = tc + np.linspace(-0.6, 0.6, 6001)
        m = M.lc(t, P, 0.0, e, w, p["i%d" % s], p["r1"], p["r2"], p["J"], u["u1"][0], u["u1"][1],
                 u["u2"][0], u["u2"][1], l3=l3)
        out[which] = 1 - m.min()
    return out


lines = []
def pr(*a):
    s = " ".join(str(x) for x in a)
    print(s); lines.append(s)

# ---------------- variants ----------------
VAR = {"l3=0": "l3_000", "l3=0.15": "l3_015", "LD u1-0.1": "ld_u1m", "LD u1+0.1": "ld_u1p", "LD u2-0.1": "ld_u2m",
       "LD u2+0.1": "ld_u2p", "LD star1!=star2": "ld_diff", "aperture r=1.5": "ap15", "aperture r=2.5": "ap25",
       "aperture r=3.0": "ap30", "QLP reduction": "qlp2", "esw per season": "esw_season",
       "l3=-0.025 (data-preferred)": "l3scan_-0.025"}
V = {}
for lab, tag in VAR.items():
    pj = json.load(open(os.path.join(H, "fit_tess_%s.json" % tag)))
    k2v = os.path.join(H, "k2_results_%s.json" % tag)
    i14 = json.load(open(k2v))["primary_by_pipeline"]["k2sff"]["i2014"] if os.path.exists(k2v) else np.nan
    V[lab] = dict(pj=pj, i2014=i14)
V["nominal"] = dict(pj=nom, i2014=k2["primary_by_pipeline"]["k2sff"]["i2014"])

pr("=" * 100)
pr("A. FIT QUALITY / NOMINAL GLOBAL SOLUTION (own TESScut r=2.0 px photometry, l3=0.035, LD (0.20,0.30) both stars)")
pr("   chi2 = %.1f for N = %d points (per-eclipse errors = OOT rms, before beta inflation)" % (nom["_chi2"], nom["_N"]))
for k in ["r1", "r2", "J", "esw", "i2021", "i2023", "i2026"]:
    pr("   %-6s = %.5f +%.5f -%.5f (stat)" % ((k,) + q(k)))
pr("   e cos w: 2021 %.4f, 2023 %.4f (2026 tied to 2023);  e(2023) = %.4f, w(2023) = %.2f deg" %
   (nom["ecw2021"], nom["ecw2023"], q("e")[0], q("w_deg")[0]))
pe = V["esw per season"]["pj"]
w21 = np.degrees(np.arctan2(pe["esw2021"], pe["ecw2021"])); w23 = np.degrees(np.arctan2(pe["esw2023"], pe["ecw2023"]))
pr("   esw-per-season variant: esw2021=%.4f esw2023=%.4f -> w2021=%.2f w2023=%.2f deg (dw=%.2f deg in %.2f yr); "
   "e=%.4f/%.4f; chi2 %.1f vs %.1f nominal"
   % (pe["esw2021"], pe["esw2023"], w21, w23, w23 - w21, yr(3258) - yr(2510), np.hypot(pe["esw2021"], pe["ecw2021"]),
      np.hypot(pe["esw2023"], pe["ecw2023"]), pe["_chi2"], nom["_chi2"]))

# ---------------- inclinations ----------------
pr("=" * 100)
pr("B. INCLINATION PER SEASON: statistical (MCMC 16-84%) and systematic (variant fits, LS) ")
seasons = [2014, 2021, 2023, 2026]
tref = {2014: k2["primary_by_pipeline"]["k2sff"]["tmin"]}
for s in (2021, 2023, 2026):
    tref[s] = E[E.season == s].tmid.mean()
iv = {}
hdr = "   %-28s" % "variant" + "".join("%10s" % ("i%d" % s) for s in seasons)
pr(hdr)
for lab in ["nominal"] + list(VAR):
    pj = V[lab]["pj"]
    vals = [V[lab]["i2014"]] + [pj.get("i%d" % s, np.nan) if not (lab == "QLP reduction" and s == 2026) else np.nan
                                for s in (2021, 2023, 2026)]
    iv[lab] = dict(zip(seasons, vals))
    pr("   %-28s" % lab + "".join("%10.3f" % v for v in vals))
k2s = k2["i2014_syst"]
groups = {"third light (0..0.15)": ["l3=0", "l3=0.15"], "limb darkening": [k for k in VAR if k.startswith("LD")],
          "aperture": [k for k in VAR if k.startswith("aperture")], "reduction (QLP)": ["QLP reduction"],
          "e sin w per season": ["esw per season"]}
syst = {s: {} for s in seasons}
for gname, labs in groups.items():
    for s in seasons:
        d = [iv[l][s] - iv["nominal"][s] for l in labs if np.isfinite(iv[l][s])]
        syst[s][gname] = max(np.abs(d)) if d else 0.0
syst[2014]["K2 LD (Kepler band)"] = max(abs(k2s[k] - iv["nominal"][2014]) for k in ("UK-0.1", "UK+0.1", "u2-0.1", "u2+0.1"))
syst[2014]["2014 e,w (ecw +-0.03; apsidal extrap.)"] = max(abs(k2s[k] - iv["nominal"][2014])
                                                           for k in ("ecw-0.03", "ecw+0.03", "apsidal-extrap"))
istat = {2014: ((k2["i2014_post_k2sff"]["p84"] - k2["i2014_post_k2sff"]["p16"]) / 2)}
for s in (2021, 2023, 2026):
    istat[s] = (q("i%d" % s)[1] + q("i%d" % s)[2]) / 2
itab = []
pr("")
pr("   season  year      i_nom    stat     " + "  ".join("%s" % g for g in list(groups) + ["K2-only terms"]))
for s in seasons:
    tot = np.sqrt(sum(v ** 2 for v in syst[s].values()))
    k2only = np.sqrt(sum(v ** 2 for k, v in syst[s].items() if k.startswith(("K2", "2014"))))
    pr("   %d  %.3f  %8.3f  %6.3f   " % (s, yr(tref[s]), iv["nominal"][s], istat[s]) +
       "  ".join("%6.3f" % syst[s].get(g, 0) for g in groups) + "   %6.3f  -> syst total %.3f" % (k2only, tot))
    itab.append(dict(season=s, year=yr(tref[s]), i=iv["nominal"][s], i_stat=istat[s], i_syst=tot, i_syst_k2only=k2only))
IT = pd.DataFrame(itab)

# ---------------- di/dt ----------------
pr("=" * 100)
pr("C. RATE OF INCLINATION CHANGE")
def wfit(x, y, s, deg=1):
    A = np.vander(x, deg + 1)
    W = 1 / s ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    c = cov @ A.T @ (W * y)
    chi2 = np.sum(((y - A @ c) / s) ** 2)
    return c, cov, chi2
x = IT.year.values - 2020.0; y = IT.i.values
# errors for the shape/constancy test: stat + non-common systematics (K2 band/epoch terms for 2014;
# aperture+reduction terms for TESS treated as season-specific is conservative)
s_rel = np.sqrt(IT.i_stat.values ** 2 + IT.i_syst_k2only.values ** 2)
s_rel = np.maximum(s_rel, 0.005)
rates = {}
for a, b in [(2014, 2021), (2021, 2023), (2023, 2026), (2014, 2026), (2021, 2026)]:
    ia, ib = IT.set_index("season").loc[a], IT.set_index("season").loc[b]
    r = (ib.i - ia.i) / (ib.year - ia.year)
    er = np.hypot(np.sqrt(ia.i_stat ** 2 + ia.i_syst_k2only ** 2), np.sqrt(ib.i_stat ** 2 + ib.i_syst_k2only ** 2)) / (ib.year - ia.year)
    # systematic range of the rate over the TESS variants (same variant applied to both seasons)
    rv = []
    for lab in VAR:
        va, vb = iv[lab][a], iv[lab][b]
        if np.isfinite(va) and np.isfinite(vb):
            rv.append((vb - va) / (ib.year - ia.year))
    rates["%d-%d" % (a, b)] = (r, er, min(rv), max(rv))
    pr("   %d -> %d : di/dt = %.4f +- %.4f deg/yr (stat+epoch-specific)   variant range %.4f .. %.4f" % (a, b, r, er, min(rv), max(rv)))
cT, covT, chiT = wfit(x[1:], y[1:], s_rel[1:], 1)
pr("   TESS-only linear (2021,2023,2026): di/dt = %.4f +- %.4f deg/yr, chi2 = %.2f (1 dof)" % (cT[0], np.sqrt(covT[0, 0]), chiT))
cL, covL, chiL = wfit(x, y, s_rel, 1)
pr("   All-epoch linear (2014..2026): di/dt = %.4f +- %.4f deg/yr, chi2 = %.1f (2 dof)  -> %s" %
   (cL[0], np.sqrt(covL[0, 0]), chiL, "NOT constant" if chiL > 13.8 else "consistent with constant"))
cQ, covQ, chiQ = wfit(x, y, s_rel, 2)
pr("   All-epoch quadratic: i = %.4f + %.4f (t-2020) + %.5f (t-2020)^2 ; d2i/dt2 = %.4f +- %.4f deg/yr^2, chi2 = %.2f (1 dof)"
   % (cQ[2], cQ[1], cQ[0], 2 * cQ[0], 2 * np.sqrt(covQ[0, 0]), chiQ))
for t_ in (2014.36, 2021.8, 2023.85, 2026.04):
    pr("      quadratic rate at %.2f: %.3f deg/yr" % (t_, cQ[1] + 2 * cQ[0] * (t_ - 2020)))
# within-season slopes from per-eclipse inclinations
for s in (2021, 2023):
    d = E[E.season == s]
    c, cov, chi = wfit(yr(d.tmid.values) - 2020, d.i.values, d.i_err.values, 1)
    dp = d[d.type == "P"]
    cp, covp, chip = wfit(yr(dp.tmid.values) - 2020, dp.i.values, dp.i_err.values, 1) if len(dp) > 2 else (None, None, None)
    pr("   within %d (per-eclipse i, P+S, %d eclipses): slope %.3f +- %.3f deg/yr, chi2/dof = %.1f/%d; per-eclipse i rms %.4f deg"
       % (s, len(d), c[0], np.sqrt(cov[0, 0]), chi, len(d) - 2, np.std(d.i)))

# ---------------- thresholds and onset/totality ----------------
pr("=" * 100)
pr("D. GEOMETRIC THRESHOLDS (2023 e,w) AND ONSET / TOTALITY EPOCHS")
th = {}
for lab in ["nominal"] + list(VAR):
    pj = V[lab]["pj"]
    th[lab] = thresholds(pj)
pr("   %-28s %8s %8s %8s %8s   %s" % ("variant", "ion_P", "ion_S", "itot_P", "itot_S", "r1, r2"))
for lab in th:
    pj = V[lab]["pj"]
    pr("   %-28s %8.3f %8.3f %8.3f %8.3f   %.4f %.4f" % (lab, th[lab]["ion_P"], th[lab]["ion_S"], th[lab]["itot_P"],
                                                       th[lab]["itot_S"], pj["r1"], pj["r2"]))
pr("   MCMC (nominal): ion_P %.3f +%.3f -%.3f ; ion_S %.3f ; itot_P %.3f +%.3f -%.3f" %
   (q("ion_P") + (q("ion_S")[0],) + q("itot_P")))

def t_cross(itarget, model):
    if model == "lin_14_21":
        r = rates["2014-2021"][0]; return yr(tref[2014]) + (itarget - iv["nominal"][2014]) / r
    if model == "lin_TESS":
        return 2020 + (itarget - cT[1]) / cT[0]
    if model == "quad":
        a, b, c = cQ
        roots = np.roots([a, b, c - itarget])
        roots = roots[np.isreal(roots)].real + 2020
        return roots
ion = th["nominal"]["ion_P"]
pr("   Primary-eclipse onset (i = ion_P = %.3f):" % ion)
pr("     linear, 2014-2021 rate (%.3f deg/yr) back from 2014:  %.1f" % (rates["2014-2021"][0], t_cross(ion, "lin_14_21")))
pr("     linear, TESS-era rate (%.3f deg/yr) back from 2014 K2 point: %.1f" %
   (cT[0], yr(tref[2014]) - (iv["nominal"][2014] - ion) / cT[0]))
qr = t_cross(ion, "quad"); pr("     quadratic all epochs: roots %s" % np.round(qr, 1))
# onset range over variants (2014-2021 rate)
on_v = []
for lab in th:
    i14 = V[lab]["i2014"]; i21 = iv[lab][2021]
    if np.isfinite(i14) and np.isfinite(i21):
        r = (i21 - i14) / (yr(tref[2021]) - yr(tref[2014]))
        on_v.append(yr(tref[2014]) - (i14 - th[lab]["ion_P"]) / r)
pr("     variant range (2014-2021 linear): %.1f .. %.1f" % (min(on_v), max(on_v)))
ionS = th["nominal"]["ion_S"]
pr("   Secondary-eclipse onset (i = ion_S = %.3f): 2014-2021 linear -> %.1f ; quadratic -> %s" %
   (ionS, t_cross(ionS, "lin_14_21"), np.round(t_cross(ionS, "quad"), 1)))
itp = th["nominal"]["itot_P"]
pr("   Primary totality (star 1 fully covered by the larger star 2; i >= itot_P = %.3f):" % itp)
pr("     TESS-era linear: %.1f ; quadratic: %s" % (t_cross(itp, "lin_TESS"), np.round(t_cross(itp, "quad"), 1)))
tv = []
for lab in th:
    if lab == "QLP reduction":
        continue
    pj = V[lab]["pj"]
    i26 = pj["i2026"];
    tv.append(yr(tref[2026]) + (th[lab]["itot_P"] - i26) / cT[0])
pr("     variant range (TESS-era linear from each variant's i2026): %.1f .. %.1f" % (min(tv), max(tv)))
pr("   i = 90 deg (deepest possible primary): TESS-era linear %.1f ; quadratic %s" %
   (t_cross(90.0, "lin_TESS"), np.round(t_cross(90.0, "quad"), 1)))
pr("   Secondary totality i >= itot_S = %.3f: TESS-era linear %.1f" % (th["nominal"]["itot_S"], t_cross(th["nominal"]["itot_S"], "lin_TESS")))

# ---------------- depths ----------------
pr("=" * 100)
pr("E. ECLIPSE DEPTHS PER SEASON")
pr("   Definition: D = 1 - F_min/F_out of the best-fit physical model at mid-eclipse (instantaneous).")
pr("   'obs' = in our r=2.0 px TESScut aperture (nominal l3=0.035);  'undil' = corrected to l3=0.")
dep_rows = []
for s in (2021, 2023, 2026):
    for which in ("P", "S"):
        key = "D%s%d" % (which, s)
        o, u_ = q(key + "_obs"), q(key + "_undil")
        d = E[(E.season == s) & (E.type == which)]
        sc = np.std(d.depth_obs) if len(d) > 1 else np.nan
        # aperture/reduction systematic on the as-observed depth
        dv = []
        for lab in ("aperture r=1.5", "aperture r=2.5", "aperture r=3.0", "QLP reduction"):
            pj = V[lab]["pj"]
            if lab == "QLP reduction" and s == 2026:
                continue
            dv.append(depths(pj, s, pj["_fitter_kw"]["l3"])[which])
        dsys = max(abs(np.array(dv) - o[0]))
        observed = len(d) > 0
        dep_rows.append(dict(season=s, type=which, D_obs=o[0], D_obs_stat=(o[1] + o[2]) / 2, D_ecl_scatter=sc,
                             D_obs_syst=dsys, D_undil=u_[0], observed=observed, n_ecl=len(d)))
        pr("   %d %s  D_obs = %.4f +- %.4f (stat) [eclipse-to-eclipse rms %s; aperture/QLP range +-%.4f]  D_undil = %.4f  %s"
           % (s, which, o[0], (o[1] + o[2]) / 2, ("%.4f" % sc) if np.isfinite(sc) else "n/a", dsys, u_[0],
              "(%d eclipse%s)" % (len(d), "s" if len(d) != 1 else "") if observed else "(NOT OBSERVED: model prediction)"))
# K2
pk = k2["primary_by_pipeline"]["k2sff"]
esw0, ecw14 = nom["esw"], k2["ecw2014_assumed"]
e14, w14 = np.hypot(ecw14, esw0), np.arctan2(esw0, ecw14)
Jk = k2["J_Kp_nominal"]
t = np.linspace(-0.5, 0.5, 8001)
mK = M.lc(t, P, 0, e14, w14, iv["nominal"][2014], nom["r1"], nom["r2"], Jk, 0.28, 0.30, 0.28, 0.30, l3=0.0)
mT = M.lc(t, P, 0, e14, w14, iv["nominal"][2014], nom["r1"], nom["r2"], nom["J"], 0.20, 0.30, 0.20, 0.30, l3=0.035)
pr("   2014 P  Kepler band (K2SFF/PDCSAP, l3~0): D_model = %.4f (instantaneous); lowest 30-min point %.4f (K2SFF), %.4f (PDCSAP)"
   % (1 - mK.min(), pk["depth_obs_min"], k2["primary_by_pipeline"]["k2pdc"]["depth_obs_min"]))
pr("           same geometry in TESS band & our aperture: D_obs-equivalent = %.4f" % (1 - mT.min()))
dep_rows.append(dict(season=2014, type="P", D_obs=1 - mK.min(), D_obs_stat=0.0015, D_undil=1 - mK.min(), observed=True,
                     n_ecl=1, note="Kepler band; TESS-band equivalent %.4f" % (1 - mT.min())))

# K2 secondary: see s15_k2_secondary_check.py (full light curves, random-epoch calibration)
ksc = json.load(open(os.path.join(H, "k2_secondary_check.json")))
for w_ in ksc["windows"]:
    pr("   K2 secondary window %-6s centre BTJD %.3f (box scan +-0.5 d, widths 0.15-0.35 d): best dip %.5f, S/N %.1f"
       % (w_["src"], w_["centre"], w_["dip"], w_["snr"]))
pr("   random off-eclipse epochs (K2SFF, N=%d): best-box S/N median %.1f, 99%% %.1f; best-box depth median %.5f, 99%% %.5f"
   % (ksc["n_rand"], ksc["rand_snr_p50"], ksc["rand_snr_p99"], ksc["rand_dep_p50"], ksc["rand_dep_p99"]))
sec_lim = ksc["rand_dep_p99"]
pr("   => no 2014 secondary; depth < %.4f (99%% false-alarm level of the box scan, Kepler band); the model predicts"
   " no secondary for i2014=%.3f < ion_S=%.3f" % (sec_lim, iv["nominal"][2014], ionS))
dep_rows.append(dict(season=2014, type="S", D_obs=np.nan, D_upper3sig=sec_lim, observed=False, n_ecl=0,
                     note="not detected in 2 K2 windows; Kepler band"))
DEP = pd.DataFrame(dep_rows)

# ---------------- timing ----------------
pr("=" * 100)
pr("F. ECLIPSE TIMES")
T0p, Pp = eph["T0P"], eph["P"]
E["cycle_phase_epoch"] = E.epoch
E["calc_lin"] = T0p + Pp * E.epoch
mean_int = E[E.type == "S"].tmid.values  # placeholder
IV = pd.read_csv(os.path.join(H, "intervals_nominal.csv"))
mI = IV.interval.mean()
E["OC_min"] = np.where(E.type == "P", (E.tmid - E.calc_lin) * 1440, (E.tmid - E.calc_lin - mI) * 1440)
# realistic ephemeris errors: rescale by sqrt(chi2/dof)
sc = np.sqrt(eph["chi2"] / eph["dof"])
pr("   Linear ephemeris (TESS primaries, weighted): T0 = %.5f BTJD, P = %.6f d" % (T0p, Pp))
pr("     formal errors %.5f / %.7f d; chi2 = %.0f for %d dof -> errors x%.0f: %.4f d / %.6f d" %
   (eph["T0P_err"], eph["P_err"], eph["chi2"], eph["dof"], sc, eph["T0P_err"] * sc, eph["P_err"] * sc))
# unweighted fit (more honest given ETVs)
pp = E[E.type == "P"]
cu = np.polyfit(pp.epoch, pp.tmid, 1); ru = pp.tmid - np.polyval(cu, pp.epoch)
covu = np.cov(np.polyfit(pp.epoch, pp.tmid, 1, cov=True)[1]) if False else np.polyfit(pp.epoch, pp.tmid, 1, cov=True)[1]
pr("   Unweighted linear fit: P = %.6f +- %.6f d, T0 = %.5f +- %.5f ; rms O-C = %.1f min" %
   (cu[0], np.sqrt(covu[0, 0]), cu[1], np.sqrt(covu[1, 1]), np.std(ru, ddof=2) * 1440))
pr("   %-5s %-5s %14s %9s %10s %9s %8s %8s" % ("ecl", "seas", "Tmid(BTJD)", "err(min)", "O-C(min)", "cover", "i", "D_obs"))
for _, r in E.sort_values("tmid").iterrows():
    pr("   %-5s %-5d %14.5f %9.2f %10.2f %9.2f %8.3f %8.4f" % (r.ecl, r.season, r.tmid, r.tmid_err * 1440, r.OC_min,
                                                            r.coverage, r.i, r.depth_obs))
pr("   (secondary O-C relative to primary ephemeris + mean S-P interval %.5f d)" % mI)
k2t = pk["tmin"]; k2e = pk["dt_err"]
oc_k2 = (k2t - (T0p + Pp * (-124))) * 1440
oc_k2u = (k2t - np.polyval(cu, -124)) * 1440
ext_err = np.sqrt(covu[1, 1] + 124 ** 2 * covu[0, 0] + 2 * (-124) * covu[0, 1]) * 1440
pr("   K2 primary (epoch -124): Tmid = %.5f BTJD (K2SFF; PDCSAP %.5f); O-C = %.1f min (weighted eph), %.1f +- %.1f min (unweighted eph, extrapolated)"
   % (k2t, k2["primary_by_pipeline"]["k2pdc"]["tmin"], oc_k2, oc_k2u, ext_err))
pr("   Secondary-minus-primary intervals (same cycle):")
for _, r in IV.iterrows():
    pr("     epoch %3d (%d): %.5f +- %.5f d %s" % (r.epoch, r.season, r.interval, r.err, r.note if isinstance(r.note, str) else ""))
pr("     within-2021 change: %.1f min/cycle; within-2023: %.1f min/cycle; 2021->2023 mean change: %.2f min/cycle" %
   ((IV.interval[1] - IV.interval[0]) * 1440, (IV.interval[3] - IV.interval[2]) * 1440,
    (IV.interval[2:].mean() - IV.interval[:2].mean()) / 35 * 1440))
for s in (2021, 2023):
    d = E[(E.season == s)]
    for typ in ("P", "S"):
        dd = d[d.type == typ].sort_values("epoch")
        if len(dd) >= 2:
            lp = np.diff(dd.tmid) / np.diff(dd.epoch)
            pr("     local %s period %d: %s d" % (typ, s, ", ".join("%.5f" % v for v in lp)))

# ---------------- output tables ----------------
res = []
for s in seasons:
    row = dict(season=s, year=round(yr(tref[s]), 3), i_deg=IT.set_index("season").loc[s, "i"],
               i_stat=IT.set_index("season").loc[s, "i_stat"], i_syst=IT.set_index("season").loc[s, "i_syst"])
    for which in ("P", "S"):
        d = DEP[(DEP.season == s) & (DEP.type == which)]
        if len(d):
            d = d.iloc[0]
            row["D%s_obs" % which] = d.D_obs
            row["D%s_obs_err" % which] = np.sqrt(np.nan_to_num(d.get("D_obs_stat", 0)) ** 2 + np.nan_to_num(d.get("D_obs_syst", 0)) ** 2 +
                                                 np.nan_to_num(d.get("D_ecl_scatter", 0)) ** 2)
            row["D%s_undil" % which] = d.get("D_undil", np.nan)
            row["D%s_observed" % which] = bool(d.observed)
            if "D_upper3sig" in d and np.isfinite(d.get("D_upper3sig", np.nan)):
                row["D%s_upper3sig" % which] = d.D_upper3sig
    ee = E[E.season == s]
    row["T_primaries_BTJD"] = ";".join("%.5f" % v for v in ee[ee.type == "P"].tmid) if s != 2014 else "%.5f" % k2t
    row["T_secondaries_BTJD"] = ";".join("%.5f" % v for v in ee[ee.type == "S"].tmid)
    row["OC_primaries_min"] = ";".join("%.1f" % v for v in ee[ee.type == "P"].OC_min) if s != 2014 else "%.1f" % oc_k2
    row["band"] = "Kepler" if s == 2014 else "TESS"
    res.append(row)
RES = pd.DataFrame(res)
RES.to_csv(os.path.join(H, "results.csv"), index=False)
E.to_csv(os.path.join(H, "timing.csv"), index=False)
json.dump(dict(rates={k: list(v) for k, v in rates.items()}, tess_linear=[cT[0], np.sqrt(covT[0, 0]), chiT],
               all_linear=[cL[0], np.sqrt(covL[0, 0]), chiL], quad=cQ.tolist(), thresholds=th),
          open(os.path.join(H, "summary.json"), "w"), indent=1, default=float)
open(os.path.join(H, "s14_summary.out.txt"), "w").write("\n".join(lines) + "\n")

# ---------------- figure ----------------
fig, ax = plt.subplots(2, 2, figsize=(13, 9))
a = ax[0, 0]
for which, col in (("P", "C0"), ("S", "C3")):
    d = DEP[(DEP.type == which) & (DEP.observed)]
    a.errorbar([yr(tref[s]) for s in d.season], d.D_obs, fmt="o", color=col, label="%s observed" % which)
    d2 = DEP[(DEP.type == which) & (~DEP.observed) & DEP.D_obs.notna()]
    a.plot([yr(tref[s]) for s in d2.season], d2.D_obs, "o", mfc="none", color=col, label="%s predicted (not observed)" % which)
a.plot(yr(tref[2014]), sec_lim, "v", color="C3", label="2014 S 3-sigma limit")
a.set_ylabel("depth (fraction)"); a.legend(fontsize=8); a.set_title("Eclipse depths (2014 = Kepler band)")
a = ax[0, 1]
a.errorbar(IT.year, IT.i, np.sqrt(IT.i_stat ** 2 + IT.i_syst ** 2), fmt="ko", capsize=3, label="season i (stat+syst)")
tt = np.linspace(2006, 2031, 300)
a.plot(tt, np.polyval(cQ, tt - 2020), "b-", lw=0.8, label="quadratic")
a.plot(tt, np.polyval(cT, tt - 2020), "g--", lw=0.8, label="TESS-only linear")
a.axhline(ion, color="r", ls=":", label="primary onset %.2f" % ion)
a.axhline(ionS, color="m", ls=":", label="secondary onset %.2f" % ionS)
a.axhline(itp, color="k", ls=":", label="primary totality %.2f" % itp)
a.set_ylim(83.5, 90.2); a.set_xlim(2006, 2031); a.legend(fontsize=7); a.set_ylabel("i (deg)")
a = ax[1, 0]
for typ, col in (("P", "C0"), ("S", "C3")):
    d = E[E.type == typ]
    a.errorbar(d.epoch, d.OC_min, d.tmid_err * 1440, fmt="o", color=col, label=typ)
a.errorbar(-124, oc_k2u, ext_err, fmt="s", color="C0", mfc="none", label="K2 P (vs unweighted eph.)")
a.set_xlabel("epoch"); a.set_ylabel("O-C (min)"); a.legend(fontsize=8); a.set_title("O-C, primary ephemeris")
a = ax[1, 1]
a.errorbar(IV.epoch, IV.interval, IV.err, fmt="ko")
a.set_xlabel("epoch"); a.set_ylabel("T_sec - T_pri (d)"); a.set_title("secondary-minus-primary interval")
plt.tight_layout(); plt.savefig(os.path.join(H, "fig_summary.png"), dpi=80)
