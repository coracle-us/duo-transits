"""Inclination history from the physical fits (data/physical_fit_*.json): i(2014.36) from K2, i(2021.75), i(2023.84), i(2026.04)
from TESS. Fit linear and quadratic i(t); find when the primary/secondary eclipses started (i = i_graze) and when i reaches 90 deg;
predict primary/secondary depths for the WASP/KELT epochs (2004-2008) and compare with the measured limits.
Systematic spread taken from the fit variants (l3 = 0, 0.02, 0.15; 2026 secondary phase; 2014 secondary phase 0.32-0.37)."""
import numpy as np, json, glob, sys, pandas as pd
sys.path.insert(0, "scripts")
from ebmodel import flux, ecosw_for_phase
from astropy.time import Time
Pa = 21.7748
def yr(btjd): return Time(btjd + 2457000, format="jd").decimalyear
T = {"2014": yr(-211.557), "2021": yr(2510.37), "2023": yr(3261.6), "2026": yr(4056.38)}
print("epochs (decimal year):", {k: round(v, 3) for k, v in T.items()})
rows = []
for f in sorted(glob.glob("data/physical_fit_*.json")):
    d = json.load(open(f)); p = d["params"]; tag = f.split("physical_fit_")[1][:-5]
    for k2 in d["k2"]:
        if k2["sec_phase_2014"] not in (0.32, 0.3322, 0.34, 0.35): continue
        t = np.array([T["2014"], T["2021"], T["2023"], T["2026"]]); i = np.array([k2["i2014"], p["i2021"], p["i2023"], p["i2026"]])
        lin = np.polyfit(t - 2020, i, 1); quad = np.polyfit(t - 2020, i, 2)
        res = dict(variant=tag, ph14=k2["sec_phase_2014"], i2014=i[0], i2021=i[1], i2023=i[2], i2026=i[3],
                   rate_14_21=(i[1]-i[0])/(t[1]-t[0]), rate_21_23=(i[2]-i[1])/(t[2]-t[1]), rate_23_26=(i[3]-i[2])/(t[3]-t[2]),
                   iP=d["i_graze_P"], iS=d["i_graze_S"])
        for nm, c in [("lin", lin), ("quad", quad)]:
            roots = lambda target: [r.real + 2020 for r in np.roots(np.polyadd(c, [0]*(len(c)-1) + [-target])) if abs(r.imag) < 1e-9]
            def first_before(target, tmax=2014.4):
                rr = [x for x in roots(target) if x < tmax]; return max(rr) if rr else np.nan
            def first_after(target, tmin=2026.0):
                rr = [x for x in roots(target) if x > tmin]; return min(rr) if rr else np.nan
            res[f"{nm}_tP_onset"] = first_before(d["i_graze_P"]); res[f"{nm}_tS_onset"] = first_before(d["i_graze_S"], tmax=2021.7)
            res[f"{nm}_t90"] = first_after(90.0)
            res[f"{nm}_i2006.9"] = np.polyval(c, 2006.9 - 2020)
            res[f"{nm}_resid_max"] = np.max(np.abs(np.polyval(c, t - 2020) - i))
        rows.append(res)
R = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(R.round(3).to_string(index=False))
R.to_csv("data/inclination_trend.csv", index=False)
print("\nSUMMARY (min..max over variants):")
for c in ["i2014", "i2021", "i2023", "i2026", "rate_14_21", "rate_21_23", "rate_23_26", "lin_tP_onset", "quad_tP_onset", "lin_tS_onset",
          "quad_tS_onset", "lin_t90", "quad_t90", "lin_i2006.9", "quad_i2006.9", "iP", "iS"]:
    print(f"  {c:15s} {R[c].min():9.3f} .. {R[c].max():9.3f}   (base, ph14=0.32: {R[(R.variant=='base')&(R.ph14==0.32)][c].values[0]:.3f})")
# depth vs inclination curve (base model, 2021 e,w), TESS band and "V-ish" band (LD 0.45,0.30; J scaled to 0.85)
d = json.load(open("data/physical_fit_base.json")); p = d["params"]
ec = ecosw_for_phase(0.33221, p["esinw"]); e = np.hypot(ec, p["esinw"]); w = np.arctan2(p["esinw"], ec)
tt = np.linspace(-0.6, 0.6, 1601); curve = []
for inc in np.arange(84.0, 90.001, 0.05):
    fp = flux(tt, 0, "P", Pa, e, w, np.radians(inc), p["r1"], p["r2"], p["J"]); fs = flux(tt, 0, "S", Pa, e, w, np.radians(inc), p["r1"], p["r2"], p["J"])
    fpv = flux(tt, 0, "P", Pa, e, w, np.radians(inc), p["r1"], p["r2"], 0.85, u=(0.45, 0.30)); fsv = flux(tt, 0, "S", Pa, e, w, np.radians(inc), p["r1"], p["r2"], 0.85, u=(0.45, 0.30))
    curve.append((inc, 1 - fp.min(), 1 - fs.min(), 1 - fpv.min(), 1 - fsv.min()))
C = pd.DataFrame(curve, columns=["inc", "dP_tess", "dS_tess", "dP_V", "dS_V"]); C.to_csv("data/depth_vs_inclination.csv", index=False)
print("\nmax depths at i=90 (TESS): P %.3f  S %.3f" % (C.dP_tess.iloc[-1], C.dS_tess.iloc[-1]))
b = R[(R.variant == "base") & (R.ph14 == 0.32)].iloc[0]
for nm in ["lin", "quad"]:
    i06 = b[f"{nm}_i2006.9"]
    print(f"{nm}: predicted i(2006.9)={i06:.2f} -> V-band primary depth {np.interp(i06, C.inc, C.dP_V):.4f}, secondary {np.interp(i06, C.inc, C.dS_V):.4f}")
