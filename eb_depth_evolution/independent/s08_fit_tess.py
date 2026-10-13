"""Step 8: global least-squares EB fit to all TESS eclipses (own photometry).

Stage A (linear ephemeris, done first) showed derivative-shaped residuals -> eclipse timing
variations.  Stage B: every eclipse gets its own time shift; P and T0 then only define the
reference grid.  e cos w per season is updated iteratively from the mean secondary-minus-primary
interval of that season (ecw2026 tied to ecw2023: no 2026 secondary was observed).
Usage: s08_fit_tess.py [tag] [l3] [u1a u1b u2a u2b] [aperture-col]
"""
import os, sys, json, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import ebdata, ebfit, ebmodel as M

H = os.path.dirname(os.path.abspath(__file__))
tag = sys.argv[1] if len(sys.argv) > 1 else "nominal"
l3 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.035
u = [float(x) for x in sys.argv[3:7]] if len(sys.argv) > 6 else [0.20, 0.30, 0.30, 0.25]
col = sys.argv[7] if len(sys.argv) > 7 else "flux"
src = sys.argv[8] if len(sys.argv) > 8 else "own"

df = ebdata.load(src, col=col, seasons=[2021, 2023, 2026])
df = df[df.ecl != "S72"]
p0 = json.load(open(os.path.join(H, "fit_tess_ls.json")))
pars = {k: v for k, v in p0.items() if not k.startswith("_")}
pars.setdefault("ecw2026", pars["ecw2023"])
pars["ecw2026"] = pars["ecw2023"]
pars.setdefault("i2014", 85.0); pars.setdefault("ecw2014", pars["ecw2021"])
present = sorted(df.season.unique())
shape = ["esw", "r1", "r2", "J"] + ["i%d" % s for s in present]
if os.environ.get("ESW_PER_SEASON"):
    pars["esw2021"] = pars["esw"]; pars["esw2023"] = pars["esw"]
    shape = ["esw2021", "esw2023", "r1", "r2", "J"] + ["i%d" % s for s in present]
    pars["esw2026"] = pars["esw"]
F = ebfit.Fitter(df, pars, [], u1=u[:2], u2=u[2:], l3=l3)
dts = ["dt_" + g["ecl"] for g in F.groups]
F = ebfit.Fitter(df, F.pars, shape + dts, u1=u[:2], u2=u[2:], l3=l3)

t0 = time.time(); F.resid(np.array([F.pars[n] for n in F.free])); print("eval time %.3f s" % (time.time() - t0))
for it in range(5):
    r = F.fit(max_nfev=3000, xtol=1e-10, ftol=1e-10)
    p = F.pars
    chi2 = np.sum(r.fun ** 2)
    msg = []
    for s in [s for s in (2021, 2023) if s in present]:
        dP = [p["dt_" + g["ecl"]] for g in F.groups if g["season"] == s and g["typ"] == "P"]
        dS = [p["dt_" + g["ecl"]] for g in F.groups if g["season"] == s and g["typ"] == "S"]
        delta = np.mean(dS) - np.mean(dP)
        sp = F.season_pars(p, s)
        off = M.secondary_offset(p["P"], sp["e"], sp["w"])
        new = ebfit.ecw_from_interval(p["P"], p.get("esw%d" % s, p["esw"]), off + delta)
        msg.append("%d: delta=%.5f ecw %.5f->%.5f" % (s, delta, p["ecw%d" % s], new))
        F.pars["ecw%d" % s] = new
        for g in F.groups:
            if g["season"] == s and g["typ"] == "S":
                F.pars["dt_" + g["ecl"]] -= delta
    F.pars["ecw2026"] = F.pars["ecw2023"]
    if "esw2023" in F.pars:
        F.pars["esw2026"] = F.pars["esw2023"]
    print("iter %d chi2=%.1f N=%d | %s" % (it, chi2, len(r.fun), " ; ".join(msg)))
r = F.fit(max_nfev=3000, xtol=1e-12, ftol=1e-12)
chi2 = np.sum(r.fun ** 2)
print("FINAL chi2=%.1f N=%d red=%.2f" % (chi2, len(r.fun), chi2 / (len(r.fun) - len(F.free))))
for n in [k for k in ("esw2021", "esw2023") if k in F.pars] + ["P", "T0", "esw", "ecw2021", "ecw2023", "r1", "r2", "J", "i2021", "i2023", "i2026"]:
    print("  %-8s %.6f" % (n, F.pars[n]))
for g in F.groups:
    print("  dt_%-4s %+.5f d" % (g["ecl"], F.pars["dt_" + g["ecl"]]))
out = dict(F.pars)
out["_fitter_kw"] = dict(u1=u[:2], u2=u[2:], l3=l3)
out["_chi2"] = chi2; out["_N"] = len(r.fun); out["_src"] = src; out["_col"] = col
json.dump(out, open(os.path.join(H, "fit_tess_%s.json" % tag), "w"), indent=1)
