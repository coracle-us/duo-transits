"""Plot data + model + residuals for every eclipse of a saved fit.  usage: s09_plotfit.py fit.json [out.png] [season_dt]"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ebdata, ebfit

H = os.path.dirname(os.path.abspath(__file__))
pj = json.load(open(sys.argv[1]))
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(H, "fig_fit.png")
pars = {k: v for k, v in pj.items() if not k.startswith("_")}
df = ebdata.load(pj.get("_src", "own"), col=pj.get("_col", "flux"), seasons=[2021, 2023, 2026])
df = df[df.ecl != "S72"]
F = ebfit.Fitter(df, pars, [], **pj.get("_fitter_kw", {}))
th = F.pars
n = len(F.groups)
fig, axes = plt.subplots(2, n, figsize=(2.4 * n, 5), gridspec_kw=dict(height_ratios=[3, 1]), sharex="col")
for k, g in enumerate(F.groups):
    m = F.model_group(th, g)
    c = F.baseline(m, g)
    b = c[0] + c[1] * (g["t"] - g["tc"])
    x = g["t"] - g["tc"]
    axes[0, k].plot(x, g["y"] / b, "k.", ms=1.5)
    tt = np.linspace(g["t"].min(), g["t"].max(), 1500)
    gg = dict(g); gg["t"] = tt
    axes[0, k].plot(tt - g["tc"], F.model_group(th, gg), "r-", lw=1)
    res = g["y"] / b - m
    axes[1, k].plot(x, res * 1e3, "k.", ms=1.5)
    axes[1, k].axhline(0, color="r")
    axes[0, k].set_title("%s (%d)\nrms=%.2f ppt" % (g["ecl"], g["season"], np.std(res) * 1e3), fontsize=8)
    axes[1, k].set_ylim(-15, 15)
plt.tight_layout()
plt.savefig(out, dpi=70)
print("chi2", np.sum(F.resid_p(th) ** 2))
