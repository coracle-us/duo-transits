"""Step 12: per-eclipse fits (TESS) with the global shape fixed.

For every eclipse: free inclination i_k, time shift dt_k and a linear baseline; shape (esw, r1, r2,
J, e cos w of that season, LD, l3) fixed at the nominal global solution.  Errors: covariance with
per-eclipse noise (rms * beta from s10).  Outputs per-eclipse mid-time (model minimum), depth,
inclination -> per_eclipse.csv, plus ephemeris / O-C -> timing.csv and fig_timing.png
usage: s12_per_eclipse.py [tag]
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ebdata, ebfit, ebmodel as M

H = os.path.dirname(os.path.abspath(__file__))
tag = sys.argv[1] if len(sys.argv) > 1 else "nominal"
pj = json.load(open(os.path.join(H, "fit_tess_%s.json" % tag)))
pars = {k: v for k, v in pj.items() if not k.startswith("_")}
chp = os.path.join(H, "data", "chain_%s.npz" % tag)
ch = np.load(chp if os.path.exists(chp) else os.path.join(H, "data", "chain_nominal.npz"))
errs = json.loads(str(ch["errs"]))
df = ebdata.load(pj["_src"], col=pj["_col"], seasons=[2021, 2023, 2026])
df = df[df.ecl != "S72"].reset_index(drop=True)
for ecl, (s1, beta) in errs.items():
    df.loc[df.ecl == ecl, "err"] = s1 * beta
F = ebfit.Fitter(df, pars, [], **pj["_fitter_kw"])
P, T0 = pars["P"], pars["T0"]

rows = []
for g in F.groups:
    s = g["season"]
    def res(x):
        p = dict(F.pars); p["dt_" + g["ecl"]] = x[1]
        m = F.model_group(p, g, inc=x[0])
        c = F.baseline(m, g)
        return (g["y"] - m * (c[0] + c[1] * (g["t"] - g["tc"]))) / g["e"]
    r = least_squares(res, [F.pars["i%d" % s], F.pars["dt_" + g["ecl"]]], bounds=([80, -0.2], [90, 0.2]),
                      x_scale=[0.05, 0.002])
    cov = np.linalg.pinv(r.jac.T @ r.jac)
    chi2 = 2 * r.cost
    # model minimum time & depth (as observed, our aperture) and undiluted depth
    p = dict(F.pars); p["dt_" + g["ecl"]] = r.x[1]
    tt = np.linspace(g["tc"] - 0.5, g["tc"] + 0.5, 20001)
    gg = dict(g); gg["t"] = tt; gg["nsub"] = 1
    m = F.model_group(p, gg, inc=r.x[0])
    tmin = tt[np.argmin(m)]
    m0 = F.model_group(p, gg, inc=r.x[0], l3=0.0)
    # depth uncertainty: propagate i error numerically
    di = np.sqrt(cov[0, 0])
    mp = F.model_group(p, gg, inc=min(r.x[0] + di, 90))
    mm = F.model_group(p, gg, inc=r.x[0] - di)
    d_err = 0.5 * abs(mm.min() - mp.min())
    # coverage: fraction of model eclipse duration covered by data
    inecl = m < 1 - 1e-4
    t_in = tt[inecl]
    cov_frac = np.mean([np.any(np.abs(g["t"] - x) < 0.01) for x in t_in[::20]])
    rows.append(dict(ecl=g["ecl"], type=g["typ"], epoch=int(g["ecl"][1:]), season=s, i=r.x[0], i_err=di,
                     dt=r.x[1], dt_err=np.sqrt(cov[1, 1]), tmid=tmin, tmid_err=np.sqrt(cov[1, 1]),
                     depth_obs=1 - m.min(), depth_obs_err=d_err, depth_undiluted=1 - m0.min(),
                     duration=t_in.max() - t_in.min() if len(t_in) else np.nan, coverage=cov_frac,
                     chi2=chi2, n=len(g["t"]), rms_ppt=errs[g["ecl"]][0] * 1e3, beta=errs[g["ecl"]][1]))
E = pd.DataFrame(rows).sort_values("tmid").reset_index(drop=True)

# linear ephemeris to TESS primaries
def linfit(n, t, s):
    A = np.vstack([np.ones_like(n), n]).T / s[:, None]
    c, *_ = np.linalg.lstsq(A, t / s, rcond=None)
    cov = np.linalg.inv(A.T @ A)
    return c, cov
pr = E[E.type == "P"]
c, cov = linfit(pr.epoch.values.astype(float), pr.tmid.values, pr.tmid_err.values)
E["oc_P_ephem_min"] = (E.tmid - (c[0] + c[1] * E.epoch)) * 1440.0
chi2 = np.sum(((pr.tmid - (c[0] + c[1] * pr.epoch)) / pr.tmid_err) ** 2)
print("Primary ephemeris: T0 = %.5f +- %.5f BTJD, P = %.6f +- %.6f d ; chi2 = %.1f for %d dof"
      % (c[0], np.sqrt(cov[0, 0]), c[1], np.sqrt(cov[1, 1]), chi2, len(pr) - 2))
sc = E[E.type == "S"]
c2, cov2 = linfit(sc.epoch.values.astype(float), sc.tmid.values, sc.tmid_err.values)
chi2s = np.sum(((sc.tmid - (c2[0] + c2[1] * sc.epoch)) / sc.tmid_err) ** 2)
print("Secondary ephemeris: T0s = %.5f +- %.5f, P = %.6f +- %.6f ; chi2 = %.1f for %d dof"
      % (c2[0], np.sqrt(cov2[0, 0]), c2[1], np.sqrt(cov2[1, 1]), chi2s, len(sc) - 2))
E["oc_Sown_ephem_min"] = np.where(E.type == "S", (E.tmid - (c2[0] + c2[1] * E.epoch)) * 1440.0, np.nan)
pd.set_option("display.width", 250)
print(E[["ecl", "season", "tmid", "tmid_err", "oc_P_ephem_min", "i", "i_err", "depth_obs", "depth_obs_err",
         "depth_undiluted", "duration", "coverage", "chi2", "n"]].to_string(float_format=lambda x: "%.5f" % x))
# secondary minus primary intervals, same cycle
iv = []
for n in sorted(set(E.epoch)):
    a = E[(E.epoch == n) & (E.type == "P")]; b = E[(E.epoch == n) & (E.type == "S")]
    if len(a) and len(b):
        iv.append(dict(epoch=n, season=int(a.season.iloc[0]), interval=b.tmid.iloc[0] - a.tmid.iloc[0],
                       err=np.hypot(a.tmid_err.iloc[0], b.tmid_err.iloc[0]),
                       note="S partial" if b.coverage.iloc[0] < 0.9 else ""))
IV = pd.DataFrame(iv)
print(IV.to_string(float_format=lambda x: "%.5f" % x))
E.to_csv(os.path.join(H, "per_eclipse_%s.csv" % tag), index=False)
IV.to_csv(os.path.join(H, "intervals_%s.csv" % tag), index=False)
json.dump(dict(T0P=c[0], T0P_err=np.sqrt(cov[0, 0]), P=c[1], P_err=np.sqrt(cov[1, 1]), chi2=chi2, dof=len(pr) - 2,
               T0S=c2[0], PS=c2[1], PS_err=np.sqrt(cov2[1, 1]), chi2S=chi2s, dofS=len(sc) - 2),
          open(os.path.join(H, "ephemeris_%s.json" % tag), "w"), indent=1)

fig, ax = plt.subplots(1, 3, figsize=(16, 4))
for typ, col in (("P", "b"), ("S", "r")):
    d = E[E.type == typ]
    ax[0].errorbar(d.tmid, (d.tmid - (c[0] + c[1] * d.epoch)) * 1440, d.tmid_err * 1440, fmt="o", color=col,
                   label=typ)
ax[0].set_ylabel("O-C vs primary linear ephemeris (min)"); ax[0].legend(); ax[0].set_xlabel("BTJD")
for s in (2021, 2023, 2026):
    d = E[E.season == s]
    for typ, col in (("P", "b"), ("S", "r")):
        dd = d[d.type == typ]
        ax[1].errorbar(dd.epoch, (dd.tmid - (c[0] + c[1] * dd.epoch)) * 1440, dd.tmid_err * 1440, fmt="o", color=col)
ax[1].set_xlabel("epoch"); ax[1].set_title("O-C (zoom by epoch)")
ax[2].errorbar(E.tmid, E.i, E.i_err, fmt="o")
ax[2].set_xlabel("BTJD"); ax[2].set_ylabel("i per eclipse (deg)")
plt.tight_layout(); plt.savefig(os.path.join(H, "fig_timing_%s.png" % tag), dpi=80)
