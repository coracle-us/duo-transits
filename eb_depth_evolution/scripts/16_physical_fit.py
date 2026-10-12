"""Joint physical fit of all TESS primary+secondary eclipses (own TESScut r=1.5 photometry):
global r1=R1/a, r2=R2/a, J=S2/S1 (TESS band), e*sin(w); inclination per season (2021, 2023, 2026); free time of conjunction per
event (absorbs the ~10-25 min ETVs); e*cos(w) per season fixed from the measured mean secondary phase. Quadratic LD (0.20, 0.27),
third light 0.02 (neighbour G=14.9 at 21"). Then: K2 2014 primary fitted for i_2014 with geometry fixed; predicted 2014 secondary.
usage: python 16_physical_fit.py [secondary phase assumed for 2026] [l3]"""
import numpy as np, pandas as pd, sys, json, glob
from scipy.optimize import least_squares
from astropy.io import fits
sys.path.insert(0, "scripts")
from ebmodel import flux, ecosw_for_phase
from common import fit_event
Pa = 21.7748
PH26 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.32817
L3 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.02
TAG = sys.argv[3] if len(sys.argv) > 3 else "base"
ev = pd.read_csv("data/eclipse_times.csv")
phase_season = {"2021": 0.33221, "2023": 0.32817, "2026": PH26}


def season(c):
    return "2021" if c <= 5 else ("2023" if c <= 40 else "2026")


events = []
for _, e in ev[ev.src.str.startswith("TESS")].iterrows():
    d = pd.read_csv(f"data/phot/tess_s{int(e.sector)}.csv"); d = d[d.q == 0]
    t = d.time.values; f = (d["f_r1.5"] / np.nanmedian(d["f_r1.5"])).values
    win = 0.55 if e.kind == "P" else 0.7
    r = fit_event(t, f, np.full(len(f), 1e-3), e.tc, win=win)
    m = np.abs(t - e.tc) < win
    tt, ff = t[m], f[m] / np.interp(t[m], r["t"], r["base"])
    texp = 600 / 86400 if e.sector < 50 else 200 / 86400
    sig = np.std(ff[np.abs(tt - e.tc) > 0.3]) if np.sum(np.abs(tt - e.tc) > 0.3) > 10 else 2e-3
    events.append(dict(kind=e.kind, cycle=int(e.cycle), season=season(e.cycle), t=tt, f=ff, tc0=e.tc, texp=texp, sig=sig))
nE = len(events)
names = ["r1", "r2", "J", "esinw", "i2021", "i2023", "i2026"] + [f"dt_{e['kind']}{e['cycle']}" for e in events]


def unpack(x):
    return dict(zip(names, x))


def ew(p, seas, ecw=None):
    es = p["esinw"]
    ec = ecw if ecw is not None else ecosw_for_phase(phase_season[seas], es)
    return np.hypot(ec, es), np.arctan2(es, ec)


def model_event(p, ev_, t=None, u=(0.2, 0.27), l3=L3, inc=None, ecw=None):
    e, w = ew(p, ev_["season"], ecw)
    i = np.radians(inc if inc is not None else p["i" + ev_["season"]])
    key = f"dt_{ev_['kind']}{ev_['cycle']}"
    return flux(ev_["t"] if t is None else t, ev_["tc0"] + p.get(key, 0.0), ev_["kind"], Pa, e, w, i,
                p["r1"], p["r2"], p["J"], u=u, l3=l3, texp=ev_["texp"], nsub=5)


def resid(x):
    p = unpack(x)
    return np.concatenate([(ev_["f"] - model_event(p, ev_)) / ev_["sig"] for ev_ in events])


best = None
starts = [(0.037, 0.037, 1.0, 0.07, 86.0), (0.045, 0.03, 0.9, -0.07, 86.0), (0.03, 0.045, 1.1, 0.07, 86.5),
          (0.04, 0.04, 1.0, -0.05, 85.5), (0.05, 0.025, 1.0, 0.1, 86.5), (0.025, 0.05, 1.0, -0.1, 86.5)]
for r1, r2, J, es, i0 in starts:
    x0 = [r1, r2, J, es, i0, i0 + 0.2, i0 + 0.4] + [0.0] * nE
    lb = [0.01, 0.01, 0.3, -0.5, 80, 80, 80] + [-0.05] * nE
    ub = [0.12, 0.12, 3.0, 0.5, 90, 90, 90] + [0.05] * nE
    try:
        r = least_squares(resid, x0, bounds=(lb, ub), x_scale="jac", max_nfev=3000)
    except Exception as ex:
        print("fail", ex); continue
    chi2 = np.sum(r.fun ** 2)
    print(f"start r1={r1} r2={r2} esinw={es}: chi2={chi2:.1f} ->", np.round(r.x[:7], 5))
    if best is None or chi2 < best[1]:
        best = (r, chi2)
r, chi2 = best
p = unpack(r.x); N = len(r.fun)
cov = np.linalg.pinv(r.jac.T @ r.jac) * chi2 / (N - len(r.x)); err = np.sqrt(np.diag(cov))
print(f"\nBEST chi2={chi2:.1f} N={N} red.chi2={chi2 / (N - len(r.x)):.3f}  (2026 sec phase assumed {PH26}, l3={L3})")
for n, v, s in zip(names, r.x, err):
    print(f"  {n:10s} {v: .6f} +- {s:.6f}" + (f"   ({v * 1440:+.1f} min)" if n.startswith("dt_") else ""))
print(f"  corr(i2021,i2023)={cov[4,5]/err[4]/err[5]:.3f}; sigma(i2023-i2021)={np.sqrt(cov[4,4]+cov[5,5]-2*cov[4,5]):.4f}; "
      f"sigma(i2026-i2023)={np.sqrt(cov[5,5]+cov[6,6]-2*cov[5,6]):.4f}")
for s_ in ["2021", "2023", "2026"]:
    e, w = ew(p, s_)
    print(f"  season {s_}: e={e:.5f} w={np.degrees(w):.2f} deg  ecosw={e*np.cos(w):.5f}")
out = {"params": dict(zip(names, map(float, r.x))), "errors": dict(zip(names, map(float, err))), "chi2": float(chi2), "N": int(N),
       "phase_season": phase_season, "l3": L3}
out["cov_i"] = cov[4:7, 4:7].tolist()
# ---- K2 2014 primary: fit i_2014 with Kepler-band LD (0.30, 0.30), l3=0, J_Kp = J (twin A stars)
evd = fits.open(glob.glob("data/k2/**/hlsp_everest*.fits", recursive=True)[0])[1].data
tk = evd["TIME"] + 2454833 - 2457000; fk = evd["FLUX"]
ok = np.isfinite(tk) & np.isfinite(fk) & (tk > -232); tk, fk = tk[ok].astype(float), fk[ok].astype(float) / np.nanmedian(fk[ok])
k2 = pd.read_csv("data/k2_primary_fits.csv").set_index("product").loc["EVEREST"]
rk = fit_event(tk, fk, np.full(len(fk), 4e-4), k2.tc, win=0.6)
m = np.abs(tk - k2.tc) < 0.6; tkk = tk[m]; fkk = fk[m] / np.interp(tkk, rk["t"], rk["base"])
evk = dict(kind="P", cycle=-124, season="2014", t=tkk, f=fkk, tc0=k2.tc, texp=1765.5 / 86400)
res_k2 = []
tt = np.linspace(-0.6, 0.6, 1201)
for ph14 in [0.32, 0.3322, 0.34, 0.35, 0.36, 0.37]:
    ec14 = ecosw_for_phase(ph14, p["esinw"])

    def rk2(x):
        q = dict(p); q["dt_P-124"] = x[1]
        return (fkk - model_event(q, evk, inc=x[0], ecw=ec14, u=(0.30, 0.30), l3=0.0)) / 4e-4
    rr = least_squares(rk2, [p["i2021"] - 0.5, 0.0], bounds=([80, -0.05], [90, 0.05]))
    i14 = rr.x[0]
    e14, w14 = np.hypot(ec14, p["esinw"]), np.arctan2(p["esinw"], ec14)
    fs = flux(tt, 0.0, "S", Pa, e14, w14, np.radians(i14), p["r1"], p["r2"], p["J"], u=(0.3, 0.3))
    fp = flux(tt, 0.0, "P", Pa, e14, w14, np.radians(i14), p["r1"], p["r2"], p["J"], u=(0.3, 0.3))
    res_k2.append(dict(sec_phase_2014=ph14, i2014=float(i14), k2_chi2=float(np.sum(rr.fun ** 2)), n=int(len(fkk)),
                       model_P_depth=float(1 - fp.min()), model_S_depth=float(1 - fs.min())))
print(pd.DataFrame(res_k2).to_string(index=False))
out["k2"] = res_k2
rows = []
for s_ in ["2021", "2023", "2026"]:
    inc = p["i" + s_]; e, w = ew(p, s_)
    fp = flux(tt, 0, "P", Pa, e, w, np.radians(inc), p["r1"], p["r2"], p["J"], l3=0.0)
    fs = flux(tt, 0, "S", Pa, e, w, np.radians(inc), p["r1"], p["r2"], p["J"], l3=0.0)
    rP = (1 - e ** 2) / (1 + e * np.sin(w)); rS = (1 - e ** 2) / (1 - e * np.sin(w))
    rows.append(dict(season=s_, inc=inc, depthP_undiluted=1 - fp.min(), depthS_undiluted=1 - fs.min(),
                     bP_over_rsum=rP * np.cos(np.radians(inc)) / (p["r1"] + p["r2"]), bS_over_rsum=rS * np.cos(np.radians(inc)) / (p["r1"] + p["r2"])))
print(pd.DataFrame(rows).to_string(index=False))
out["seasons"] = rows
# grazing limits: inclination at which primary / secondary eclipses vanish (b = 1) for 2021 e,w
e, w = ew(p, "2021")
for kind, rr_ in [("P", (1 - e ** 2) / (1 + e * np.sin(w))), ("S", (1 - e ** 2) / (1 - e * np.sin(w)))]:
    out[f"i_graze_{kind}"] = float(np.degrees(np.arccos((p["r1"] + p["r2"]) / rr_)))
    print(f"  eclipse {kind} vanishes for i < {out[f'i_graze_{kind}']:.3f} deg (sky separation at conjunction = r1+r2)")
json.dump(out, open(f"data/physical_fit_{TAG}.json", "w"), indent=1, default=float)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, axs = plt.subplots(3, 4, figsize=(16, 9)); axs = axs.ravel()
best14 = min(res_k2, key=lambda z: z["k2_chi2"]); ec14 = ecosw_for_phase(best14["sec_phase_2014"], p["esinw"])
for ax, ev_ in zip(axs, events + [evk]):
    if ev_["season"] == "2014":
        mod = model_event(dict(p), ev_, inc=best14["i2014"], ecw=ec14, u=(0.3, 0.3), l3=0.0)
    else:
        mod = model_event(p, ev_)
    ax.plot(ev_["t"] - ev_["tc0"], ev_["f"], "k.", ms=1.5); ax.plot(ev_["t"] - ev_["tc0"], mod, "r-")
    ax.set_title(f"{ev_['kind']} cyc {ev_['cycle']} ({ev_['season']})", fontsize=9); ax.set_ylim(0.6, 1.02)
plt.tight_layout(); plt.savefig(f"figures/physical_fit_events_{TAG}.png", dpi=80)
