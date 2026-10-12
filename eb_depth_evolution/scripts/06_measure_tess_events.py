"""Per-event primary/secondary depth & time from own TESScut photometry (r=1.0,1.5,2.0,2.5 px) and SPOC/QLP products."""
import glob, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import lightkurve as lk, warnings; warnings.filterwarnings("ignore")
from common import *
rows = []; plots = []
def events_in(tmin, tmax):
    out = []
    for k in range(-400, 200):
        for kind, ph in [("P", 0.0), ("S", 0.3300)]:
            tt = T00 + (k+ph)*P0
            if tmin - 0.3 < tt < tmax + 0.3: out.append((kind, k, tt))
    return out
def run(t, f, e, label, sector):
    for kind, k, tt in events_in(np.nanmin(t), np.nanmax(t)):
        # refine guess for secondary: search the deepest point within +-0.4 d
        m = np.abs(t - tt) < 0.5
        if m.sum() < 20: continue
        tg = t[m][np.argmin(pd.Series(f[m]).rolling(5, center=True, min_periods=1).median().values)]
        r = fit_event(t, f, e, tg, win=0.8 if kind=="P" else 0.9)
        if r is None: continue
        # completeness: fraction of +-1.2*w window that has data
        tw = np.abs(t - r["tc"]) < 1.5*r["w"]
        comp = tw.sum()*np.median(np.diff(t))/(3*r["w"])
        rows.append(dict(sector=sector, product=label, kind=kind, cycle=k, tc=r["tc"], e_tc=r["e_tc"], depth=r["D"], e_depth=r["e_D"],
                         w=r["w"], p=r["p"], beta=r["beta"], rms=r["rms"], n=r["n"], coverage=round(min(comp, 1.0), 2), covers_min=r["covers_min"]))
        plots.append((label, kind, k, r))
for fn in sorted(glob.glob("data/phot/tess_s*.csv")):
    sec = int(fn.split("_s")[1][:-4]); d = pd.read_csv(fn); d = d[d.q == 0]
    for ap in ["1.0", "1.5", "2.0", "2.5"]:
        y = d[f"f_r{ap}"].values; y = y/np.nanmedian(y); t = d.time.values
        ok = np.isfinite(y); e = np.full(ok.sum(), np.std(np.diff(y[ok]))/np.sqrt(2))
        run(t[ok], y[ok], e, f"own r={ap}", sec)
# pipeline products
for fn in glob.glob("data/tess_lc/**/*.fits", recursive=True):
    lc = lk.read(fn)
    auth = "QLP" if "qlp" in fn.lower() else "TESS-SPOC"
    col = "sap_flux" if auth == "QLP" else "pdcsap_flux"
    try: lc = lk.read(fn, flux_column=col)
    except Exception: pass
    lc = lc.remove_nans(); t = np.asarray(lc.time.value, float); y = np.asarray((lc.flux/np.nanmedian(lc.flux)).value, float)
    e = np.full(len(y), np.std(np.diff(y))/np.sqrt(2))
    run(t, y, e, f"{auth} {col}", lc.sector)
df = pd.DataFrame(rows).sort_values(["kind", "tc", "product"])
df.to_csv("data/tess_events_all_products.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print(df.to_string(index=False, float_format=lambda x: f"{x:.5f}"))
# plot own r=1.5
sel = [p for p in plots if p[0] == "own r=1.5"]
fig, axs = plt.subplots(int(np.ceil(len(sel)/4)), 4, figsize=(16, 3*np.ceil(len(sel)/4)))
for ax, (lab, kind, k, r) in zip(axs.flat, sel):
    ax.plot(r["t"]-r["tc"], r["f"]/r["base"], "k.", ms=2); ax.plot(r["t"]-r["tc"], r["model"]/r["base"], "r-")
    ax.set_title(f"{kind} cyc {k} D={r['D']:.4f}"); 
plt.tight_layout(); plt.savefig("figures/tess_event_fits_r1.5.png", dpi=80)
