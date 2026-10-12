"""Compile the depth history: model-free depth (median flux within +-15 min of minimum vs local out-of-eclipse median) for every
event, own TESScut r=1.5 (TESS) and EVEREST/K2SFF (K2); generalized-Gaussian depths from all products for the spread; S/P ratios."""
import numpy as np, pandas as pd, glob
from astropy.io import fits
from astropy.time import Time
ev = pd.read_csv("data/eclipse_times.csv"); allp = pd.read_csv("data/tess_events_all_products.csv")
rows = []
def mf(t, f, tc, half_out=0.35):
    out = (np.abs(t - tc) > half_out) & (np.abs(t - tc) < 0.85)
    inn = np.abs(t - tc) < 0.0105
    base = np.median(f[out]); mn = np.median(f[inn]); e = np.std(f[out])/np.sqrt(max(inn.sum(), 1))
    return 1 - mn/base, e/base, int(inn.sum())
for _, e in ev.iterrows():
    if e.src.startswith("TESS"):
        d = pd.read_csv(f"data/phot/tess_s{int(e.sector)}.csv"); d = d[d.q == 0]
        t = d.time.values; f = (d["f_r1.5"]/np.nanmedian(d["f_r1.5"])).values
        D, eD, n = mf(t, f, e.tc, 0.35 if e.kind == "P" else 0.45)
        g = allp[(allp.kind == e.kind) & (allp.cycle == e.cycle) & (allp.covers_min)]
        g = g[~g["product"].str.contains("r=2.5")]
        spread = (g.depth.min(), g.depth.max())
        prods = ", ".join(f"{p_}:{x:.4f}" for p_, x in zip(g["product"], g.depth))
    else:
        ef = fits.open(glob.glob("data/k2/**/hlsp_everest*.fits", recursive=True)[0])[1].data
        t = (ef["TIME"] + 2454833 - 2457000).astype(float); f = ef["FLUX"].astype(float); ok = np.isfinite(t) & np.isfinite(f); t, f = t[ok], f[ok]/np.median(f[ok])
        D, eD, n = mf(t, f, e.tc)
        k2 = pd.read_csv("data/k2_primary_fits.csv"); spread = (k2.depth[:2].min(), k2.depth[:2].max())
        prods = ", ".join(f"{p_}:{x:.4f}" for p_, x in zip(k2["product"], k2.depth))
    rows.append(dict(year=round(Time(e.tc + 2457000, format="jd").decimalyear, 3), btjd=round(e.tc, 5), kind=e.kind, cycle=e.cycle,
                     source="TESS S%d own r=1.5" % e.sector if e.src.startswith("TESS") else "K2 C0 EVEREST",
                     depth_modelfree=round(D, 4), err=round(eD, 4), n_min=n, gg_depth_range=f"{spread[0]:.4f}-{spread[1]:.4f}",
                     gg_width_d=round(e.w, 4), coverage=e.coverage, products=prods))
T = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print(T.drop(columns=["products"]).to_string(index=False))
T.to_csv("results_depth_events.csv", index=False)
# season means and dilution-independent S/P ratio
T["season"] = np.where(T.year < 2015, "2014 K2", np.where(T.year < 2022, "2021 S43-45", np.where(T.year < 2024, "2023 S71-72", "2026 S1751")))
S = T[T.coverage > 0.8].groupby(["season", "kind"]).depth_modelfree.agg(["mean", "std", "size"]).unstack("kind")
print(S.round(4))
for s in S.index:
    try:
        rP = S.loc[s, ("mean", "P")]; rS = S.loc[s, ("mean", "S")]
        print(f"{s}: S/P = {rS/rP:.3f}")
    except Exception: pass
