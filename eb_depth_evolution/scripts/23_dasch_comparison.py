"""DASCH faint-tail test: is the fraction of plate points >0.3 / >0.4 mag below the median larger for EPIC 202126877 than for
comparison stars of similar magnitude/colour within 40 arcmin (same plates, same quality cuts)? An eclipsing configuration like
2023-2026 would put ~1% of random-epoch points >0.3 mag faint. Also checks whether the target's faint points cluster in time."""
import json, io, os, numpy as np, pandas as pd, requests
from scipy.stats import binomtest
os.makedirs("data/dasch/comps", exist_ok=True)
C = pd.read_csv("data/dasch/comparison_candidates.csv")
C = C[(C.color > -0.1) & (C.color < 0.9) & (C.stdmag > 10.85) & (C.stdmag < 11.75) & (C.ref_number != 406315521252856)]


def get(ref, gsc):
    fn = f"data/dasch/comps/{ref}.json"
    if not os.path.exists(fn):
        r = requests.post("https://api.starglass.cfa.harvard.edu/public/dasch/dr7/lightcurve",
                          json={"refcat": "apass", "gsc_bin_index": int(gsc), "ref_number": int(ref)}, timeout=300)
        open(fn, "w").write(r.text)
    return pd.read_csv(io.StringIO("\n".join(json.load(open(fn)))))


def clean(d):
    d = d[d.magcal_magdep.notna()]
    q = d[(d.magcal_local_rms < 0.3) & (d.limiting_mag_local > d.magcal_magdep + 0.75) & (d.aflags == 0) & (d.reject_flag == 0)
          & (d.plate_center_dist_deg < 15) & (d.blended_mag.isna())].copy()
    q["dm"] = q.magcal_magdep - np.median(q.magcal_magdep); return q


tgt = clean(pd.read_csv(io.StringIO("\n".join(json.load(open("data/dasch/lc_raw.json"))))))
rows = []
for _, c in C.iterrows():
    try:
        q = clean(get(c.ref_number, c.gsc_bin_index))
    except Exception as e:
        print("fail", c.ref_text, e); continue
    if len(q) < 300: continue
    s = 1.48 * np.median(np.abs(q.dm))
    rows.append(dict(star=c.ref_text, mag=c.stdmag, color=c.color, n=len(q), rsig=s, f03=(q.dm > 0.3).mean(), f04=(q.dm > 0.4).mean(),
                     n03=int((q.dm > 0.3).sum()), fbright03=(q.dm < -0.3).mean()))
R = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print(R.round(4).to_string(index=False))
s_t = 1.48 * np.median(np.abs(tgt.dm))
print(f"\nTARGET: n={len(tgt)} rsig={s_t:.3f} f(>0.3)={(tgt.dm>0.3).mean():.4f} ({(tgt.dm>0.3).sum()}) f(>0.4)={(tgt.dm>0.4).mean():.4f} "
      f"f(< -0.3)={(tgt.dm<-0.3).mean():.4f}")
pool03 = R.n03.sum() / R.n.sum(); pool04 = (R.f04 * R.n).sum() / R.n.sum()
print(f"comparison pool: N_stars={len(R)} N_pts={R.n.sum()} f(>0.3)={pool03:.4f} f(>0.4)={pool04:.4f}; median rsig {R.rsig.median():.3f}")
k = int((tgt.dm > 0.3).sum())
print(f"binomial P(>= {k} faint | pooled rate) = {binomtest(k, len(tgt), pool03, alternative='greater').pvalue:.3g}")
k4 = int((tgt.dm > 0.4).sum())
print(f"binomial P(>= {k4} faint>0.4 | pooled rate) = {binomtest(k4, len(tgt), pool04, alternative='greater').pvalue:.3g}")
# normalised by each star's own scatter: count points > 3.2 rsig
for nm, q, s in [("target", tgt, s_t)]:
    print(f"target f(>3.2 rsig)={(q.dm > 3.2*s).mean():.4f}")
print("comparison f(>3.2 rsig) median:", np.median([ (clean(get(r_.ref_number, r_.gsc_bin_index)).dm > 3.2*rs).mean()
                                                     for (_, r_), rs in zip(C[C.ref_text.isin(R.star)].iterrows(), R.rsig)]))
R.to_csv("data/dasch/comparison_faint_tail.csv", index=False)
