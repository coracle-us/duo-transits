"""Step 15: closer look at the K2 secondary windows (are the ~0.05% box 'dips' real?).

Plots K2SFF and PDCSAP flux over +-2 d around each predicted 2014 secondary, plus the K2SFF
arclength (roll) series, and runs the same box scan at 200 random off-eclipse epochs to
calibrate the false-alarm level of the box statistic.
"""
import os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from astropy.io import fits
import glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import ebmodel as M

H = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(H, "data")
nom = json.load(open(os.path.join(H, "fit_tess_nominal.json")))
k2 = json.load(open(os.path.join(H, "k2_results_nominal.json")))
P = nom["P"]
ecw14 = k2["ecw2014_assumed"]; esw = nom["esw"]
e14, w14 = np.hypot(ecw14, esw), np.arctan2(esw, ecw14)
dsec = M.secondary_offset(P, e14, w14)
tP = k2["primary_by_pipeline"]["k2sff"]["tmin"]

f = glob.glob(os.path.join(D, "mastDownload", "HLSP", "*k2sff*", "*.fits"))[0]
d = fits.open(f)[1].data
ts, fs, arc = d["T"] - 2167.0, d["FCOR"] / np.nanmedian(d["FCOR"]), d["ARCLENGTH"]
f = glob.glob(os.path.join(D, "mastDownload", "K2", "*", "ktwo*llc.fits"))[0]
d = fits.open(f)[1].data
m = (d["SAP_QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"])
tp, fp = d["TIME"][m] - 2167.0, d["PDCSAP_FLUX"][m] / np.nanmedian(d["PDCSAP_FLUX"][m])


def scan(tt, ff, tc, half=0.5):
    best = (-np.inf, 0, 0, 0)
    for wdt in (0.15, 0.25, 0.35):
        for c0 in np.arange(-half, half + 1e-9, 0.01):
            inn = np.abs(tt - tc - c0) < wdt / 2
            out = (np.abs(tt - tc - c0) > wdt / 2 + 0.05) & (np.abs(tt - tc - c0) < wdt / 2 + 0.4)
            if inn.sum() < 3 or out.sum() < 6:
                continue
            cc = np.polyfit(tt[out] - tc, ff[out], 1)
            dep = np.mean(np.polyval(cc, tt[inn] - tc) - ff[inn])
            err = np.std(ff[out] - np.polyval(cc, tt[out] - tc)) / np.sqrt(inn.sum())
            if dep / err > best[0]:
                best = (dep / err, dep, err, c0)
    return best

centres = [tP + dsec - P, tP + dsec]
OUT = {"windows": []}
fig, ax = plt.subplots(3, 1, figsize=(12, 10))
for k, tc in enumerate(centres):
    for tt, ff, lab in ((ts, fs, "K2SFF"), (tp, fp, "PDCSAP")):
        b = scan(tt, ff, tc)
        print("%s secondary centre %.3f: best dip %.5f (S/N %.1f) at offset %+.2f d" % (lab, tc, b[1], b[0], b[3]))
        OUT["windows"].append(dict(src=lab, centre=tc, dip=b[1], snr=b[0], offset=b[3], err=b[2]))
    a = ax[k]
    w = np.abs(ts - tc) < 2
    a.plot(ts[w] - tc, fs[w], "b.-", lw=0.3, ms=3, label="K2SFF")
    w = np.abs(tp - tc) < 2
    a.plot(tp[w] - tc, fp[w] - 0.001, "r.-", lw=0.3, ms=3, label="PDCSAP - 0.001")
    a.axvspan(-0.25, 0.25, color="y", alpha=0.2)
    a.set_title("predicted 2014 secondary at BTJD %.3f (window +-0.25 shaded)" % tc); a.legend(fontsize=8)
# false-alarm calibration on K2SFF away from eclipses
rng = np.random.default_rng(5)
snr = []; deps = []
ecl_times = [tP] + centres
for _ in range(300):
    tc = rng.uniform(ts.min() + 1, ts.max() - 1)
    if min(abs(tc - x) for x in ecl_times) < 1.2:
        continue
    bb = scan(ts, fs, tc)
    snr.append(bb[0]); deps.append(bb[1])
snr = np.array(snr)
print("K2SFF random-epoch box scan S/N: median %.1f, 90%% %.1f, 99%% %.1f, max %.1f (N=%d)" %
      (np.median(snr), np.percentile(snr, 90), np.percentile(snr, 99), snr.max(), len(snr)))
deps = np.array(deps)
print("K2SFF random-epoch best-box depth: median %.5f, 99%% %.5f, max %.5f" % (np.median(deps), np.percentile(deps, 99), deps.max()))
OUT.update(rand_snr_p50=np.median(snr), rand_snr_p99=np.percentile(snr, 99), rand_dep_p50=np.median(deps),
           rand_dep_p99=np.percentile(deps, 99), rand_dep_max=deps.max(), n_rand=len(snr))
json.dump(OUT, open(os.path.join(H, "k2_secondary_check.json"), "w"), indent=1, default=float)
ax[2].hist(snr, 40); ax[2].set_xlabel("best box S/N at random epochs (K2SFF)")
plt.tight_layout(); plt.savefig(os.path.join(H, "fig_k2_secondary_check.png"), dpi=80)
