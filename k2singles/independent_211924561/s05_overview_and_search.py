"""Overview plots of every dataset and a single-event matched-filter search.

For each light curve: running-median detrend (1.0 d window, robust), then slide a box of
duration D (several trial durations) in steps of 0.01 d and compute depth and depth/sigma, where
sigma comes from the local out-of-box scatter scaled by sqrt(N_in) and inflated by a red-noise
factor beta measured empirically from binned residual scatter.
"""
import os, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lcio

DATASETS = {
    # name: loader
    "K2C05_EVEREST_noQmask": None,  # special-cased below
    "K2C05_K2SFF": lambda: lcio.k2_k2sff(5),
    "K2C16_EVEREST": lambda: lcio.k2_everest(16),
    "K2C16_K2SFF": lambda: lcio.k2_k2sff(16),
    "K2C18_EVEREST": lambda: lcio.k2_everest(18),
    "K2C18_K2SFF": lambda: lcio.k2_k2sff(18),
    "S44_SPOC2min": lambda: lcio.tess_spoc2min(44),
    "S45_SPOC2min": lambda: lcio.tess_spoc2min(45),
    "S46_SPOC2min": lambda: lcio.tess_spoc2min(46),
    "S72_SPOC2min": lambda: lcio.tess_spoc2min(72),
    "S44_TESSSPOC": lambda: lcio.tess_tessspoc(44),
    "S45_TESSSPOC": lambda: lcio.tess_tessspoc(45),
    "S46_TESSSPOC": lambda: lcio.tess_tessspoc(46),
    "S72_TESSSPOC": lambda: lcio.tess_tessspoc(72),
    "S1751_TESSSPOC": lambda: lcio.tess_tessspoc(1751),
    "S44_QLP": lambda: lcio.tess_qlp(44),
    "S45_QLP": lambda: lcio.tess_qlp(45),
    "S46_QLP": lambda: lcio.tess_qlp(46),
    "S72_QLP": lambda: lcio.tess_qlp(72),
    "S44_TARS": lambda: lcio.tess_tars(44),
    "S45_TARS": lambda: lcio.tess_tars(45),
    "S46_TARS": lambda: lcio.tess_tars(46),
    "S72_TARS": lambda: lcio.tess_tars(72),
}


def load_everest_noq(c):
    from astropy.io import fits
    d = fits.getdata(lcio._find(f"hlsp_everest_k2_llc_211924561-c{c:02d}_kepler_v2.0_lc.fits"), 1)
    t = d["TIME"] - 2167.0
    # keep everything except thruster-firing (2^20) and NaN; EVEREST's own outlier bit (2^24) is NOT applied
    q = (d["QUALITY"] & (1 << 20)) == 0
    return lcio._norm(t, np.asarray(d["FLUX"], float), np.asarray(d["FRAW_ERR"], float), q)


def bin_to(t, f, dt):
    """Bin to cadence dt (days) to equalise high-cadence data."""
    if np.median(np.diff(t)) >= dt * 0.9:
        return t, f
    edges = np.arange(t.min(), t.max() + dt, dt)
    idx = np.digitize(t, edges)
    tb, fb = [], []
    for i in np.unique(idx):
        s = idx == i
        if s.sum() >= 2:
            tb.append(t[s].mean()); fb.append(np.median(f[s]))
    return np.array(tb), np.array(fb)


def box_search(t, r, durations=(0.06, 0.09, 0.12, 0.16, 0.20), step=0.005):
    sig_pt = 1.4826 * np.median(np.abs(r - np.median(r)))
    cad = np.median(np.diff(t))
    out = []
    grid = np.arange(t.min(), t.max(), step)
    for D in durations:
        # red-noise beta: scatter of means in bins of length D vs white expectation
        nb = max(int(round(D / cad)), 1)
        nbins = len(r) // nb
        means = r[: nbins * nb].reshape(nbins, nb).mean(1)
        beta = max(1.0, np.std(means) / (sig_pt / np.sqrt(nb)))
        for tc in grid:
            inb = np.abs(t - tc) < D / 2
            n = inb.sum()
            if n < max(3, 0.6 * D / cad):  # require 60% coverage of the box
                continue
            dep = -np.mean(r[inb])
            snr = dep / (beta * sig_pt / np.sqrt(n))
            out.append((tc, D, dep, snr, n, beta))
    return np.array(out), sig_pt


results = {}
fig, axes = plt.subplots(len(DATASETS), 1, figsize=(16, 2.0 * len(DATASETS)))
for ax, (name, loader) in zip(axes, DATASETS.items()):
    d = load_everest_noq(5) if name == "K2C05_EVEREST_noQmask" else loader()
    t, f = d["t"], d["f"]
    if name.endswith("SPOC2min") or "K2C18" in name and False:
        t, f = bin_to(t, f, 10 / 1440.)
    if name.startswith("S72") or name.startswith("S1751"):
        t, f = bin_to(t, f, 10 / 1440.)
    o = np.argsort(t); t, f = t[o], f[o]
    flat, trend = lcio.flatten(t, f, window=1.0)
    r = flat - 1
    good = np.isfinite(r)
    t, r, f, trend = t[good], r[good], f[good], trend[good]
    s, sig = box_search(t, r)
    ax.plot(t, 1e3 * r, ".", ms=1.5, color="k")
    ax.set_ylim(-4, 2)
    ax.text(0.005, 0.05, f"{name}  pt-rms={sig*1e6:.0f} ppm", transform=ax.transAxes, fontsize=8)
    # top candidates: local maxima in snr, merge within 0.3 d
    cand = []
    if len(s):
        order = np.argsort(-s[:, 3])
        for k in order:
            tc = s[k, 0]
            if any(abs(tc - c[0]) < 0.3 for c in cand):
                continue
            cand.append(tuple(s[k]))
            if len(cand) >= 5:
                break
    for c in cand:
        ax.axvline(c[0], color="r" if c[3] > 7 else "orange", alpha=0.5)
    results[name] = dict(rms_ppm=sig * 1e6, tmin=float(t.min()), tmax=float(t.max()), n=int(len(t)),
                         top=[dict(tc=round(c[0], 4), dur=c[1], depth_ppt=round(1e3 * c[2], 3), snr=round(c[3], 2), n=int(c[4]), beta=round(c[5], 2)) for c in cand])
    print(name, f"rms={sig*1e6:.0f}ppm  span={t.min():.2f}-{t.max():.2f}")
    for c in cand:
        print(f"    tc={c[0]:.3f} D={c[1]:.2f} depth={1e3*c[2]:.2f} ppt snr={c[3]:.1f} beta={c[5]:.2f}")
    np.save(os.path.join(lcio.HERE, "cache", f"{name}.npy"), np.vstack([t, f, trend, r]))
plt.tight_layout()
plt.savefig(os.path.join(lcio.HERE, "fig_overview_all_detrended.png"), dpi=80)
json.dump(results, open(os.path.join(lcio.HERE, "search_top_events.json"), "w"), indent=1)
