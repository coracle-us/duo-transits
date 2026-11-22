"""Single-event search: detrend every sector, slide a box of the paper's duration across
all data, and list dips with depth comparable to the paper event.

Per sector the search source is TESS-SPOC if present, else QLP, else TESScut.
(QLP's S94+ clock offset of a few minutes is irrelevant for detection; timing uses TESScut.)
Writes out/search_<tic>.txt and plots out/search_<tic>.png
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import TARGETS, OUT, load_lc_files, load_tesscut_phot, detrend


def pick_sources(tic):
    lcs = load_lc_files(tic) + load_tesscut_phot(tic)
    rank = {"TESS-SPOC": 0, "QLP": 1, "TESScut": 2}
    by = {}
    for lc in lcs:
        s = lc["sector"]
        if len(lc["t"]) < 100:
            continue
        if s not in by or rank[lc["source"]] < rank[by[s]["source"]]:
            by[s] = lc
    return [by[s] for s in sorted(by)], lcs


def box_scan(t, f, e, dur, step=10 / 1440):
    """Return grid tc, depth, depth_err, snr (white), n_in for a box of width dur."""
    tcs = np.arange(t.min() + dur / 2, t.max() - dur / 2, step)
    w = 1 / e**2
    out = []
    for tc in tcs:
        inn = np.abs(t - tc) < dur / 2
        n = inn.sum()
        # require >= 60% coverage of the box
        cad = np.nanmedian(np.diff(t))
        if n < 0.6 * dur / cad:
            out.append((tc, np.nan, np.nan, n)); continue
        d = np.sum(w[inn] * (1 - f[inn])) / np.sum(w[inn])
        de = 1 / np.sqrt(np.sum(w[inn]))
        out.append((tc, d, de, n))
    a = np.array(out)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def red_factor(t, f, dur):
    """Ratio of the scatter of dur-binned flux to the white-noise expectation."""
    cad = np.nanmedian(np.diff(t))
    nb = max(1, int(round(dur / cad)))
    k = (len(f) // nb) * nb
    if k < 3 * nb:
        return 1.0
    s1 = np.std(f)
    sb = np.std(f[:k].reshape(-1, nb).mean(axis=1))
    return max(1.0, sb / (s1 / np.sqrt(nb)))


def run(tic):
    P = TARGETS[tic]
    dur, dep = P["dur"], P["depth"]
    secs, _ = pick_sources(tic)
    lines, allev = [], []
    fig, axs = plt.subplots(len(secs), 1, figsize=(14, 2.2 * len(secs)), squeeze=False)
    for ax, lc in zip(axs[:, 0], secs):
        t, f, e = lc["t"], lc["f"], lc["e"]
        win = max(0.75, 3 * dur)
        flat, tr = detrend(t, f, win)
        # iterative: mask strongest dips >0.4*dep, re-detrend
        _, d1, _, _ = box_scan(t, flat, e, dur)
        tcs, d, de, n = box_scan(t, flat, e, dur)
        m = np.zeros(len(t), bool)
        for tc_, dd in zip(tcs, d):
            if np.isfinite(dd) and dd > 0.4 * dep:
                m |= np.abs(t - tc_) < 0.75 * dur
        if m.any():
            flat, tr = detrend(t, f, win, mask=m)
            tcs, d, de, n = box_scan(t, flat, e, dur)
        # empirical error: scatter of out-of-dip points
        oot = ~m
        sig = 1.4826 * np.nanmedian(np.abs(flat[oot] - np.nanmedian(flat[oot])))
        beta = red_factor(t[oot], flat[oot], dur)
        cad = np.nanmedian(np.diff(t))
        de_emp = sig / np.sqrt(np.maximum(n, 1)) * beta
        snr = d / de_emp
        lines.append(f"S{lc['sector']:3d} {lc['source']:<9} n={len(t):6d} cad={cad*1440:5.1f}m "
                     f"rms={sig*1e6:6.0f}ppm beta={beta:4.2f} box-depth-err={np.nanmedian(de_emp)*1e6:5.0f}ppm "
                     f"span {t.min():.2f}-{t.max():.2f}")
        # local maxima events
        ok = np.isfinite(snr)
        idx = np.where(ok)[0]
        used = np.zeros(len(tcs), bool)
        order = idx[np.argsort(-snr[idx])]
        for i in order:
            if snr[i] < 5 or used[i]:
                continue
            used |= np.abs(tcs - tcs[i]) < dur
            allev.append((lc["sector"], lc["source"], tcs[i], d[i], de_emp[i], snr[i]))
        ax.plot(t, flat, ",k", alpha=0.5)
        nb = max(1, int(round((30 / 1440) / cad)))
        k = (len(t) // nb) * nb
        ax.plot(t[:k].reshape(-1, nb).mean(1), flat[:k].reshape(-1, nb).mean(1), ".", color="C0", ms=2)
        ax.axhline(1 - dep, color="r", lw=0.5)
        for ev in allev:
            if ev[0] == lc["sector"] and ev[3] > 0.4 * dep:
                ax.axvline(ev[2], color="C1", alpha=0.4)
        ax.set_ylim(1 - 2.2 * dep, 1 + 1.0 * dep)
        ax.text(0.01, 0.05, f"S{lc['sector']} {lc['source']}", transform=ax.transAxes)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"search_{tic}.png"), dpi=80)
    plt.close(fig)
    allev.sort(key=lambda x: -x[5])
    with open(os.path.join(OUT, f"search_{tic}.txt"), "w") as fh:
        fh.write(f"TIC {tic} paper: t0={P['t0']} depth={dep*1e6:.0f} dur={dur*24:.2f}h\n")
        fh.write("\n".join(lines) + "\n\nEvents SNR>5 (box of paper duration), sorted by SNR:\n")
        fh.write(" sec source        tc        depth[ppm]  err  SNR  depth/paper\n")
        for ev in allev[:40]:
            fh.write(f" {ev[0]:3d} {ev[1]:<9} {ev[2]:11.4f} {ev[3]*1e6:8.0f} {ev[4]*1e6:6.0f} {ev[5]:5.1f}  {ev[3]/dep:5.2f}\n")
    print(open(os.path.join(OUT, f"search_{tic}.txt")).read())


if __name__ == "__main__":
    for tic in ([int(a) for a in sys.argv[1:]] or TARGETS):
        run(tic)
