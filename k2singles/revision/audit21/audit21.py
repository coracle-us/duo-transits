"""
Audit of the 21 catalogued K2 single events (LaCourse & Jacobs 2018) that are NOT 'supported' at their catalogued
epoch in k2singles/revision/catalogued_event_check.csv. For each event, every available K2 product is examined in a
+-1.5 d window around the catalogued epoch:

  PDC_q     K2 mission PDCSAP_FLUX with the repeat-search mask (DEFAULT | possible thruster firing) = what the check used
  PDC_all   K2 mission PDCSAP_FLUX, all finite cadences (no quality mask)
  SAP_q     K2 mission SAP_FLUX with the same mask
  SAP_all   K2 mission SAP_FLUX, all finite cadences (no quality mask)
  EVEREST   EVEREST v2 FLUX (cadences with K2 hard flags in the repeat-search mask dropped; EVEREST's own
            outlier/mask bits 23 and 25, which EVEREST sets on in-transit cadences, are kept)
  K2SFF     K2SFF v1 BESTAPER FCOR (thruster-firing cadences already removed by K2SFF)
  K2SFF_raw K2SFF v1 BESTAPER FRAW (K2SFF aperture photometry before the arclength correction)
  K2SC      K2SC PDC extension, flux - trposi + median(trposi) (only where the file exists)

For each product: k2lib.test_window_k2 at the catalogued epoch (depth, depth/catalogued, SNR, verdict), a scan of
the same test over offsets -1.5..+1.5 d in steps of D/8 (to detect epoch errors), and the in-window cadence census
(how many cadences inside +-0.35 D of the epoch are present / NaN / quality-flagged in the K2 LC file).
One diagnostic plot per event goes to plots/. Metrics go to audit21_metrics.csv (one row per event x product).
Classification is NOT done here; it is done by hand in classify21.py from these metrics and the plots.

Usage:  python k2singles/revision/audit21/audit21.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "k2singles" / "scripts"))
import k2lib  # noqa: E402

warnings.filterwarnings("ignore")
DATA = HERE / "data"
PLOTS = HERE / "plots"
BK = 2454833 - 2457000          # BKJD -> BTJD
HALF_WIN = 1.5
PRODUCTS = ["PDC_q", "PDC_all", "SAP_q", "SAP_all", "EVEREST", "K2SFF", "K2SFF_raw", "K2SC"]


def targets():
    d = pd.read_csv(ROOT / "k2singles" / "revision" / "catalogued_event_check.csv")
    c = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv")
    c["event_idx"] = c.groupby("EPIC").cumcount()
    bad = d[d.verdict != "supported"][["EPIC", "event_idx", "verdict", "depth_ratio", "snr"]]
    return bad.merge(c, on=["EPIC", "event_idx"])


def _seg(t):
    return np.concatenate([[0], np.cumsum(np.diff(t) > 0.5)]) if len(t) else np.zeros(0, int)


def _norm(t, f, flag, T, D):
    """Sort, keep finite, normalise by the median out of the event window (+-D) within +-3 d."""
    o = np.argsort(t)
    t, f, flag = t[o], f[o], flag[o]
    ok = np.isfinite(t) & np.isfinite(f)
    t, f, flag = t[ok], f[ok], flag[ok]
    near = (np.abs(t - T) < 3) & (np.abs(t - T) > D)
    med = np.nanmedian(f[near]) if near.sum() > 10 else np.nanmedian(f)
    return t, f / med, flag


def load_products(epic, camp, T, D):
    """Return {product: (t, f, flagged_bool)} restricted to |t - T| < 4 d, plus a K2-LC cadence census."""
    out, census = {}, {}
    cdir = ROOT / "k2singles" / "cache" / str(epic) / "dl"
    fns = sorted(cdir.glob(f"ktwo{epic}-c{camp:02d}*_llc.fits")) + sorted(cdir.glob(f"ktwo{epic}-c{camp}?*_llc.fits"))
    fns = list(dict.fromkeys(fns))
    from astropy.io import fits
    from lightkurve.utils import KeplerQualityFlags
    mask = KeplerQualityFlags.DEFAULT_BITMASK | k2lib.THRUSTER_POSSIBLE
    for fn in fns:
        with fits.open(fn) as h:
            d = h[1].data
            t = np.asarray(d["TIME"], float) + BK
            sap, pdc = np.asarray(d["SAP_FLUX"], float), np.asarray(d["PDCSAP_FLUX"], float)
            q = np.asarray(d["SAP_QUALITY"], int)
        w = np.isfinite(t) & (np.abs(t - T) < 4)
        if not w.any():
            continue
        t, sap, pdc, q = t[w], sap[w], pdc[w], q[w]
        bad = (q & mask) != 0
        out["PDC_q"] = _norm(t[~bad], pdc[~bad], np.zeros((~bad).sum(), bool), T, D)
        out["PDC_all"] = _norm(t, pdc, bad, T, D)
        out["SAP_q"] = _norm(t[~bad], sap[~bad], np.zeros((~bad).sum(), bool), T, D)
        out["SAP_all"] = _norm(t, sap, bad, T, D)
        core = np.abs(t - T) < 0.35 * D
        half = np.abs(t - T) < 0.5 * D
        bits = sorted({b for qq in q[half] for b in range(24) if qq & (1 << b)})
        census = dict(lc_file=fn.name, n_cad_core=int(core.sum()), n_pdc_finite_core=int(np.isfinite(pdc[core]).sum()),
                      n_sap_finite_core=int(np.isfinite(sap[core]).sum()),
                      n_flagged_core=int(bad[core].sum()), n_pdc_used_core=int((np.isfinite(pdc) & ~bad)[core].sum()),
                      qual_bits_in_half=";".join(str(b + 1) for b in bits))
    hl = DATA / "mastDownload" / "HLSP"
    for p in hl.glob(f"hlsp_everest_k2_llc_{epic}-c{camp:02d}*/*.fits"):
        with fits.open(p) as h:
            a = h["ARRAYS"].data
            t = np.asarray(a["TIME"], float) + BK
            f, q = np.asarray(a["FLUX"], float), np.asarray(a["QUALITY"], np.int64)
        w = np.abs(t - T) < 4
        # flag only K2 hard flags (repeat-search mask); EVEREST's own bits 23 (outlier) and 25 (masked) are set on
        # in-transit cadences by EVEREST itself and are NOT treated as bad here
        out["EVEREST"] = _norm(t[w], f[w], (q[w] & mask) != 0, T, D)
    for p in hl.glob(f"hlsp_k2sff_k2_lightcurve_{epic}-c{camp:02d}*/*.fits"):
        with fits.open(p) as h:
            a = h["BESTAPER"].data
            t = np.asarray(a["T"], float) + BK
            f = np.asarray(a["FCOR"], float)
            fr = np.asarray(a["FRAW"], float)
        w = np.abs(t - T) < 4
        out["K2SFF"] = _norm(t[w], f[w], np.zeros(w.sum(), bool), T, D)
        out["K2SFF_raw"] = _norm(t[w], fr[w], np.zeros(w.sum(), bool), T, D)
    for p in hl.glob(f"hlsp_k2sc_k2_llc_{epic}-c{camp:02d}*/*.fits"):
        with fits.open(p) as h:
            a = h["PDC"].data
            t = np.asarray(a["time"], float) + BK
            f = np.asarray(a["flux"], float) - np.asarray(a["trposi"], float) + np.nanmedian(a["trposi"])
            q = np.asarray(a["quality"], int)
        w = np.abs(t - T) < 4
        out["K2SC"] = _norm(t[w], f[w], (q[w] & mask) != 0, T, D)
    return out, census


def window(t, f, T, D, delta):
    if len(t) < 5:
        return {"verdict": "no-data"}
    return k2lib.test_window_k2(t, f, T, D, delta, _seg(t))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    PLOTS.mkdir(exist_ok=True)
    rows = []
    for _, r in targets().iterrows():
        epic, camp, T = int(r.EPIC), int(r.campaign), float(r.epoch_btjd)
        D, delta = r.dur_h / 24, r.depth_ppm * 1e-6
        prods, census = load_products(epic, camp, T, D)
        offs = np.arange(-HALF_WIN, HALF_WIN + 1e-9, D / 8)
        fig, axs = plt.subplots(len(PRODUCTS), 1, figsize=(9, 1.9 * len(PRODUCTS)), sharex=True)
        for ax, name in zip(axs, PRODUCTS):
            row = dict(EPIC=epic, event_idx=int(r.event_idx), check_verdict=r.verdict, product=name, campaign=camp,
                       epoch=T, dur_h=r.dur_h, depth_ppm=r.depth_ppm, comment=r.comment, **census)
            if name not in prods:
                ax.text(0.5, 0.5, f"{name}: not available", transform=ax.transAxes, ha="center")
                rows.append({**row, "available": False})
                continue
            t, f, flag = prods[name]
            use = ~flag if name in ("PDC_all", "SAP_all", "EVEREST") else np.ones(len(t), bool)
            # PDC_all / SAP_all: the test uses ALL finite cadences, flagged included. EVEREST: K2-flagged cadences
            # are dropped from the test (they carry EVEREST-unconstrained outliers) but still plotted as red x.
            tt, ff = (t[use], f[use]) if name == "EVEREST" else (t, f)
            w0 = window(tt, ff, T, D, delta)
            scan = [window(tt, ff, T + o, D, delta) for o in offs]
            dep = np.array([s.get("depth", np.nan) for s in scan])
            snr = np.array([s.get("snr", np.nan) for s in scan])
            sup = np.array([s["verdict"] == "supported" for s in scan])
            # deepest supported window in the scan (depth, not SNR, so a deep event is not beaten by a
            # longer shallow systematic)
            if sup.any():
                k = np.nanargmax(np.where(sup, dep, -np.inf))
                best = dict(scan_best_offset_h=offs[k] * 24, scan_best_ratio=dep[k] / delta, scan_best_snr=snr[k],
                            scan_n_supported=int(sup.sum()))
            else:
                best = dict(scan_best_offset_h=np.nan, scan_best_ratio=np.nanmax(dep) / delta if np.isfinite(dep).any()
                            else np.nan, scan_best_snr=np.nan, scan_n_supported=0)
            core = np.abs(t - T) < 0.35 * D
            # model-free "minimum depth": local median baseline minus mean of the 3 lowest unflagged cadences
            # inside +-D/2 (meaningless for SAP, which carries the 6-h roll sawtooth)
            bl = use & (np.abs(t - T) > D / 2 + 0.05) & (np.abs(t - T) < D / 2 + 0.6)
            inw = np.sort(f[use & (np.abs(t - T) < D / 2)])
            min3 = (np.median(f[bl]) - np.mean(inw[:3])) if (bl.sum() > 5 and len(inw) >= 3) else np.nan
            rows.append({**row, "available": True, "n_core": int(core.sum()), "n_core_flagged": int(flag[core].sum()), "min3_depth": min3, "min3_ratio": min3 / delta,
                         "verdict": w0["verdict"], "depth": w0.get("depth", np.nan),
                         "ratio": w0.get("depth", np.nan) / delta, "snr": w0.get("snr", np.nan),
                         "beta": w0.get("beta", np.nan), "sigma_white": w0.get("sigma_white", np.nan), **best})
            win = np.abs(t - T) < HALF_WIN + 0.3
            ax.plot(t[win & use] - T, f[win & use], ".", ms=3, color="k")
            if (~use & win).any():
                ax.plot(t[win & ~use] - T, f[win & ~use], "x", ms=4, color="tab:red", label="quality-flagged")
                ax.legend(loc="lower left", fontsize=7)
            ax.axvspan(-D / 2, D / 2, color="tab:blue", alpha=0.12)
            ax.axhline(1 - delta, ls="--", color="tab:orange", lw=1)
            ax.axhline(1 - 0.5 * delta, ls=":", color="tab:orange", lw=1)
            ax.axhline(1, color="0.6", lw=0.8)
            ax.plot(offs, 1 - dep, "-", color="tab:green", lw=1, alpha=0.7)
            if np.isfinite(best["scan_best_offset_h"]):
                ax.axvline(best["scan_best_offset_h"] / 24, color="tab:purple", lw=1, ls="-.")
            fw = f[win & use]
            lo = min(1 - 1.3 * delta, np.nanpercentile(fw, 1)) if len(fw) else 1 - 1.3 * delta
            hi = np.nanpercentile(fw, 99) if len(fw) else 1.01
            ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.1 * (hi - lo))
            ax.set_ylabel(name, fontsize=8)
            ax.text(0.01, 0.97, f"at epoch: {w0['verdict']}, ratio {w0.get('depth', np.nan) / delta:.2f}, "
                    f"SNR {w0.get('snr', np.nan):.1f}; deepest supported scan offset "
                    f"{best['scan_best_offset_h']:.1f} h", transform=ax.transAxes, va="top", fontsize=7)
        axs[-1].set_xlabel(f"BTJD - {T:.4f} (catalogued epoch) [d]")
        axs[0].set_title(f"EPIC {epic} C{camp}: depth {r.depth_ppm:.0f} ppm, D {r.dur_h:.2f} h, "
                         f"comment: {r.comment if isinstance(r.comment, str) else '-'}; check: {r.verdict}\n"
                         "blue = catalogued duration; orange dashed/dotted = 1 and 0.5 x catalogued depth; "
                         "green = scan depth; purple = deepest supported scan window", fontsize=8)
        axs[0].set_xlim(-HALF_WIN, HALF_WIN)
        fig.tight_layout()
        fig.savefig(PLOTS / f"{epic}_{int(r.event_idx)}.png", dpi=110)
        plt.close(fig)
        print(epic, "done", flush=True)
    m = pd.DataFrame(rows)
    m.to_csv(HERE / "audit21_metrics.csv", index=False)
    pd.set_option("display.width", 250)
    print(m[["EPIC", "product", "n_core", "n_core_flagged", "verdict", "ratio", "min3_ratio", "snr", "beta", "scan_best_offset_h",
             "scan_best_ratio", "scan_best_snr"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
