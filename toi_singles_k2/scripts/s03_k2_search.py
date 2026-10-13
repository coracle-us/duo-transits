"""
Step 2b: K2 epoch test for each TOI on a K2 star (rows of k2_overlap.csv).

Reductions searched (each normalised per campaign):
  K2PDC   K2 mission PDCSAP (k2lib.k2_lightcurves: lightkurve DEFAULT mask + 'possible thruster firing')
  K2SAP   K2 mission SAP flux, same mask (no systematics correction; the local polynomial baseline in
          the window test absorbs slow trends, not the 6-h roll sawtooth)
  EVEREST EVEREST v2 (QUALITY == 0)
  K2SFF   K2SFF FCOR (best aperture)
Search: k2lib.grid_windows_k2 (box test of the TOI's duration and depth at every D/4 step) and
k2lib.find_repeats (n_in >= 4, beta <= 2 or SNR >= 15, depth within pilot.DMIN..DMAX of the TOI depth,
baseline on both sides). The TESS epoch is the reference (no K2 window is 'near' it).
Candidates from any reduction are merged within one duration, and for each the window test is re-run in
every reduction at that time, so a hit lists depth/SNR/verdict per reduction. Candidates also record
the masked-cadence count within +-D of the event and the distance to the campaign edge.
Also reported: the 5 highest-SNR windows per reduction regardless of thresholds (sub-threshold look),
and, for TOIs with a maximum period, the K2 verdict for every alias whose predicted transit falls
in K2 data.

Outputs: toi_singles_k2/k2_hits.csv, toi_singles_k2/derived/k2_top_windows.csv,
         toi_singles_k2/derived/k2_alias_verdicts.csv, toi_singles_k2/plots/*.png
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "toi_singles_k2"
DATA = OUT / "data"
sys.path.insert(0, str(ROOT / "k2singles" / "scripts"))
import k2lib  # noqa: E402

k2lib.CACHE = DATA / "k2cache"
from lightkurve.utils import KeplerQualityFlags  # noqa: E402

MASK = KeplerQualityFlags.DEFAULT_BITMASK | k2lib.THRUSTER_POSSIBLE


def _segments(df):
    out = []
    for sec, g in df.groupby("sector"):
        g = g.sort_values("time").copy()
        g["flux"] = g.flux / np.nanmedian(g.flux)
        g["segment"] = sec * 1000 + np.concatenate([[0], np.cumsum(np.diff(g.time.to_numpy()) > 0.5)])
        out.append(g)
    return pd.concat(out, ignore_index=True).sort_values("time").reset_index(drop=True)


def load_sap(epic):
    from astropy.io import fits
    parts = []
    for fn in sorted((k2lib.CACHE / str(epic) / "dl").glob("ktwo*_llc.fits")):
        camp = int(fn.stem.split("-c")[1].split("_")[0])
        with fits.open(fn) as h:
            d = h[1].data
            t = np.asarray(d["TIME"], float) - 2167.0
            f = np.asarray(d["SAP_FLUX"], float)
            q = np.asarray(d["SAP_QUALITY"], int)
        ok = np.isfinite(t) & np.isfinite(f) & ((q & MASK) == 0)
        if ok.sum() > 50:
            parts.append(pd.DataFrame(dict(time=t[ok], flux=f[ok], err=np.nan, sector=1000 + camp, product="K2SAP")))
    return _segments(pd.concat(parts)) if parts else None


def raw_quality(epic):
    """(time, quality) for all cadences, to count masked cadences near an event."""
    from astropy.io import fits
    ts, qs = [], []
    for fn in sorted((k2lib.CACHE / str(epic) / "dl").glob("ktwo*_llc.fits")):
        with fits.open(fn) as h:
            d = h[1].data
            ts.append(np.asarray(d["TIME"], float) - 2167.0)
            qs.append(np.asarray(d["SAP_QUALITY"], int))
    if not ts:
        return np.array([]), np.array([], int)
    t, q = np.concatenate(ts), np.concatenate(qs)
    ok = np.isfinite(t)
    return t[ok], q[ok]


def load_hlsp(epic, author):
    import lightkurve as lk
    from astropy.io import fits
    ddir = DATA / "hlsp" / str(epic)
    ddir.mkdir(parents=True, exist_ok=True)
    try:
        sr = lk.search_lightcurve(f"EPIC {epic}", mission="K2", author=author)
    except Exception as e:  # noqa: BLE001
        print("  search failed", author, e)
        return None
    parts = []
    for i in range(len(sr)):
        try:
            lcf = sr[i].download(download_dir=str(ddir), quality_bitmask="none")
            fn = lcf.meta.get("FILENAME") or lcf.filename
        except Exception as e:  # noqa: BLE001
            print("  download failed", author, i, e)
            continue
        with fits.open(fn) as h:
            d = h[1].data
            camp = int(h[0].header.get("CAMPAIGN", -1))
            if author == "EVEREST":
                t = np.asarray(d["TIME"], float) - 2167.0
                f = np.asarray(d["FLUX"], float)
                # Kepler bits only (EVEREST sets its own bits >= 2**23, incl. an outlier flag that
                # catches real transits, e.g. EPIC 211924561 C16 in the K2-singles project)
                ok = (np.asarray(d["QUALITY"], int) & MASK & (2**23 - 1)) == 0
            else:  # K2SFF
                t = np.asarray(d["T"], float) - 2167.0
                f = np.asarray(d["FCOR"], float)
                ok = np.ones(len(t), bool)
        ok &= np.isfinite(t) & np.isfinite(f)
        if ok.sum() > 50:
            parts.append(pd.DataFrame(dict(time=t[ok], flux=f[ok], err=np.nan, sector=1000 + camp, product=author)))
    if not parts:
        return None
    df = pd.concat(parts).drop_duplicates(subset=["sector", "time"])
    return _segments(df)


def window_at(lc, T, D, delta):
    r = k2lib.test_window_k2(lc.time.to_numpy(), lc.flux.to_numpy(), T, D, delta, lc.segment.to_numpy())
    return r


def plot_hit(toi, epic, T, D, lcs, fn):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = [k for k in lcs if lcs[k] is not None]
    fig, axs = plt.subplots(len(names), 1, figsize=(8, 2.2 * len(names)), sharex=True)
    axs = np.atleast_1d(axs)
    for ax, k in zip(axs, names):
        lc = lcs[k]
        w = (lc.time - T).abs() < max(2.0, 4 * D)
        ax.plot(lc.time[w] - T, lc.flux[w], ".", ms=3, color="k")
        ax.axvspan(-D / 2, D / 2, color="tab:red", alpha=0.15)
        ax.set_ylabel(k)
    axs[-1].set_xlabel(f"BTJD - {T:.3f} (d)")
    axs[0].set_title(f"TOI-{toi} / EPIC {epic}: K2 candidate at BTJD {T:.3f} (shaded = TOI duration)")
    fig.tight_layout()
    fig.savefig(fn, dpi=110)
    plt.close(fig)


def main():
    ov = pd.read_csv(OUT / "k2_overlap.csv")
    ov = ov[ov.campaigns.notna() & (ov.campaigns.astype(str) != "")]
    hits, tops, aliasv = [], [], []
    (OUT / "plots").mkdir(exist_ok=True)
    for r in ov.itertuples():
        epic, D, delta = int(r.EPIC), r.dur_h / 24, r.depth_ppm / 1e6
        print(f"TOI-{r.TOI} EPIC {epic}: D={r.dur_h:.2f} h, depth={r.depth_ppm:.0f} ppm, C{r.campaigns}", flush=True)
        pdc = k2lib.k2_lightcurves(epic)
        lcs = {"K2PDC": _segments(pdc) if pdc is not None else None, "K2SAP": load_sap(epic),
               "EVEREST": load_hlsp(epic, "EVEREST"), "K2SFF": load_hlsp(epic, "K2SFF")}
        tq, qq = raw_quality(epic)
        cands = []
        for name, lc in lcs.items():
            if lc is None:
                continue
            g = k2lib.grid(lc, D, delta, r.epoch_btjd)
            if "snr" in g:
                top = g[g.snr.notna()].sort_values("snr", ascending=False).head(5)
                for _, w in top.iterrows():
                    tops.append(dict(TOI=r.TOI, EPIC=epic, product=name, T=round(w.T_pred, 4), snr=round(w.snr, 2),
                                     depth_ppm=w.depth * 1e6, verdict=w.verdict, beta=round(w.beta, 2),
                                     n_in=w.n_in, cover=round(w.cover, 2)))
            reps = k2lib.find_repeats(g, delta, D)
            print(f"  {name}: {len(lc)} cadences, {len(reps)} candidate events", flush=True)
            for x in reps:
                cands.append(dict(x, product=name))
            if np.isfinite(r.max_period_d) and name == "K2PDC":
                Pmax = r.max_period_d
                t = lc.time.to_numpy()
                for n in range(1, int(Pmax // 5) + 1):
                    P = Pmax / n
                    k = np.arange(np.floor((t.min() - r.epoch_btjd) / P) - 1, np.ceil((t.max() - r.epoch_btjd) / P) + 2)
                    for T in r.epoch_btjd + k * P:
                        if T < t.min() - D or T > t.max() + D:
                            continue
                        w = window_at(lc, T, D, delta)
                        if w["verdict"] == "no-data":
                            continue
                        aliasv.append(dict(TOI=r.TOI, EPIC=epic, n=n, P=round(P, 5), T_pred=round(T, 4),
                                           verdict=w["verdict"], depth_ppm=w.get("depth", np.nan) * 1e6,
                                           snr=w.get("snr", np.nan)))
        # merge candidates across reductions
        merged = []
        for c in sorted(cands, key=lambda c: -c["snr"]):
            m = next((m for m in merged if abs(m["T"] - c["T"]) < D), None)
            if m is None:
                merged.append(dict(T=c["T"], found_in=[c["product"]], best_snr=c["snr"]))
            elif c["product"] not in m["found_in"]:
                m["found_in"].append(c["product"])
        for m in merged:
            T = m["T"]
            row = dict(TOI=r.TOI, TIC=r.TIC, EPIC=epic, T_btjd=T, campaign=None, toi_depth_ppm=round(r.depth_ppm),
                       toi_dur_h=round(r.dur_h, 2), found_in=" ".join(m["found_in"]), n_found=len(m["found_in"]))
            for name, lc in lcs.items():
                if lc is None:
                    row[f"{name}_verdict"] = "no-product"
                    continue
                w = window_at(lc, T, D, delta)
                row[f"{name}_verdict"] = w["verdict"]
                row[f"{name}_depth_ppm"] = w.get("depth", np.nan) * 1e6
                row[f"{name}_snr"] = round(w.get("snr", np.nan), 1) if "snr" in w else np.nan
                if name == "K2PDC" and "depth" in w:
                    k = np.argmin(np.abs(lc.time.to_numpy() - T))
                    row["campaign"] = int(lc.sector.iloc[k]) - 1000
                    seg = lc[lc.sector == lc.sector.iloc[k]]
                    row["d_edge_d"] = round(float(min(T - seg.time.min(), seg.time.max() - T)), 2)
            near = np.abs(tq - T) < D
            row["n_cad_within_D"] = int(near.sum())
            row["n_masked_within_D"] = int(((qq[near] & MASK) != 0).sum())
            row["n_thruster_within_D"] = int(((qq[near] & (1048576 | 524288 | 32)) != 0).sum())
            hits.append(row)
            plot_hit(r.TOI, epic, T, D, lcs, OUT / "plots" / f"toi{r.TOI}_epic{epic}_T{T:.2f}.png")
    pd.DataFrame(hits).to_csv(OUT / "k2_hits.csv", index=False)
    pd.DataFrame(tops).to_csv(OUT / "derived" / "k2_top_windows.csv", index=False)
    pd.DataFrame(aliasv).to_csv(OUT / "derived" / "k2_alias_verdicts.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(hits).to_string() if hits else "no K2 candidate events")


if __name__ == "__main__":
    main()
