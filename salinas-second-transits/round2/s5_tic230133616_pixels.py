#!/usr/bin/env python3
"""Round 2, section 5: pixel-level check of TIC 230133616 at P = 23.7535 d.

- Gaia DR3 neighbours within 2 arcmin brighter than T + 6.
- Every predicted S47-S60 epoch (from the alias step's epoch file): covered and
  flat, in a gap, or showing a shallower dip.
- 15x15 TESScut cutouts for the S40 paper transit, three S47-S60 predicted
  epochs and three S73+ transits: out-of-transit image, difference image
  (out minus in), its flux-weighted centroid, and the in- vs out-of-transit
  centroid shift of the aperture.

Run after the alias step (needs work/aliases/230133616_*_epochs.csv).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import st2core as core  # noqa: E402
import tesscut_phot as tp  # noqa: E402

TIC = 230133616
WORK, OUT = HERE / "work", HERE / "results" / "round2"
PIX_ARCSEC = 21.0


def tmag_from_gaia(g, bp_rp):
    """Stassun et al. 2019 (TIC-8) G -> T relation; G - 0.43 without a colour."""
    c = np.where(np.isfinite(bp_rp), np.clip(bp_rp, -0.2, 3.0), np.nan)
    t = g - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(t), t, g - 0.43)


def gaia_cone(ra, dec, radius_arcmin=2.0):
    """Gaia DR3 sources in a cone: ESA archive, or the CDS VizieR copy if ESA fails."""
    try:
        return gaia_esa(ra, dec, radius_arcmin)
    except Exception as e:
        print(f"ESA Gaia archive failed ({e}); using VizieR I/355/gaiadr3")
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astroquery.vizier import Vizier
    v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag", "BP-RP", "Plx", "RUWE", "_r"], row_limit=-1)
    r = v.query_region(SkyCoord(ra, dec, unit="deg"), radius=radius_arcmin * u.arcmin, catalog="I/355/gaiadr3")[0]
    r = r.to_pandas().rename(columns={"Source": "source_id", "RA_ICRS": "ra", "DE_ICRS": "dec", "Gmag": "phot_g_mean_mag",
                                      "BP-RP": "bp_rp", "Plx": "parallax", "RUWE": "ruwe", "_r": "sep_arcsec"})
    r["sep_arcsec"] *= 60
    return r.sort_values("phot_g_mean_mag")


def gaia_esa(ra, dec, radius_arcmin):
    from astroquery.gaia import Gaia
    q = f"""SELECT source_id, ra, dec, phot_g_mean_mag, bp_rp, parallax, ruwe,
            DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec})) * 3600 AS sep_arcsec
            FROM gaiadr3.gaia_source
            WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius_arcmin / 60}))
            ORDER BY phot_g_mean_mag"""
    return Gaia.launch_job(q).get_results().to_pandas()


def gaia_neighbours(ra, dec, tmag):
    r = gaia_cone(ra, dec)
    r["tmag_est"] = tmag_from_gaia(r["phot_g_mean_mag"].to_numpy(float), r["bp_rp"].to_numpy(float))
    r["dT"] = r["tmag_est"] - tmag
    r["max_depth_if_source_ppm"] = np.nan  # filled by caller from the observed depth
    return r[r["tmag_est"] < tmag + 6].reset_index(drop=True)


def epoch_cutout(ra, dec, sector, tc, dur, gaia):
    tpf = tp.cutout(ra, dec, sector, WORK / "tesscut")
    ap, bkg_mask, edge_px = tp.apertures(tpf, ra, dec)
    fl = tpf.flux.value
    bkg = np.nanmedian(fl[:, bkg_mask], axis=1)
    ok = (tpf.quality == 0) & np.isfinite(tpf.time.value) & np.isfinite(bkg)
    sub = fl - bkg[:, None, None]
    diff, img_out, n_in, n_out = tp.difference_image(tpf, tc, dur, ok)
    x0, y0 = tp.target_pixel(tpf, ra, dec)
    row = dict(sector=int(sector), tc=tc, n_in=n_in, n_out=n_out, target_col=x0, target_row=y0,
               target_to_detector_edge_px=edge_px, n_ap=int(ap.sum()))
    if diff is None:
        row["note"] = "too few in- or out-of-transit cadences"
        return row, None
    diff = diff  # background cancels in the difference
    img_out = img_out - np.nanmedian(bkg[ok])
    # aperture photometry depth and its noise from out-of-transit scatter
    t = tpf.time.value
    f_ap = np.nansum(sub[:, ap], axis=1)
    dt = np.abs(t - tc)
    inn, out = ok & (dt < 0.4 * dur), ok & (dt > 0.75 * dur) & (dt < 2.5 * dur)
    f0 = np.nanmedian(f_ap[out])
    depth = 1 - np.nanmean(f_ap[inn]) / f0
    depth_err = core.robust_std(f_ap[out] / f0) / np.sqrt(inn.sum())
    # centroids: of the difference image (where the flux went missing) within 3 px
    yy, xx = np.indices(diff.shape)
    near = np.hypot(xx - x0, yy - y0) <= 3
    dcx, dcy = tp.centroid(diff, near)
    # in- vs out-of-transit flux-weighted centroid of the aperture
    def cen(m):
        im = np.nanmean(sub[m], axis=0)
        return tp.centroid(im, ap)
    icx, icy = cen(inn)
    ocx, ocy = cen(out)
    peak = np.unravel_index(np.nanargmax(np.where(near, diff, -np.inf)), diff.shape)
    row.update(ap_depth_ppm=depth * 1e6, ap_depth_err_ppm=depth_err * 1e6, ap_depth_snr=depth / depth_err,
               diff_centroid_col=dcx, diff_centroid_row=dcy,
               diff_centroid_offset_arcsec=np.hypot(dcx - x0, dcy - y0) * PIX_ARCSEC,
               diff_peak_col=int(peak[1]), diff_peak_row=int(peak[0]),
               incen_minus_outcen_col_mpix=(icx - ocx) * 1e3, incen_minus_outcen_row_mpix=(icy - ocy) * 1e3)
    gx, gy = tpf.wcs.all_world2pix(gaia[["ra", "dec"]].to_numpy(), 0).T if len(gaia) else ([], [])
    return row, dict(img_out=img_out, diff=diff, ap=ap, x0=x0, y0=y0, gx=gx, gy=gy, t=t[ok],
                     f=f_ap[ok] / f0, tc=tc)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cov = pd.read_csv(WORK / "coverage.csv").set_index("tic").loc[TIC]
    ev = pd.read_csv(WORK / "events_all.csv")
    ev = ev[(ev.tic == TIC) & (ev.passes == True)]
    t1, dur = ev.t1.iloc[0], ev.dur_ref_h.iloc[0] / 24
    depth_new = ev.depth_ppm.median()

    gaia = gaia_neighbours(cov.ra, cov.dec, cov.tmag)
    # a neighbour fainter by dT can make an observed depth d only if d * 10^(0.4 dT) <= 1
    gaia["needed_own_depth"] = depth_new * 1e-6 * 10 ** (0.4 * gaia["dT"].clip(lower=0))
    gaia.to_csv(OUT / "tic230133616_gaia_neighbours.csv", index=False)
    print(f"Gaia DR3 within 2' brighter than T+6 (T={cov.tmag:.2f}): {len(gaia)}")
    print(gaia[["source_id", "sep_arcsec", "phot_g_mean_mag", "bp_rp", "tmag_est", "dT", "ruwe",
                "needed_own_depth"]].to_string(index=False))

    ep = pd.concat(pd.read_csv(p) for p in (WORK / "aliases").glob(f"{TIC}*_epochs.csv"))
    hyp = ep[ep.period_source.str.startswith("hypothesis")]
    P = hyp.period_d.iloc[0]
    k_all = np.arange(int(np.floor((1950 - t1) / P)), int(np.ceil((3300 - t1) / P)) + 1)
    pred = pd.DataFrame(dict(k=k_all, epoch_btjd=t1 + k_all * P))
    old = pred[(pred.epoch_btjd > 2600) & (pred.epoch_btjd < 3100)].merge(
        hyp[["k", "sector", "pipeline", "coverage", "edge_dist_d", "box_depth_ppm", "box_depth_err_ppm",
             "verdict", "verdict_any_coverage"]], on="k", how="left")
    old["verdict"] = old["verdict"].fillna("in a gap (no data)")
    old.to_csv(OUT / "tic230133616_S47-S60_epochs.csv", index=False)
    print(f"\nPredicted S47-S60 epochs at P = {P:.6f} d (box depth vs {depth_new:.0f} ppm in S73+):")
    print(old.to_string(index=False))

    # cutout epochs: paper transit, the three best-covered S47-S60 epochs, three S73+ transits
    sel = [("paper transit", 40, t1)]
    cand = old.dropna(subset=["sector"]).sort_values(["coverage", "edge_dist_d"], ascending=False)
    sel += [(f"predicted ({r.verdict})", int(r.sector), r.epoch_btjd) for r in cand.head(3).itertuples()]
    clean = ev[(ev.edge_dist_d >= 0.5) & (ev.frac >= 0.8)].sort_values("snr", ascending=False)
    sel += [(f"new transit SNR {r.snr:.1f}", int(r.sector), r.t2) for r in clean.head(3).itertuples()]
    rows, pans = [], []
    for lab, s, tc in sel:
        try:
            row, pan = epoch_cutout(cov.ra, cov.dec, s, tc, dur, gaia)
        except Exception as e:
            row, pan = dict(sector=s, tc=tc, note=str(e)), None
        row["label"] = lab
        rows.append(row)
        pans.append((lab, s, pan))
        print(f"[cutout] S{s} {tc:.3f} {lab}: {row.get('note', 'ok')}")
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "tic230133616_pixel_check.csv", index=False)
    cols = ["label", "sector", "tc", "ap_depth_ppm", "ap_depth_err_ppm", "diff_centroid_offset_arcsec",
            "diff_peak_col", "diff_peak_row", "target_col", "target_row",
            "incen_minus_outcen_col_mpix", "incen_minus_outcen_row_mpix"]
    print(res[[c for c in cols if c in res]].to_string(index=False))
    plot(pans, dur)


def plot(pans, dur):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(pans)
    fig, axs = plt.subplots(n, 3, figsize=(12, 3 * n), squeeze=False)
    for (lab, s, p), ax in zip(pans, axs):
        if p is None:
            ax[0].set_title(f"S{s} {lab}: no cutout")
            for a in ax:
                a.axis("off")
            continue
        for a, im, ttl in ((ax[0], p["img_out"], "out of transit"), (ax[1], p["diff"], "out - in")):
            a.imshow(im, origin="lower", cmap="viridis")
            a.contour(p["ap"], levels=[0.5], colors="w", linewidths=0.8)
            a.plot(p["x0"], p["y0"], "r+", ms=12)
            a.scatter(p["gx"], p["gy"], s=12, facecolors="none", edgecolors="orange")
            a.set_xlim(-0.5, im.shape[1] - 0.5)
            a.set_ylim(-0.5, im.shape[0] - 0.5)
            a.set_title(f"S{s} {lab}\n{ttl}", fontsize=8)
        m = np.abs(p["t"] - p["tc"]) < 1.0
        x = (p["t"][m] - p["tc"]) * 24
        ax[2].plot(x, (p["f"][m] - 1) * 1e6, ".", ms=2, color="0.6")
        e = np.arange(-24, 24.5, 0.5)
        idx = np.digitize(x, e)
        ax[2].plot([x[idx == j].mean() for j in np.unique(idx)],
                   [(p["f"][m][idx == j].mean() - 1) * 1e6 for j in np.unique(idx)], "o", ms=3, color="C0")
        ax[2].axvspan(-dur * 12, dur * 12, color="C3", alpha=0.1)
        ax[2].set_title("TESScut aperture flux (ppm)", fontsize=8)
        ax[2].set_xlabel("hours")
    fig.tight_layout()
    fig.savefig(OUT / "tic230133616_pixel_check.png", dpi=90)
    plt.close(fig)


if __name__ == "__main__":
    main()
