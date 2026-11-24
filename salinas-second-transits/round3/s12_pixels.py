#!/usr/bin/env python3
"""Round 3, section 12: pixel-level and Gaia vetting of 275647095, 241550976 and 107113048.

The same check round 2 ran on TIC 230133616 (round2/s5_tic230133616_pixels.py), on every clean transit
from round 3 (results/round3/tic*_transits.csv, mid-times from the fitted ephemeris):

- Gaia DR3 sources within 2' brighter than T + 6, the depth each would need to make the observed dip,
  and the target's own binarity indicators (RUWE, ipd_frac_multi_peak, non_single_star, RV error and
  amplitude).
- Per transit, a 15x15 TESScut cutout: out-of-transit image, difference image (out minus in), its
  flux-weighted centroid within 3 px of the target, the brightest difference pixel anywhere in the
  cutout, and the aperture's in- minus out-of-transit centroid shift.
- Stacked: each transit's difference-centroid offset converted to sky (arcsec east, north of the
  target) and averaged, with the error from the scatter between transits. Each neighbour that could
  make the dip (needs under 100% of its own flux) is compared with that mean offset, in sigma. A transit
  whose aperture centroid shift is over 3x the star's median shift (a pointing or scattered-light
  excursion: its difference image is dominated by a neighbour) is flagged, and the stack is given with
  and without flagged transits.
- Gaia DR3 epoch RV scatter. For stars this faint DR3 gives the median of the epoch RVs, with error
  sqrt(pi / 2N) times their scatter (Katz et al. 2023), so that scatter is radial_velocity_error *
  sqrt(2N / pi). It includes the noise, so it is only an upper bound on any orbital signal.

  python round3/s12_pixels.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "round2"))
import tesscut_phot as tp  # noqa: E402
from s5_tic230133616_pixels import epoch_cutout, gaia_cone, tmag_from_gaia  # noqa: E402

WORK, OUT = HERE / "work", HERE / "results" / "round3"
STARS = [275647095, 241550976, 107113048]
GAIA_EXTRA = ["ipd_frac_multi_peak", "non_single_star", "radial_velocity", "radial_velocity_error",
              "rv_amplitude_robust", "rv_nb_transits"]


def target_binarity(ra, dec):
    """Binarity columns for the Gaia source nearest the target (ESA archive only)."""
    from astroquery.gaia import Gaia
    cols = ", ".join(["source_id", "phot_g_mean_mag", "ruwe"] + GAIA_EXTRA)
    q = f"""SELECT {cols}, DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec})) * 3600 AS sep_arcsec
            FROM gaiadr3.gaia_source
            WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {3 / 3600}))
            ORDER BY phot_g_mean_mag"""
    try:
        r = Gaia.launch_job(q).get_results().to_pandas()
        return r.iloc[0].to_dict() if len(r) else {}
    except Exception as e:
        return dict(error=str(e)[:200])


def sky_offset(tpf, col, row, ra0, dec0):
    """Pixel position -> arcsec east and north of (ra0, dec0)."""
    ra, dec = tpf.wcs.all_pix2world([[col, row]], 0)[0]
    return (ra - ra0) * np.cos(np.radians(dec0)) * 3600, (dec - dec0) * 3600


def one_star(tic, cov):
    c = cov.loc[tic]
    tr = pd.read_csv(OUT / f"tic{tic}_transits.csv")
    summ = pd.read_csv(OUT / "oddeven_summary.csv").set_index("tic").loc[tic]
    P, T0, dur, depth = summ.P_d, summ.T0_btjd, summ.t14_h / 24, summ.depth_ppm * 1e-6
    clean = tr[tr["clean"]].copy()
    clean["tc_eph"] = T0 + clean["epoch"] * P

    gaia = gaia_cone(c.ra, c.dec)
    gaia["tmag_est"] = tmag_from_gaia(gaia["phot_g_mean_mag"].to_numpy(float), gaia["bp_rp"].to_numpy(float))
    gaia["dT"] = gaia["tmag_est"] - c.tmag
    gaia = gaia[gaia["tmag_est"] < c.tmag + 6].reset_index(drop=True)
    # the dip in the blended light is depth; a source fainter by dT needs depth * (1 + 10^(0.4 dT)) of its own
    # flux if it were the only other light, so depth * 10^(0.4 dT) is a lower bound
    gaia["needed_own_depth"] = depth * 10 ** (0.4 * gaia["dT"].clip(lower=0))
    gaia["east_arcsec"] = (gaia["ra"] - c.ra) * np.cos(np.radians(c.dec)) * 3600
    gaia["north_arcsec"] = (gaia["dec"] - c.dec) * 3600
    is_target = gaia["sep_arcsec"] < 3
    gaia.insert(0, "tic", tic)
    gaia["is_target"] = is_target

    rows, pans = [], []
    for r in clean.itertuples():
        try:
            row, pan = epoch_cutout(c.ra, c.dec, int(r.sector), r.tc_eph, dur, gaia)
            if pan is not None and np.isfinite(row.get("diff_centroid_col", np.nan)):
                tpf = tp.cutout(c.ra, c.dec, int(r.sector), WORK / "tesscut")
                e, n = sky_offset(tpf, row["diff_centroid_col"], row["diff_centroid_row"], c.ra, c.dec)
                row.update(diff_centroid_east_arcsec=e, diff_centroid_north_arcsec=n)
                d = pan["diff"]
                pk = np.unravel_index(np.nanargmax(d), d.shape)
                row.update(diff_global_peak_col=int(pk[1]), diff_global_peak_row=int(pk[0]),
                           diff_global_peak_dist_px=float(np.hypot(pk[1] - row["target_col"], pk[0] - row["target_row"])))
        except Exception as ex:
            row, pan = dict(sector=int(r.sector), tc=r.tc_eph, note=str(ex)[:200]), None
        row.update(tic=tic, epoch=int(r.epoch), lc_source=r.source)
        rows.append(row)
        pans.append((f"E{int(r.epoch)}", int(r.sector), pan))
        print(f"  TIC {tic} E{int(r.epoch)} S{int(r.sector)}: {row.get('note', 'ok')}", flush=True)
    res = pd.DataFrame(rows)
    shift = np.hypot(res["incen_minus_outcen_col_mpix"], res["incen_minus_outcen_row_mpix"])
    res["excursion"] = shift > 3 * np.nanmedian(shift)

    def stack(ok):
        n = len(ok)
        me, mn = ok["diff_centroid_east_arcsec"].mean(), ok["diff_centroid_north_arcsec"].mean()
        if n >= 2:
            se_e = ok["diff_centroid_east_arcsec"].std(ddof=1) / np.sqrt(n)
            se_n = ok["diff_centroid_north_arcsec"].std(ddof=1) / np.sqrt(n)
        else:
            se_e = se_n = np.nan
        return n, me, mn, float(np.hypot(se_e, se_n) / np.sqrt(2))  # per-axis error, averaged

    good = res.dropna(subset=["diff_centroid_east_arcsec"])
    n_all, e_all, n_all_n, se_all = stack(good)
    n, mean_e, mean_n, se = stack(good[~good["excursion"]])
    off = float(np.hypot(mean_e, mean_n))
    cand = gaia[~gaia["is_target"] & (gaia["needed_own_depth"] < 1)].copy()
    cand["sep_from_centroid_arcsec"] = np.hypot(cand["east_arcsec"] - mean_e, cand["north_arcsec"] - mean_n)
    cand["sigma_from_centroid"] = cand["sep_from_centroid_arcsec"] / se if se > 0 else np.nan
    gaia = gaia.merge(cand[["source_id", "sep_from_centroid_arcsec", "sigma_from_centroid"]], on="source_id", how="left")
    binar = target_binarity(c.ra, c.dec)
    nrv = binar.get("rv_nb_transits", np.nan)
    rv_scatter = binar.get("radial_velocity_error", np.nan) * np.sqrt(2 * nrv / np.pi)
    summary = dict(tic=tic, tmag=c.tmag, n_transits_imaged=len(good),
                   excursion_epochs=" ".join(str(int(x)) for x in res.loc[res["excursion"], "epoch"]),
                   n_transits_stacked=n,
                   diff_centroid_east_arcsec=mean_e, diff_centroid_north_arcsec=mean_n,
                   diff_centroid_offset_arcsec=off, diff_centroid_err_arcsec=se,
                   offset_sigma=off / se if se > 0 else np.nan,
                   all_transits_offset_arcsec=float(np.hypot(e_all, n_all_n)), all_transits_err_arcsec=se_all,
                   max_diff_global_peak_dist_px=float(res.get("diff_global_peak_dist_px", pd.Series([np.nan])).max()),
                   max_abs_centroid_shift_mpix=float(np.nanmax(np.abs(res[["incen_minus_outcen_col_mpix",
                                                                           "incen_minus_outcen_row_mpix"]].to_numpy())))
                   if "incen_minus_outcen_col_mpix" in res else np.nan,
                   ap_depth_ppm_median=float(res["ap_depth_ppm"].median()) if "ap_depth_ppm" in res else np.nan,
                   n_gaia_within_2arcmin=int(len(gaia)),
                   n_neighbours_able=int(len(cand)),
                   nearest_able_neighbour_sep_arcsec=float(cand["sep_arcsec"].min()) if len(cand) else np.nan,
                   min_sigma_able_neighbour=float(cand["sigma_from_centroid"].min()) if len(cand) else np.nan,
                   target_ruwe=binar.get("ruwe", np.nan),
                   **{f"target_{k}": binar.get(k, np.nan) for k in GAIA_EXTRA},
                   target_rv_epoch_scatter_kms=rv_scatter,
                   gaia_error=binar.get("error", ""))
    res.to_csv(OUT / f"tic{tic}_pixel_check.csv", index=False)
    gaia.to_csv(OUT / f"tic{tic}_gaia_neighbours.csv", index=False)
    plot(tic, pans, dur, gaia, mean_e, mean_n, se)
    return summary


def plot(tic, pans, dur, gaia, me, mn, se):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(pans)
    fig, axs = plt.subplots(n + 1, 3, figsize=(12, 3 * (n + 1)), squeeze=False)
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
            a.set_title(f"S{s} {lab}: {ttl}", fontsize=8)
        m = np.abs(p["t"] - p["tc"]) < 1.0
        x = (p["t"][m] - p["tc"]) * 24
        ax[2].plot(x, (p["f"][m] - 1) * 1e6, ".", ms=2, color="0.6")
        e = np.arange(-24, 24.5, 0.5)
        idx = np.digitize(x, e)
        u = np.unique(idx)
        ax[2].plot([x[idx == j].mean() for j in u], [(p["f"][m][idx == j].mean() - 1) * 1e6 for j in u],
                   "o", ms=3, color="C0")
        ax[2].axvspan(-dur * 12, dur * 12, color="C3", alpha=0.1)
        ax[2].set_title("TESScut aperture flux (ppm)", fontsize=8)
        ax[2].set_xlabel("hours")
    a = axs[-1, 0]
    g = gaia
    a.scatter(g["east_arcsec"], g["north_arcsec"], s=np.clip(80 - 10 * g["dT"], 5, 120), c="orange", alpha=0.7)
    able = g["needed_own_depth"] < 1
    a.scatter(g.loc[able & ~g["is_target"], "east_arcsec"], g.loc[able & ~g["is_target"], "north_arcsec"],
              s=60, facecolors="none", edgecolors="C3", label="could make the dip")
    a.plot(0, 0, "r+", ms=14, label="target")
    if np.isfinite(se):
        a.errorbar(me, mn, se, se, fmt="kx", ms=8, label="mean diff centroid")
    a.set(xlabel="arcsec east", ylabel="arcsec north", title=f"TIC {tic}: Gaia DR3 within 2' (T < T*+6)",
          xlim=(130, -130), ylim=(-130, 130), aspect="equal")
    a.legend(fontsize=7)
    for a in axs[-1, 1:]:
        a.axis("off")
    fig.tight_layout()
    fig.savefig(OUT / f"tic{tic}_pixel_check.png", dpi=80)
    plt.close(fig)


def main():
    cov = pd.read_csv(WORK / "coverage.csv").set_index("tic")
    rows = []
    for tic in STARS:
        print(f"TIC {tic}", flush=True)
        rows.append(one_star(tic, cov))
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "pixel_summary.csv", index=False)
    print(s.T.to_string())


if __name__ == "__main__":
    main()
