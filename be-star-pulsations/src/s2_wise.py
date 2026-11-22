"""Stage 2: WISE disk phase at each TESS sector for the pilot stars.

W1 curve = China-VO per-visit curve (same processing as Jia et al.), extended past its last
visit with IRSA NEOWISE-R single exposures (cuts as Jian et al. 2024, binned by visit).
IRSA is queried once per 25 stars with a TAP upload join (3 arcsec).
"""
import sys
import numpy as np, pandas as pd
import pyvo
from astropy.io import fits
from astropy.table import Table
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import splev
from common import CAT_DIR, WISE_DIR, TESS_DIR, PILOT, RESULTS, FIGURES, BTJD_TO_MJD, figpath
from wise import irsa_cuts, bin_visits, DiskCurve

TAP = "https://irsa.ipac.caltech.edu/TAP"
Q = """SELECT u.tic, n.mjd, n.w1mpro, n.w1sigmpro, n.w2mpro, n.w2sigmpro, n.qual_frame,
n.cc_flags, n.nb, n.saa_sep, n.moon_masked
FROM neowiser_p1bs_psd AS n, TAP_UPLOAD.pos AS u
WHERE CONTAINS(POINT('ICRS', n.ra, n.dec), CIRCLE('ICRS', u.ra0, u.dec0, 0.000833333)) = 1"""
VIS_DIR = WISE_DIR / "visits"
VIS_DIR.mkdir(exist_ok=True)


def fetch_irsa(stars, chunk=25):
    path = WISE_DIR / "neowise_irsa.parquet"
    have = pd.read_parquet(path) if path.exists() else pd.DataFrame(columns=["tic"])
    done = set(pd.read_csv(WISE_DIR / "irsa_done.csv").tic) if (WISE_DIR / "irsa_done.csv").exists() else set()
    todo = stars[~stars.tic.isin(done)]
    svc = pyvo.dal.TAPService(TAP)
    import time
    for i in range(0, len(todo), chunk):
        t0 = time.time()
        ch = todo.iloc[i:i + chunk]
        up = Table({"tic": ch.tic.values.astype("int64"), "ra0": ch.ra.values, "dec0": ch.dec.values})
        r = svc.run_async(Q, uploads={"pos": up}).to_table().to_pandas()
        for c in ("cc_flags", "moon_masked"):
            r[c] = r[c].astype(str)
        have = pd.concat([have, r], ignore_index=True)
        have.to_parquet(path, index=False)
        done |= set(ch.tic)
        pd.DataFrame({"tic": sorted(done)}).to_csv(WISE_DIR / "irsa_done.csv", index=False)
        print(f"IRSA chunk {i//chunk+1}/{-(-len(todo)//chunk)}: {len(ch)} stars, {len(r)} rows, "
              f"{time.time()-t0:.0f} s", flush=True)
    return have


def sector_mid(path):
    with fits.open(TESS_DIR / path, memmap=False) as h:
        t = np.asarray(h[1].data["TIME"], float)
    t = t[np.isfinite(t)]
    return 0.5 * (t.min() + t.max())


def build_visits(tic, gaia, irsa):
    f = CAT_DIR / "wise_lc" / "wise_lc" / f"{gaia}.csv"
    cvo = pd.read_csv(f) if f.exists() else None
    ir = irsa[irsa.tic == tic]
    ir = irsa_cuts(ir) if len(ir) else ir
    iv = bin_visits(ir.mjd.values, ir.w1mpro.values, ir.w2mpro.values) if len(ir) else pd.DataFrame()
    info = dict(cvo_visits=0 if cvo is None else len(cvo), irsa_visits=len(iv), ext_visits=0,
                overlap_n=0, overlap_dw1=np.nan)
    if cvo is None:
        v = iv.assign(source="irsa")
    else:
        v = cvo.assign(source="chinavo")
        if len(iv):
            # agreement where both have the same visit
            j = np.abs(iv.mjd.values[:, None] - cvo.mjd.values[None, :])
            near = j.min(axis=1) < 5
            if near.any():
                k = j.argmin(axis=1)[near]
                info["overlap_n"] = int(near.sum())
                info["overlap_dw1"] = float(np.median(iv.w1mpro_ep.values[near] - cvo.w1mpro_ep.values[k]))
            ext = iv[iv.mjd > cvo.mjd.max() + 10].assign(source="irsa")
            info["ext_visits"] = len(ext)
            v = pd.concat([v, ext], ignore_index=True)
    v = v.sort_values("mjd").reset_index(drop=True)
    v.to_csv(VIS_DIR / f"{tic}.csv", index=False)
    return v, info


def main():
    cat = pd.read_parquet(RESULTS / "catalogue.parquet").set_index("tic")
    samp = pd.read_csv(PILOT / "sample.csv")
    samp = samp[samp.status == "pass"]
    sec = pd.read_csv(PILOT / "sectors.csv")
    stars = cat.loc[samp.tic, ["ra", "dec", "gaia_dr3_id"]].reset_index()
    irsa = fetch_irsa(stars)
    rows, sinfo, disks = [], [], {}
    for _, s in stars.iterrows():
        v, info = build_visits(int(s.tic), int(s.gaia_dr3_id), irsa)
        w1med = float(np.median(v.w1mpro_ep))
        d = DiskCurve(v.mjd.values, v.w1mpro_ep.values, v.w1_w2.values)
        disks[int(s.tic)] = d
        info.update(tic=int(s.tic), n_visits=len(v), wise_first_mjd=v.mjd.min(),
                    wise_last_mjd=v.mjd.max(), w1_median=w1med, bright_w1lt8=w1med < 8.0,
                    w1_quiet=d.w1_quiet, w1_spline_range=d.w1_quiet - float(np.min(splev(
                        np.linspace(d.t0, d.t1, 4000), d.tck))))
        sinfo.append(info)
        ss = sec[sec.tic == s.tic].copy()
        ss["mid_btjd"] = [sector_mid(p) for p in ss.file]
        ss["mid_mjd"] = ss.mid_btjd + BTJD_TO_MJD
        o = d.at(ss.mid_mjd.values)
        for k in ("w1_fit", "disk_excess", "dw1dt", "w1w2", "wise_gap", "phase_label", "wise_usable"):
            ss[k] = o[k]
        ss["bright_w1lt8"] = w1med < 8.0
        rows.append(ss)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(PILOT / "sector_wise.csv", index=False)
    si = pd.DataFrame(sinfo)
    si.to_csv(PILOT / "star_wise.csv", index=False)
    print("stars", len(si), "sectors", len(out), "WISE-usable", int(out.wise_usable.sum()))
    print("last NEOWISE visit over pilot (MJD):", si.wise_last_mjd.max())
    print("China-VO vs IRSA overlap dW1 median of star medians:", si.overlap_dw1.median(),
          "n stars", si.overlap_dw1.notna().sum(), "| ext visits per star median", si.ext_visits.median())
    print("W1<8 stars:", int(si.bright_w1lt8.sum()))
    u = out[out.wise_usable]
    print("phase labels (usable):", u.phase_label.value_counts().to_dict())
    print("unusable reasons: outside baseline",
          int(((out.mid_mjd < out.tic.map(si.set_index('tic').wise_first_mjd)) |
               (out.mid_mjd > out.tic.map(si.set_index('tic').wise_last_mjd))).sum()),
          "gap>250", int(((out.wise_gap > 250) & np.isfinite(out.wise_gap)).sum()))
    spline_check(stars.tic.values[:10], disks, out)


def spline_check(tics, disks, sec):
    fig, axs = plt.subplots(5, 2, figsize=(12, 13), sharex=True)
    for ax, tic in zip(axs.flat, tics):
        d = disks[int(tic)]
        v = pd.read_csv(VIS_DIR / f"{int(tic)}.csv")
        g = np.linspace(d.t0, d.t1, 1500)
        for src, c in (("chinavo", "C0"), ("irsa", "C1")):
            m = v.source == src
            ax.errorbar(v.mjd[m], v.w1mpro_ep[m], v.w1sigmpro_ep[m], fmt="o", ms=3, c=c,
                        label=f"{src} visits")
        ax.plot(g, splev(g, d.tck), "k-", lw=1, label="spline s=0.04")
        alt = DiskCurve(v.mjd.values, v.w1mpro_ep.values, s=0.04 * len(v) / 25)
        ax.plot(g, splev(g, alt.tck), "g--", lw=1, label="s=0.04 N/25")
        for t in sec[sec.tic == tic].mid_mjd:
            ax.axvline(t, c="r", alpha=0.4, lw=1)
        ax.invert_yaxis()
        ax.set_title(f"TIC {int(tic)}  (N={len(v)})", fontsize=9)
        ax.set_ylabel("W1")
    axs.flat[0].legend(fontsize=7)
    for ax in axs[-1]:
        ax.set_xlabel("MJD (red: TESS sector mid-times)")
    fig.tight_layout()
    fig.savefig(figpath("spline_check"), dpi=110)


if __name__ == "__main__":
    main()
