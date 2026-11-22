"""Stage 4: within-star test of ln(metric) against WISE disk excess."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import PILOT, FIGURES, WISE_DIR, SEED, figpath
from wise import DiskCurve
import stats4

METRICS = ["A1", "B1", "B2", "B3", "red_noise", "flicker_rms", "white_noise"]  # last = data-quality control
KS = [0.0, 0.25, 0.5]
N_NULL, N_BOOT = 10000, 2000


def evaluate(d, t):
    o = d.at(t)
    return o["disk_excess"], o["phase_label"].astype(object), o["wise_usable"]


def load():
    sm = pd.read_csv(PILOT / "sector_metrics.csv")
    sw = pd.read_csv(PILOT / "star_wise.csv").set_index("tic")
    samp = pd.read_csv(PILOT / "sample.csv").set_index("tic")
    tics = sorted(sm.tic.unique())
    disks, dos = [], {}
    for j, tic in enumerate(tics):
        v = pd.read_csv(WISE_DIR / "visits" / f"{tic}.csv")
        disks.append(DiskCurve(v.mjd.values, v.w1mpro_ep.values, v.w1_w2.values))
        dos[tic] = j
    sm["preferred"] = sm.tic.map(samp.preferred)
    return sm, sw, disks, dos


def variants(sm):
    base = sm.metrics_ok.values.astype(bool)
    crowd = sm.crowd_ok.fillna(False).values.astype(bool)
    spoc_star = sm.preferred.isin(["SPOC120", "TESS-SPOC"]).values
    qlp_star = (sm.preferred == "QLP").values
    main = base & crowd
    return {
        "main (crowd-cut)": main,
        "no crowding cut": base,
        "drop mixed_pipeline": main & ~sm.mixed_pipeline.values.astype(bool),
        "drop W1<8 stars": main & ~sm.bright_w1lt8.values.astype(bool),
        "SPOC/TESS-SPOC stars only": main & spoc_star,
        "QLP stars only": main & qlp_star,
    }


def noise_floor(sm, sel, metric):
    d = sm[sel & sm.wise_usable.values.astype(bool)]
    dy = []
    for _, g in d.groupby("tic"):
        x, y = g.disk_excess.values, np.log(g[metric].values)
        for a in range(len(g)):
            for b in range(a + 1, len(g)):
                if abs(x[a] - x[b]) < 0.03:
                    dy.append(y[a] - y[b])
    dy = np.array(dy)
    if len(dy) < 2:
        return len(dy), np.nan, np.nan
    # pairs differences: sign is arbitrary, so the scatter is the rms of dy
    return len(dy), float(np.sqrt(np.mean(dy ** 2))), float(1.4826 * np.median(np.abs(dy)))


def main():
    sm, sw, disks, dos = load()
    sid = sm.tic.values
    t = sm.mid_mjd.values
    x = sm.disk_excess.values
    lab = sm.phase_label.values.astype(object)
    wu = sm.wise_usable.values.astype(bool)
    rows, keep = [], {}
    table = stats4.build_swap_table(t, disks, evaluate)   # every sector under every disk
    for vname, sel in variants(sm).items():
        for metric in METRICS:
            y0 = np.where(sel, np.log(sm[metric].values), np.nan)
            nf = noise_floor(sm, sel, metric)
            for k in KS:
                y = y0 + 0.921 * k * x
                r = stats4.run_test(sid, x, y, lab, wu, t, disks, dos, evaluate,
                                    n_null=N_NULL, n_boot=N_BOOT, seed=SEED, swap_table=table,
                                    keep_nulls=(vname.startswith("main") and metric == "A1" and k == 0))
                if "_perm_null" in r:
                    keep = r
                r = {kk: vv for kk, vv in r.items() if not kk.startswith("_")}
                r.update(variant=vname, metric=metric, k=k, noise_pairs=nf[0],
                         noise_rms_dlny=nf[1], noise_mad_dlny=nf[2])
                rows.append(r)
                print(vname, metric, k, {kk: (round(vv, 4) if isinstance(vv, float) else vv)
                                         for kk, vv in r.items() if kk in
                                         ("n_stars", "n_sectors", "beta", "p_perm", "p_swap")}, flush=True)
    cols = ["variant", "metric", "k", "n_stars", "n_sectors", "beta", "beta_lo", "beta_hi", "rho",
            "contrast", "n_build", "n_diss", "p_perm", "p_swap", "p_perm_contrast",
            "p_swap_contrast", "swap_null_mean", "noise_pairs", "noise_rms_dlny", "noise_mad_dlny"]
    out = pd.DataFrame(rows)[cols]
    out.to_csv(PILOT / "test_results.csv", index=False)
    gonogo(sm)
    pooled_fig(keep)


def gonogo(sm):
    u = sm[sm.usable.astype(bool)]
    uc = u[u.crowd_ok.fillna(False).astype(bool)]
    def count(d):
        g = d.groupby("tic").disk_excess
        n2 = g.count() >= 2
        rng = (g.max() - g.min())
        return int(n2.sum()), int((n2 & (rng >= 0.1)).sum())
    print("GO/NO-GO stars with >=2 usable sectors, and of those with disk range >= 0.1 mag:")
    print("  usable (WISE+metrics):", count(u), "  usable & crowd_ok:", count(uc))


def pooled_fig(r):
    if not r:
        return
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.scatter(r["_xd"], r["_yd"], s=12, alpha=0.7)
    g = np.linspace(r["_xd"].min(), r["_xd"].max(), 10)
    ax.plot(g, r["beta"] * g, "r-", label=f"beta = {r['beta']:.2f} "
            f"[{r['beta_lo']:.2f}, {r['beta_hi']:.2f}] per mag")
    ax.set_xlabel("disk_excess - star mean (mag)")
    ax.set_ylabel("ln A1 - star mean")
    ax.set_title(f"Pilot, main sample: {r['n_stars']} stars, {r['n_sectors']} sectors")
    ax.legend(loc="upper left", fontsize=9)
    ins = ax.inset_axes([0.62, 0.08, 0.35, 0.3])
    ins.hist(r["_perm_null"], bins=60, histtype="step", density=True, label="within-star")
    ins.hist(r["_swap_null"][np.isfinite(r["_swap_null"])], bins=60, histtype="step",
             density=True, label="cross-star")
    ins.axvline(r["beta"], c="r")
    ins.set_title(f"null beta; p = {r['p_perm']:.3f}, {r['p_swap']:.3f}", fontsize=7)
    ins.tick_params(labelsize=6)
    ins.legend(fontsize=5)
    fig.tight_layout()
    fig.savefig(figpath("pooled"), dpi=120)


if __name__ == "__main__":
    main()
