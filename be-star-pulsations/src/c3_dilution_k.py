"""Check 3: measure the dilution factor k from SAP flux.

Stage 4's dilution correction assumes ln(F_TESS) rises by 0.921*k per mag of W1 excess, i.e.
Tmag falls by k per mag: slope dTmag/d(disk_excess) = -k. A positive slope means the disk
dims the optical (absorption, e.g. shell stars).

For stars whose preferred product is SPOC120 or TESS-SPOC, main crowd-cut WISE-usable
sectors of that product only, stars with >= 3 such sectors: per-sector median SAP_FLUX
(QUALITY == 0), corrected as SAP * CROWDSAP / FLFRCSAP, Tmag = -2.5 log10(F) + 20.44.
Within-star pooled slope of Tmag on disk_excess with a bootstrap over stars; per-star OLS
slopes; stars with slope > 0 and two-sided p < 0.05 are shell candidates. Then the stage 4
main test (A1, B1, flicker_rms; k = 0, both nulls) without the candidates.
Outputs: results/checks/dilution_k.csv (summary), dilution_k_per_star.csv,
shell_candidates.csv. Run with BE_RUN=full.
"""
import numpy as np, pandas as pd
from astropy.io import fits
from scipy import stats
from common import RESULTS, TESS_DIR, SEED
import s4_test as T, stats4, c1_phase_terms as C1

OUT = RESULTS / "checks"
ZP = 20.44


def sap_mag(path):
    with fits.open(TESS_DIR / path, memmap=False) as h:
        d, hd = h[1].data, h[1].header
        f = np.asarray(d["SAP_FLUX"], float)
        q = np.asarray(d["QUALITY"])
        ok = np.isfinite(f) & (q == 0) & (f > 0)
        crowd, flfr = hd.get("CROWDSAP", np.nan), hd.get("FLFRCSAP", np.nan)
    if ok.sum() < 100 or not (np.isfinite(crowd) and np.isfinite(flfr) and flfr > 0):
        return np.nan, crowd, flfr
    F = np.median(f[ok]) * crowd / flfr
    return -2.5 * np.log10(F) + ZP, crowd, flfr


def main():
    sub, vsub, disks, dos = C1.main_subset()
    m = (sub.wise_usable.astype(bool) & sub.preferred.isin(["SPOC120", "TESS-SPOC"])
         & (sub["product"] == sub.preferred)).values
    d = sub[m].copy()
    d = d[d.groupby("tic").tic.transform("size") >= 3]
    res = [sap_mag(p) for p in d.file]
    d["tmag_sap"] = [r[0] for r in res]
    d["crowdsap_hdr"] = [r[1] for r in res]
    d["flfrcsap_hdr"] = [r[2] for r in res]
    d = d[np.isfinite(d.tmag_sap)]
    d = d[d.groupby("tic").tic.transform("size") >= 3]
    # pooled within-star slope + bootstrap over stars
    xd = d.disk_excess - d.groupby("tic").disk_excess.transform("mean")
    yd = d.tmag_sap - d.groupby("tic").tmag_sap.transform("mean")
    g = pd.DataFrame(dict(tic=d.tic, sxy=xd * yd, sxx=xd * xd)).groupby("tic").sum()
    slope = g.sxy.sum() / g.sxx.sum()
    rng = np.random.default_rng(SEED)
    W = rng.multinomial(len(g), np.full(len(g), 1 / len(g)), size=5000)
    bs = (W @ g.sxy.values) / (W @ g.sxx.values)
    lo68, hi68 = np.percentile(bs, [16, 84])
    lo95, hi95 = np.percentile(bs, [2.5, 97.5])
    # noise: Tmag differences between sector pairs of one star with |dx| < 0.03 mag
    pair = []
    for _, s in d.groupby("tic"):
        x, y = s.disk_excess.values, s.tmag_sap.values
        for a in range(len(s)):
            for b in range(a + 1, len(s)):
                if abs(x[a] - x[b]) < 0.03:
                    pair.append(y[a] - y[b])
    pair = np.array(pair)
    rows = []
    for tic, s in d.groupby("tic"):
        x, y = s.disk_excess.values, s.tmag_sap.values
        n = len(s)
        if np.ptp(x) <= 0:
            continue
        lr = stats.linregress(x, y)
        rows.append(dict(tic=tic, n_sec=n, product=s.preferred.iloc[0], x_range=np.ptp(x),
                         tmag_range=np.ptp(y), slope=lr.slope, slope_se=lr.stderr,
                         p=lr.pvalue, tmag_median=np.median(y)))
    ps = pd.DataFrame(rows)
    ps["shell_candidate"] = (ps.slope > 0) & (ps.p < 0.05)
    ps.to_csv(OUT / "dilution_k_per_star.csv", index=False)
    shell = ps[ps.shell_candidate].sort_values("p")
    shell.to_csv(OUT / "shell_candidates.csv", index=False)
    q = np.percentile(ps.slope, [10, 25, 50, 75, 90])
    summ = dict(n_stars=len(g), n_sectors=len(d), pooled_slope=slope, slope_lo68=lo68,
                slope_hi68=hi68, slope_lo95=lo95, slope_hi95=hi95, k=-slope, k_lo95=-hi95,
                k_hi95=-lo95, per_star_p10=q[0], per_star_p25=q[1], per_star_median=q[2],
                per_star_p75=q[3], per_star_p90=q[4], frac_positive=float(np.mean(ps.slope > 0)),
                n_sig_neg=int(((ps.slope < 0) & (ps.p < 0.05)).sum()),
                n_shell_candidates=len(shell), pair_noise_n=len(pair),
                pair_noise_rms_mag=float(np.sqrt(np.mean(pair ** 2))) if len(pair) else np.nan,
                pair_noise_mad_mag=float(1.4826 * np.median(np.abs(pair))) if len(pair) else np.nan,
                within_star_x_sd=float(np.sqrt(np.mean(xd ** 2))))
    table = stats4.build_swap_table(sub.mid_mjd.values, disks, T.evaluate)
    # the same pooled slope through stats4 (within-star permutation and cross-star swap
    # nulls), so a time-dependent SAP systematic shared across stars cannot fake it
    ty = np.full(len(sub), np.nan)
    ty[d.index.values] = d.tmag_sap.values
    r = stats4.run_test(sub.tic.values, sub.disk_excess.values, ty,
                        sub.phase_label.values.astype(object), sub.wise_usable.values.astype(bool),
                        sub.mid_mjd.values, disks, dos, T.evaluate, n_null=T.N_NULL,
                        n_boot=T.N_BOOT, seed=SEED, swap_table=table)
    summ.update(stats4_slope=r["beta"], stats4_n_stars=r["n_stars"], p_perm=r["p_perm"],
                p_swap=r["p_swap"], swap_null_mean=r["swap_null_mean"],
                n_candidate_sectors=int(m.sum()))
    print(summ, flush=True)
    # stage 4 main test without the candidates
    main = vsub["main (crowd-cut)"] & ~sub.tic.isin(shell.tic).values
    for metric in ["A1", "B1", "flicker_rms"]:
        y = np.where(main, np.log(sub[metric].values), np.nan)
        r = stats4.run_test(sub.tic.values, sub.disk_excess.values, y,
                            sub.phase_label.values.astype(object),
                            sub.wise_usable.values.astype(bool), sub.mid_mjd.values, disks, dos,
                            T.evaluate, n_null=T.N_NULL, n_boot=T.N_BOOT, seed=SEED,
                            swap_table=table)
        for kk in ("n_stars", "beta", "beta_lo", "beta_hi", "p_perm", "p_swap"):
            summ[f"noshell_{metric}_{kk}"] = r[kk]
        print(metric, {kk: r[kk] for kk in ("n_stars", "beta", "beta_lo", "beta_hi", "p_perm",
                                             "p_swap")}, flush=True)
    pd.DataFrame([summ]).T.rename(columns={0: "value"}).to_csv(OUT / "dilution_k.csv",
                                                               index_label="quantity")


if __name__ == "__main__":
    main()
