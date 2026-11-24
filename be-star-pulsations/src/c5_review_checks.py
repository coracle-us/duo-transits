"""Reruns asked for by the Mythos review of the RNAAS draft (rnaas/REVIEW_mythos.md, 2026-11-23).

Every row in review_checks.csv goes through stats4.run_test (within-star permutation and
cross-star swap nulls, 10,000 draws; bootstrap over stars, 2,000), main sample, k = 0 unless
stated. Rows:
  M1  confirmed subset: stars Gaia DR3 ESP-ELS labels 'beStar' (H-alpha emission in the RP
      spectrum; cached gaia_dr3_params.parquet), and main minus ESP-ELS Herbig/WR stars.
      (CDS was down, so Jian et al. 2024 could not be cross-matched.)
  M2  stable-f1 stars (every main sector's f1 within 10% of the star's median f1) vs the rest.
  M3  SPOC/TESS-SPOC stars with the measured dilution correction, k = 0.213 (results/checks/
      dilution_k.csv). QLP is not diluted (Huang et al. 2020: difference flux + TIC flux).
  M4  terciles of the star-median white_noise level: white_noise, A1, B1.
  M5  drop the 10, 50, 100 stars with the largest within-star sum of squared excess.
  S3  star x cadence groups (only changes within one star and one FFI cadence).
  S1  relative-frequency bins recentred so one bin spans f/f1 = 0.891-1.122 (main, SPOC,
      stable-f1).
Bootstrap-only extras (no nulls) go to review_extras.csv: two-way star+sector fixed effects,
ln white_noise on excess + ln B1, QLP - SPOC beta gaps. review_signs.csv: per-star slope sign
by the star's excess range.
"""
import numpy as np, pandas as pd
from common import PILOT, CACHE, RESULTS, SEED
from metrics import LOG_GRID
import s4_test as T, stats4
import s6_freq_beta as F

OUT = RESULTS / "checks"
N_NULL, N_BOOT = 10000, 2000
K_MEAS = 0.213
METRICS = ["A1", "B1", "red_noise", "flicker_rms", "white_noise"]


def per_star(sm, sel, metric):
    d = sm[sel & sm.wise_usable.astype(bool)].copy()
    d["y"] = np.log(d[metric])
    d = d[np.isfinite(d.y) & np.isfinite(d.disk_excess)]
    d = d[d.groupby("tic").tic.transform("size") >= 2]
    xd = d.disk_excess - d.groupby("tic").disk_excess.transform("mean")
    yd = d.y - d.groupby("tic").y.transform("mean")
    g = pd.DataFrame(dict(tic=d.tic, sxx=xd * xd, sxy=xd * yd,
                          x=d.disk_excess)).groupby("tic").agg(
        sxx=("sxx", "sum"), sxy=("sxy", "sum"), xmin=("x", "min"), xmax=("x", "max"))
    g["slope"] = g.sxy / g.sxx
    g["xr"] = g.xmax - g.xmin
    return g


def boot_beta(gid_sxy, gid_sxx, rng, n=N_BOOT):
    G = len(gid_sxy)
    W = rng.multinomial(G, np.full(G, 1 / G), size=n)
    return (W @ gid_sxy) / (W @ gid_sxx)


def main():
    OUT.mkdir(exist_ok=True)
    sm, sw, disks, dos = T.load()
    samp = pd.read_csv(PILOT / "sample.csv").set_index("tic")
    sid, t, x = sm.tic.values, sm.mid_mjd.values, sm.disk_excess.values
    lab, wu = sm.phase_label.values.astype(object), sm.wise_usable.values.astype(bool)
    vs = T.variants(sm)
    main_m = vs["main (crowd-cut)"]
    spoc_m = vs["SPOC/TESS-SPOC stars only"]
    table = stats4.build_swap_table(t, disks, T.evaluate)
    rows = []

    def run(tag, sel, metric, k=0.0, ids=None, dmap=None, y=None):
        if y is None:
            y = np.where(sel, np.log(sm[metric].values), np.nan) + 0.921 * k * x
        r = stats4.run_test(sid if ids is None else ids, x, y, lab, wu, t, disks,
                            dos if dmap is None else dmap, T.evaluate, n_null=N_NULL,
                            n_boot=N_BOOT, seed=SEED, swap_table=table)
        r = {kk: vv for kk, vv in r.items() if not kk.startswith("_")}
        r.update(check=tag, metric=metric, k=k)
        rows.append(r)
        print(tag, metric, k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in
                               r.items() if kk in ("n_stars", "n_sectors", "beta", "beta_lo",
                                                   "beta_hi", "p_perm", "p_swap")}, flush=True)
        return r

    # ---- M1: Gaia ESP-ELS
    g = pd.read_parquet(CACHE / "gaia" / "gaia_dr3_params.parquet").set_index("source_id")
    lbl = sm.tic.map(samp.gaia_dr3_id).map(g.classlabel_espels).fillna("").values
    be = main_m & (lbl == "beStar")
    nonherb = main_m & ~np.isin(lbl, ["HerbigStar", "wN", "wC"])
    for m in METRICS:
        run("M1 Gaia ESP-ELS beStar", be, m)
    for m in ["A1", "B1", "flicker_rms"]:
        run("M1 main minus ESP-ELS Herbig/WR", nonherb, m)

    # ---- M2: stable f1
    u = sm[main_m & wu & sm.metrics_ok.astype(bool)]
    u = u[u.groupby("tic").tic.transform("size") >= 2]
    rel = (u.f1 / u.groupby("tic").f1.transform("median") - 1).abs()
    stable_tics = set(rel.groupby(u.tic).max().loc[lambda s: s <= 0.10].index)
    stable = main_m & sm.tic.isin(stable_tics).values
    moving = main_m & ~sm.tic.isin(stable_tics).values
    for m in ["A1", "B1", "white_noise"]:
        run("M2 stable-f1 stars", stable, m)
        run("M2 moving-f1 stars", moving, m)
    f1med = u.groupby("tic").f1.median()
    stable_hi = stable & sm.tic.map(f1med).gt(0.8).values
    run("M2 stable-f1 stars, f1 > 0.8 c/d", stable_hi, "A1")

    # ---- M3: SPOC with measured k
    for m in ["A1", "B1", "red_noise", "flicker_rms", "white_noise"]:
        run("M3 SPOC/TESS-SPOC, k measured", spoc_m, m, k=K_MEAS)

    # ---- M4: white-noise terciles (star median over main usable sectors)
    wmed = u.groupby("tic").white_noise.median()
    q = wmed.quantile([1 / 3, 2 / 3]).values
    terc = pd.cut(wmed, [-np.inf, q[0], q[1], np.inf], labels=["quietest", "middle", "noisiest"])
    tstar = sm.tic.map(terc).astype(object).values
    for name in ["quietest", "middle", "noisiest"]:
        sel = main_m & (tstar == name)
        for m in ["white_noise", "A1", "B1", "flicker_rms"]:
            r = run(f"M4 white-noise tercile: {name}", sel, m)
            r["tercile_median_white_ppt"] = float(wmed[terc == name].median())

    # ---- M5: drop top-N leverage stars (within-star sum of squared excess)
    lev = per_star(sm, main_m, "A1").sxx.sort_values(ascending=False)
    for n in (10, 50, 100):
        sel = main_m & ~sm.tic.isin(lev.index[:n]).values
        for m in ["A1", "B1", "flicker_rms"]:
            run(f"M5 drop top {n} leverage stars", sel, m)

    # ---- S3: star x cadence groups
    cad = sm.cadence_s.round().fillna(-1).astype(int).values
    comp = sid.astype(np.int64) * 10000 + cad
    dmap = {c: dos[s] for c, s in zip(comp, sid)}
    for m in ["A1", "B1", "flicker_rms", "white_noise"]:
        run("S3 within star x cadence", main_m, m, ids=comp, dmap=dmap)

    # ---- S1: recentred relative bins
    S = F.load_spectra(sm)
    lg = np.log10(LOG_GRID)
    f1s = sm[sm.usable.astype(bool)].groupby("tic").f1.median()
    lrel = lg[None, :] - np.log10(sm.tic.map(f1s).values)[:, None]
    e_rel = np.round(np.arange(-1.25, 1.2501, 0.1), 4)
    Y = F.rebin(S, lrel, e_rel)
    for vname, sel in (("main", main_m), ("SPOC", spoc_m), ("stable-f1", stable)):
        for b in range(Y.shape[1]):
            y = np.where(sel, np.log(Y[:, b]), np.nan)
            r = run(f"S1 relative bin {vname}", sel, "rel_amp", y=y)
            r.update(lo=10 ** e_rel[b], hi=10 ** e_rel[b + 1],
                     centre=10 ** (0.5 * (e_rel[b] + e_rel[b + 1])))
    out = pd.DataFrame(rows)
    front = ["check", "metric", "k", "n_stars", "n_sectors", "beta", "beta_lo", "beta_hi",
             "p_perm", "p_swap", "swap_null_mean"]
    out = out[front + [c for c in out.columns if c not in front]]
    out.to_csv(OUT / "review_checks.csv", index=False)

    # ---- bootstrap-only extras
    rng = np.random.default_rng(SEED)
    ext = []
    d = sm[main_m & wu].copy()
    d = d[d.groupby("tic").tic.transform("size") >= 2]
    for m in ["A1", "B1", "B2", "B3", "red_noise", "flicker_rms", "white_noise"]:
        dd = d[np.isfinite(np.log(d[m]))]
        xx, yy = dd.disk_excess.values.copy(), np.log(dd[m].values)
        for _ in range(60):   # two-way within transformation, alternating projections
            for gcol in ("tic", "sector"):
                xx = xx - pd.Series(xx).groupby(dd[gcol].values).transform("mean").values
                yy = yy - pd.Series(yy).groupby(dd[gcol].values).transform("mean").values
        gs = pd.DataFrame(dict(tic=dd.tic.values, sxy=xx * yy, sxx=xx * xx)).groupby("tic").sum()
        bs = boot_beta(gs.sxy.values, gs.sxx.values, rng)
        ext.append(dict(extra="two-way star+sector FE", metric=m,
                        beta=gs.sxy.sum() / gs.sxx.sum(), beta_lo=np.percentile(bs, 16),
                        beta_hi=np.percentile(bs, 84), n_stars=len(gs)))
    # ln white on excess + ln B1 (within star)
    dd = d.copy()
    for c, src in (("yw", "white_noise"), ("lb", "B1")):
        dd[c] = np.log(dd[src])
    for c in ("yw", "lb", "disk_excess"):
        dd[c + "_d"] = dd[c] - dd.groupby("tic")[c].transform("mean")
    Xm = np.c_[dd.disk_excess_d, dd.lb_d]
    coef = np.linalg.lstsq(Xm, dd.yw_d, rcond=None)[0]
    tic_codes, gi = np.unique(dd.tic.values, return_inverse=True)
    bs = []
    for _ in range(N_BOOT):
        w = np.bincount(rng.integers(0, len(tic_codes), len(tic_codes)),
                        minlength=len(tic_codes))[gi]
        sw_ = np.sqrt(w)
        bs.append(np.linalg.lstsq(Xm * sw_[:, None], dd.yw_d * sw_, rcond=None)[0])
    bs = np.array(bs)
    ext.append(dict(extra="ln white ~ excess + ln B1: excess coef", metric="white_noise",
                    beta=coef[0], beta_lo=np.percentile(bs[:, 0], 16),
                    beta_hi=np.percentile(bs[:, 0], 84), n_stars=len(tic_codes)))
    ext.append(dict(extra="ln white ~ excess + ln B1: ln B1 coef", metric="white_noise",
                    beta=coef[1], beta_lo=np.percentile(bs[:, 1], 16),
                    beta_hi=np.percentile(bs[:, 1], 84), n_stars=len(tic_codes)))
    ext.append(dict(extra="within-star corr(ln white, ln B1)", metric="white_noise",
                    beta=float(np.corrcoef(dd.yw_d, dd.lb_d)[0, 1])))
    # QLP - SPOC gaps (bootstrap over stars, independent samples)
    qlp_m = vs["QLP stars only"]
    for m in ["A1", "B1", "flicker_rms"]:
        a, b = per_star(sm, qlp_m, m), per_star(sm, spoc_m, m)
        ba, bb = boot_beta(a.sxy.values, a.sxx.values, rng), boot_beta(b.sxy.values, b.sxx.values, rng)
        gap = ba - bb
        ext.append(dict(extra="QLP - SPOC beta (k=0)", metric=m,
                        beta=a.sxy.sum() / a.sxx.sum() - b.sxy.sum() / b.sxx.sum(),
                        beta_lo=np.percentile(gap, 16), beta_hi=np.percentile(gap, 84),
                        sd=gap.std(), predicted_by_k=0.921 * K_MEAS))
    pd.DataFrame(ext).to_csv(OUT / "review_extras.csv", index=False)
    print(pd.DataFrame(ext).round(3).to_string(index=False))

    # ---- per-star slope signs by excess range
    srows = []
    bins = [0, 0.05, 0.1, 0.2, 0.4, np.inf]
    labels = ["<0.05", "0.05-0.1", "0.1-0.2", "0.2-0.4", ">0.4"]
    gA = per_star(sm, main_m, "A1")
    gb = pd.cut(gA.xr, bins, labels=labels, right=False)
    tot_sxx = gA.sxx.sum()
    for m in ["A1", "B1", "flicker_rms"]:
        gm = per_star(sm, main_m, m)
        for lb_ in labels:
            s = gm.loc[gb.index[gb == lb_].intersection(gm.index)]
            srows.append(dict(metric=m, x_range=lb_, n_stars=len(s),
                              frac_positive=float(np.mean(s.slope > 0)),
                              share_of_sxx=float(gA.loc[s.index].sxx.sum() / tot_sxx)))
    lev_c = gA.sxx.sort_values(ascending=False).cumsum() / tot_sxx
    signs = pd.DataFrame(srows)
    signs.to_csv(OUT / "review_signs.csv", index=False)
    print(signs.round(3).to_string(index=False))
    print("stars with x range < 0.1 mag:", int((gA.xr < 0.1).sum()), "of", len(gA),
          "; median x range", round(float(gA.xr.median()), 3),
          "; stars carrying 50% / 90% of sum x^2:", int((lev_c < 0.5).sum()) + 1,
          int((lev_c < 0.9).sum()) + 1)
    print("Gaia ESP-ELS labels among main stars:",
          pd.Series(lbl[main_m]).groupby(sm.tic.values[main_m]).first().value_counts().to_dict())
    print("stable-f1 stars:", len(stable_tics))


if __name__ == "__main__":
    main()
