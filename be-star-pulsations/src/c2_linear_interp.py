"""Check 2: disk excess without the spline.

w1_fit at each sector mid-time = straight line between the two bracketing WISE visits
(wise.LinearDiskCurve); a sector is dropped outside the WISE span or if either bracketing
visit is > 200 d away; w1_quiet = faintest visit. Visit values are those of the stage 2
visit tables (CACHE/wise/visits/<tic>.csv, w1mpro_ep: the China-VO per-epoch value, and for
the IRSA extension the per-visit median of single exposures). Reruns the stage 4 main test
(k = 0, both nulls) and the Check 1 phase-term fit on this disk_excess, main sample.
Output: results/checks/linear_interp.csv. Run with BE_RUN=full.
"""
import time
import numpy as np, pandas as pd
from common import RESULTS, PILOT, WISE_DIR, SEED
from wise import LinearDiskCurve
import s4_test as T, stats4, c1_phase_terms as C1

OUT = RESULTS / "checks"


def main():
    sub, vsub, _, dos = C1.main_subset()
    tics = sorted(dos, key=dos.get)
    lin = []
    for tic in tics:
        v = pd.read_csv(WISE_DIR / "visits" / f"{tic}.csv")
        lin.append(LinearDiskCurve(v.mjd.values, v.w1mpro_ep.values, v.w1_w2.values))
    x = np.full(len(sub), np.nan)
    lab = np.empty(len(sub), dtype=object)
    use = np.zeros(len(sub), bool)
    for tic, g in sub.groupby("tic").groups.items():
        o = lin[dos[tic]].at(sub.mid_mjd.values[g])
        x[g], lab[g], use[g] = o["disk_excess"], o["phase_label"], o["wise_usable"]
    print(f"linear: usable {use.sum()} of {len(sub)} main rows "
          f"(spline usable {int(sub.wise_usable.sum())}); corr(x_lin, x_spline) on both-usable "
          f"{np.corrcoef(x[use & sub.wise_usable.values], sub.disk_excess.values[use & sub.wise_usable.values])[0, 1]:.3f}",
          flush=True)
    table = stats4.build_swap_table(sub.mid_mjd.values, lin, T.evaluate)
    spl = pd.read_csv(PILOT / "test_results.csv")
    spl = spl[(spl.variant == "main (crowd-cut)") & (spl.k == 0)].set_index("metric")
    main = vsub["main (crowd-cut)"]
    rows = []
    for metric in T.METRICS:
        y = np.where(main, np.log(sub[metric].values), np.nan)
        r = stats4.run_test(sub.tic.values, x, y, lab, use, sub.mid_mjd.values, lin, dos,
                            T.evaluate, n_null=T.N_NULL, n_boot=T.N_BOOT, seed=SEED,
                            swap_table=table)
        s = spl.loc[metric]
        rows.append(dict(test="stage4 main", metric=metric, term="disk_excess",
                         n_stars=r["n_stars"], n_sectors=r["n_sectors"], coef=r["beta"],
                         coef_lo=r["beta_lo"], coef_hi=r["beta_hi"], p_perm=r["p_perm"],
                         p_swap=r["p_swap"], spline_coef=s.beta, spline_lo=s.beta_lo,
                         spline_hi=s.beta_hi, spline_p_perm=s.p_perm, spline_p_swap=s.p_swap,
                         ratio_lin_spline=r["beta"] / s.beta))
        print(f"{metric}: beta linear {r['beta']:.3f} [{r['beta_lo']:.3f}, {r['beta_hi']:.3f}] "
              f"p {r['p_perm']:.4f} {r['p_swap']:.4f} | spline {s.beta:.3f} "
              f"({r['n_stars']} stars, {r['n_sectors']} sectors)", flush=True)
    ph = pd.DataFrame(C1.phase_rows(sub, {"main (crowd-cut)": main}, lin, dos, x, lab, use,
                                    table, "linear"))
    sp = pd.read_csv(OUT / "phase_terms.csv")
    sp = sp[(sp.disk_model == "spline") & (sp.variant == "main (crowd-cut)")]
    ph = ph.merge(sp[["metric", "term", "coef", "coef_lo", "coef_hi", "p_perm", "p_swap"]]
                  .rename(columns=lambda c: c if c in ("metric", "term") else "spline_" + c),
                  on=["metric", "term"])
    ph["test"] = "check1 phase-term fit"
    ph["ratio_lin_spline"] = ph.coef / ph.spline_coef
    out = pd.concat([pd.DataFrame(rows), ph.drop(columns=["disk_model", "variant", "n_build",
                                                          "n_diss", "swap_null_mean"])],
                    ignore_index=True)
    out.to_csv(OUT / "linear_interp.csv", index=False)
    print(out[["test", "metric", "term", "coef", "spline_coef", "p_perm", "p_swap"]]
          .round(4).to_string(index=False))


if __name__ == "__main__":
    main()
