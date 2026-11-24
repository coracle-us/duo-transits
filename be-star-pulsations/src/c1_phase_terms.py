"""Check 1: cross-star swap null for the disk-phase terms.

For each metric: y = ln(metric) on [disk_excess, building, dissipating], all demeaned within
star, main crowd-cut usable sample (and SPOC/TESS-SPOC-only, QLP-only). Within-star
permutation and cross-star swap nulls, 10,000 draws each; the swap re-evaluates the donor's
WISE spline at this star's sector times with the stage 2 code (wise.DiskCurve.at, through
s4_test.evaluate). Output: results/checks/phase_terms.csv. Run with BE_RUN=full.
"""
import time
import numpy as np, pandas as pd
from common import RESULTS, SEED
import s4_test as T, stats4, phasefit

OUT = RESULTS / "checks"
OUT.mkdir(exist_ok=True)
METRICS = T.METRICS
N_NULL, N_BOOT = 10000, 2000


def main_subset():
    """Stage-4 inputs restricted to the main crowd-cut rows (all variants here are subsets)."""
    sm, sw, disks, dos = T.load()
    vs = T.variants(sm)
    keep = vs["main (crowd-cut)"]
    sub = sm[keep].reset_index(drop=True)
    vsub = {k: v[keep] for k, v in vs.items()
            if k in ("main (crowd-cut)", "SPOC/TESS-SPOC stars only", "QLP stars only")}
    return sub, vsub, disks, dos


def phase_rows(sub, vsub, disks, dos, x, lab, usable, table, tag):
    Yall = np.log(sub[METRICS].values.astype(float))
    rows = []
    for vname, sel in vsub.items():
        t0 = time.time()
        Y = np.where(sel[:, None], Yall, np.nan)
        r = phasefit.run(sub.tic.values, x, lab, usable, Y, table, dos, n_null=N_NULL,
                         n_boot=N_BOOT, seed=SEED)
        for m, metric in enumerate(METRICS):
            for a, term in enumerate(phasefit.TERMS):
                rows.append(dict(disk_model=tag, variant=vname, metric=metric, term=term,
                                 n_stars=r["n_stars"], n_sectors=r["n_sectors"],
                                 n_build=r["n_build"], n_diss=r["n_diss"],
                                 coef=r["coef"][a, m], coef_lo=r["lo"][a, m],
                                 coef_hi=r["hi"][a, m],
                                 p_perm=r["p_perm"][a, m], p_swap=r["p_swap"][a, m],
                                 swap_null_mean=r["swap_mean"][a, m]))
        print(tag, vname, f"{r['n_stars']} stars, {time.time()-t0:.0f} s", flush=True)
    return rows


def main():
    sub, vsub, disks, dos = main_subset()
    t0 = time.time()
    table = stats4.build_swap_table(sub.mid_mjd.values, disks, T.evaluate)
    print(f"swap table {table[0].shape}, {time.time()-t0:.0f} s", flush=True)
    rows = phase_rows(sub, vsub, disks, dos, sub.disk_excess.values,
                      sub.phase_label.values.astype(object),
                      sub.wise_usable.values.astype(bool), table, "spline")
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "phase_terms.csv", index=False)
    show = out[out.term != "disk_excess"]
    print(show[["variant", "metric", "term", "coef", "p_perm", "p_swap"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
