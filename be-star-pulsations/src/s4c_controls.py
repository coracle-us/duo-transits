"""Stage 4c (post hoc, full run): is the white_noise 'control' really instrumental?

In the full run the white_noise metric (median amplitude at 15-24 c/d of the 30-min-binned
spectrum) also rises with disk excess, so it is not a clean data-quality control: Be-star
stochastic variability reaches that band. This stage adds:

1. err_ppt: the pipeline's formal relative flux error per sector (median over quality-0
   points of PDCSAP_FLUX_ERR / KSPSAP_FLUX_ERR / FLUX_ERR, over the median flux), a purely
   instrumental control (photon, read and background noise; tracks brightness).
2. Each metric over white_noise, i.e. whether the metric rises faster than the
   high-frequency floor (beta of the ratio = beta_metric - beta_white; the nulls are new).

Same within-star and cross-star nulls as stage 4, k = 0 only.
"""
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd
from astropy.io import fits
from common import PILOT, TESS_DIR, SEED
import s4_test as T, stats4

ERR_COL = {"SPOC120": ("PDCSAP_FLUX", "PDCSAP_FLUX_ERR"),
           "TESS-SPOC": ("PDCSAP_FLUX", "PDCSAP_FLUX_ERR"),
           "QLP": ("KSPSAP_FLUX", "KSPSAP_FLUX_ERR"),
           "GSFC-ELEANOR-LITE": ("CORR_FLUX", "FLUX_ERR")}


def err_ppt(args):
    path, product = args
    fcol, ecol = ERR_COL[product]
    try:
        with fits.open(TESS_DIR / path, memmap=False) as h:
            d = h[1].data
            f, e = np.asarray(d[fcol], float), np.asarray(d[ecol], float)
            q = np.asarray(d["QUALITY"]) if "QUALITY" in d.columns.names else np.zeros(len(f))
        ok = np.isfinite(f) & np.isfinite(e) & (q == 0) & (e > 0)
        return 1e3 * np.median(e[ok]) / np.median(f[ok]) if ok.sum() > 50 else np.nan
    except Exception:
        return np.nan


def main():
    sm, sw, disks, dos = T.load()
    with ProcessPoolExecutor(6) as ex:
        sm["err_ppt"] = list(ex.map(err_ppt, zip(sm.file, sm["product"]), chunksize=50))
    sm[["tic", "sector", "product", "err_ppt"]].to_csv(PILOT / "sector_controls.csv", index=False)
    print("err_ppt finite:", int(np.isfinite(sm.err_ppt).sum()), "of", len(sm))
    print(sm.groupby("product").err_ppt.median().round(3).to_dict())
    sid, t, x = sm.tic.values, sm.mid_mjd.values, sm.disk_excess.values
    lab, wu = sm.phase_label.values.astype(object), sm.wise_usable.values.astype(bool)
    table = stats4.build_swap_table(t, disks, T.evaluate)
    ys = {"err_ppt (instrumental)": np.log(sm.err_ppt.values)}
    for m in ["A1", "B1", "B2", "B3", "red_noise", "flicker_rms"]:
        ys[f"{m}/white_noise"] = np.log(sm[m].values / sm.white_noise.values)
    vs = T.variants(sm)
    rows = []
    for vname in ["main (crowd-cut)", "no crowding cut", "SPOC/TESS-SPOC stars only",
                  "QLP stars only"]:
        sel = vs[vname]
        for name, y0 in ys.items():
            y = np.where(sel, y0, np.nan)
            r = stats4.run_test(sid, x, y, lab, wu, t, disks, dos, T.evaluate, n_null=10000,
                                n_boot=2000, seed=SEED, swap_table=table)
            r.update(variant=vname, metric=name)
            rows.append(r)
            print(vname, name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()
                                if k in ("n_stars", "beta", "beta_lo", "beta_hi", "p_perm", "p_swap")},
                  flush=True)
    cols = ["variant", "metric", "n_stars", "n_sectors", "beta", "beta_lo", "beta_hi", "rho",
            "contrast", "p_perm", "p_swap", "p_perm_contrast", "p_swap_contrast", "swap_null_mean"]
    pd.DataFrame(rows)[cols].to_csv(PILOT / "control_tests.csv", index=False)


if __name__ == "__main__":
    main()
