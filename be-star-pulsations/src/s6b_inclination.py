"""Stage 6b (follow-up, post hoc): does beta depend on viewing angle?

Proxy: Gaia DR3 ESP-HS vsini (gaia_dr3_astrophys_param at IRSA, fetched by source_id;
cache ~/astro-cache/be-star-pulsations/gaia/gaia_dr3_params.parquet). Be stars rotate near
critical, so v sin i tracks sin i. Gaia vbroad is not used: 85 of 2,541 stars, fitted with
cool templates.

Quality: vsini_esphs depends on brightness in this sample (median 244 km/s at G < 9, about
65-70 km/s at G > 12, 30% below 10 km/s there), a signal-to-noise failure. Primary split:
stars with G <= 11, at the median vsini of the stars in the test. Sensitivity: all stars
with vsini, split at the median within each 1-mag G bin.

Prediction, written before running: if the disk acts on the photometry mainly through its
own optical light (diluting the pulsation, low i) or by absorbing starlight (shell stars,
high i), the observed beta should be larger for high v sin i. If the coupling is intrinsic
to the star, no difference.

Test: beta_high - beta_low. Null: v sin i group labels shuffled among the stars in the test
(10,000 draws); 68% interval: bootstrap over stars within each group. Each group's own beta
also gets the stage-4 nulls. Output: results/<run>/inclination_split.csv.
"""
import os
import numpy as np, pandas as pd
from common import PILOT, SEED, CACHE
import s4_test as T, stats4

METRICS = ["A1", "B1", "B2", "B3", "red_noise", "flicker_rms", "white_noise"]
N_NULL = 10000


def star_sums(sid, x, y, keep):
    """Per-star sums of demeaned x*y and x*x over kept sectors (stars with >= 2)."""
    ids, gid = stats4.compact(sid, keep & np.isfinite(y) & np.isfinite(x))
    if len(ids) == 0:
        return np.array([]), np.array([]), np.array([])
    _, _, xd, yd = stats4._stats(gid, x[ids], y[ids], np.zeros(len(ids), np.int8))
    G = gid.max() + 1
    stars = np.unique(sid[ids])          # compact sorts by star, so gid order = sorted ids
    return stars, np.bincount(gid, xd * yd, G), np.bincount(gid, xd * xd, G)


def diff_test(hi, sxy, sxx, rng):
    def b(m):
        return sxy[m].sum() / sxx[m].sum()
    d = b(hi) - b(~hi)
    null = np.empty(N_NULL)
    for k in range(N_NULL):
        p = rng.permutation(hi)
        null[k] = b(p) - b(~p)
    boot = np.empty(2000)
    ih, il = np.where(hi)[0], np.where(~hi)[0]
    for k in range(2000):
        bh, bl = rng.choice(ih, len(ih)), rng.choice(il, len(il))
        boot[k] = sxy[bh].sum() / sxx[bh].sum() - sxy[bl].sum() / sxx[bl].sum()
    p = (1 + np.sum(np.abs(null - null.mean()) >= abs(d - null.mean()))) / (1 + N_NULL)
    return d, np.percentile(boot, [16, 84]), p


def main():
    sm, sw, disks, dos = T.load()
    samp = pd.read_csv(PILOT / "sample.csv").set_index("tic")
    g = pd.read_parquet(CACHE / "gaia" / "gaia_dr3_params.parquet").set_index("source_id")
    gid_of = samp.gaia_dr3_id
    vs_star = gid_of.map(g.vsini_esphs)
    gmag = gid_of.map(g.phot_g_mean_mag)
    sid, t, x = sm.tic.values, sm.mid_mjd.values, sm.disk_excess.values
    lab, wu = sm.phase_label.values.astype(object), sm.wise_usable.values.astype(bool)
    table = stats4.build_swap_table(t, disks, T.evaluate)
    vs = T.variants(sm)
    rng = np.random.default_rng(SEED)
    rows = []
    for split in ("G<=11, median", "all G, median within 1-mag G bins"):
        for vname in ("main (crowd-cut)", "no crowding cut"):
            for m in METRICS:
                y0 = np.log(sm[m].values)
                keep = vs[vname] & wu
                stars, sxy, sxx = star_sums(sid, x, y0, keep)
                v = vs_star.reindex(stars).values
                gm = gmag.reindex(stars).values
                ok = np.isfinite(v) & (sxx > 0)
                if split.startswith("G<=11"):
                    ok &= gm <= 11
                    hi = v >= np.median(v[ok])
                else:
                    gb = np.floor(np.clip(gm, 8, 14))
                    med = pd.Series(v[ok]).groupby(gb[ok]).median()
                    hi = v >= pd.Series(gb).map(med).values
                stars, sxy, sxx, hi, v = stars[ok], sxy[ok], sxx[ok], hi[ok], v[ok]
                d, (dlo, dhi), p = diff_test(hi, sxy, sxx, rng)
                row = dict(split=split, variant=vname, metric=m, n_hi=int(hi.sum()),
                           n_lo=int((~hi).sum()), vsini_med_hi=float(np.median(v[hi])),
                           vsini_med_lo=float(np.median(v[~hi])), diff=d, diff_lo=dlo,
                           diff_hi=dhi, p_diff=p)
                for name, grp in (("highv", stars[hi]), ("lowv", stars[~hi])):
                    y = np.where(vs[vname] & np.isin(sid, grp), y0, np.nan)
                    r = stats4.run_test(sid, x, y, lab, wu, t, disks, dos, T.evaluate,
                                        n_null=N_NULL, n_boot=2000, seed=SEED,
                                        swap_table=table)
                    for k in ("beta", "beta_lo", "beta_hi", "p_perm", "p_swap"):
                        row[f"{k}_{name}"] = r.get(k)
                rows.append(row)
                print(split, vname, m, {k: (round(vv, 3) if isinstance(vv, float) else vv)
                                        for k, vv in row.items() if k in
                                        ("n_hi", "n_lo", "beta_highv", "beta_lowv", "diff", "p_diff")},
                      flush=True)
    pd.DataFrame(rows).to_csv(PILOT / "inclination_split.csv", index=False)


if __name__ == "__main__":
    main()
