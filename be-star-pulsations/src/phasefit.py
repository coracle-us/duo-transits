"""Within-star OLS of ln(metric) on [disk_excess, building, dissipating] (checks/c1, c2).

All columns are demeaned within each star (stars need >= 2 kept sectors). Nulls:
  - within-star permutation: shuffle (disk_excess, label) jointly among a star's sectors;
  - cross-star swap: give each star another star's WISE curve, evaluated at its own sector
    times (the stats4 swap table), keeping only sectors usable under the donor curve.
Several metrics (columns of Y) are fitted at once; every metric must be finite on the
same sectors (true for sector_metrics.csv: metrics_ok rows are all finite).
"""
import numpy as np
import stats4

TERMS = ["disk_excess", "building", "dissipating"]
B, D = stats4.LABELS.index("building"), stats4.LABELS.index("dissipating")


def design(x, lc):
    return np.column_stack([x, lc == B, lc == D]).astype(float)


def fit(gid, X, Y, w=None, G=None, per_star=False):
    """Coefficients (3 x n_metric) of the within-star demeaned OLS with 0/1 weights w."""
    G = G if G is not None else gid.max() + 1
    w = np.ones(len(gid)) if w is None else w.astype(float)
    n = np.bincount(gid, weights=w, minlength=G)
    w = w * (n[gid] >= 2)
    n = np.maximum(np.bincount(gid, weights=w, minlength=G), 1)

    def dm(M):
        mean = np.stack([np.bincount(gid, weights=w * M[:, i], minlength=G)
                         for i in range(M.shape[1])], 1) / n[:, None]
        return (M - mean[gid]) * w[:, None]
    Xd, Yd = dm(X), dm(Y)
    xtx, xty = Xd.T @ Xd, Xd.T @ Yd
    try:
        coef = np.linalg.solve(xtx, xty)
    except np.linalg.LinAlgError:
        coef = np.full((X.shape[1], Y.shape[1]), np.nan)
    if not per_star:
        return coef
    # per-star cross-product sums for the bootstrap over stars
    P = Xd[:, :, None] * Xd[:, None, :]
    Q = Xd[:, :, None] * Yd[:, None, :]
    sxx = np.stack([np.bincount(gid, weights=P[:, a, b], minlength=G)
                    for a in range(3) for b in range(3)], 1).reshape(G, 3, 3)
    sxy = np.stack([np.bincount(gid, weights=Q[:, a, m], minlength=G)
                    for a in range(3) for m in range(Y.shape[1])], 1).reshape(G, 3, Y.shape[1])
    return coef, sxx, sxy


def bootstrap(sxx, sxy, n, rng):
    G = len(sxx)
    W = rng.multinomial(G, np.full(G, 1.0 / G), size=n).astype(float)
    A = np.einsum("kg,gab->kab", W, sxx)
    C = np.einsum("kg,gam->kam", W, sxy)
    out = np.full((n,) + sxy.shape[1:], np.nan)
    ok = np.abs(np.linalg.det(A)) > 0
    out[ok] = np.linalg.solve(A[ok], C[ok])
    return np.nanpercentile(out, [16, 84], axis=0)


def run(star_ids, x, lab, usable, Y, table, disk_of_star, n_null=10000, n_boot=2000,
        seed=20261104):
    """Y: (n_sector, n_metric) with NaN rows outside the variant. Returns dict of arrays."""
    rng = np.random.default_rng(seed)
    lc_all = stats4.code(lab)
    fin = np.all(np.isfinite(Y), axis=1)
    ids, gid = stats4.compact(star_ids, usable & fin & np.isfinite(x))
    G = int(gid.max() + 1)
    X = design(x[ids], lc_all[ids])
    Ys = Y[ids]
    coef, sxx, sxy = fit(gid, X, Ys, G=G, per_star=True)
    lo, hi = bootstrap(sxx, sxy, n_boot, rng)
    xs, ls = x[ids], lc_all[ids]
    perm = np.empty((n_null,) + coef.shape)
    for k in range(n_null):
        o = np.lexsort((rng.random(len(gid)), gid))
        perm[k] = fit(gid, design(xs[o], ls[o]), Ys, G=G)
    # swap: all finite sectors of the test stars (as stats4.run_test), donor-usable only
    XT, LT, UT = table
    Dn = XT.shape[1]
    stars = np.unique(star_ids[ids])
    sel = np.where(np.isin(star_ids, stars) & fin)[0]
    _, sidx = np.unique(star_ids[sel], return_inverse=True)
    o = np.argsort(sidx, kind="stable")
    sel, sidx = sel[o], sidx[o]
    own = np.array([disk_of_star[s] for s in stars])
    Ysel = Y[sel]
    swap = np.empty_like(perm)
    for k in range(n_null):
        donors = rng.integers(0, Dn - 1, size=G)
        donors += donors >= own
        dj = donors[sidx]
        Xk = design(XT[sel, dj].astype(float), LT[sel, dj])
        swap[k] = fit(sidx, Xk, Ysel, w=UT[sel, dj], G=G)
    p_perm = np.full(coef.shape, np.nan)
    p_swap = p_perm.copy()
    for a in range(coef.shape[0]):
        for m in range(coef.shape[1]):
            p_perm[a, m] = stats4.pval(coef[a, m], perm[:, a, m])
            p_swap[a, m] = stats4.pval(coef[a, m], swap[:, a, m])
    return dict(n_stars=G, n_sectors=len(ids), n_build=int(np.sum(ls == B)),
                n_diss=int(np.sum(ls == D)), coef=coef, lo=lo, hi=hi, p_perm=p_perm,
                p_swap=p_swap, swap_mean=np.nanmean(swap, axis=0))
