"""Within-star test of a metric against disk excess (stage 4 helpers), vectorised.

Per-sector arrays: star id, x = disk_excess, y = ln(metric), lab = phase label, t = mid MJD.
The cross-star null re-evaluates a donor star's DiskCurve at this star's own sector times;
`evaluate(disk, t)` must return (x, lab, usable) arrays. All (sector, donor) evaluations are
computed once (build_swap_table) and reused across metrics and variants.
"""
import numpy as np
from scipy.stats import spearmanr

LABELS = ["building", "dissipating", "plateau", "quiescent"]
A_CODE, B_CODE = 0, 1   # contrast = building - dissipating


def code(lab):
    lab = np.asarray(lab, dtype=object)
    out = np.full(len(lab), -1, np.int8)
    for i, l in enumerate(LABELS):
        out[lab == l] = i
    return out


def compact(star_ids, keep):
    """Indices of kept sectors of stars with >= 2 kept sectors, sorted by star; group ids."""
    ids = np.where(keep)[0]
    sid = star_ids[ids]
    _, inv, cnt = np.unique(sid, return_inverse=True, return_counts=True)
    ok = cnt[inv] >= 2
    ids, sid = ids[ok], sid[ok]
    o = np.argsort(sid, kind="stable")
    ids, sid = ids[o], sid[o]
    _, gid = np.unique(sid, return_inverse=True)
    return ids, gid


def _stats(gid, x, y, lc, w=None, G=None):
    """beta and contrast with optional 0/1 sector mask w (stars need >= 2 kept sectors)."""
    G = G if G is not None else gid.max() + 1
    if w is None:
        w = np.ones(len(x))
    n = np.bincount(gid, weights=w, minlength=G)
    w = w * (n[gid] >= 2)
    n = np.bincount(gid, weights=w, minlength=G)
    nz = np.maximum(n, 1)
    xd = (x - (np.bincount(gid, weights=w * x, minlength=G) / nz)[gid]) * w
    yd = (y - (np.bincount(gid, weights=w * y, minlength=G) / nz)[gid]) * w
    sxx = np.sum(xd * xd)
    beta = np.sum(xd * yd) / sxx if sxx > 0 else np.nan
    a, b = (lc == A_CODE) & (w > 0), (lc == B_CODE) & (w > 0)
    contrast = yd[a].mean() - yd[b].mean() if a.any() and b.any() else np.nan
    return beta, contrast, xd, yd


def pval(obs, null):
    """Two-sided, about the null's own centre (so a shifted null is handled)."""
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or len(null) == 0:
        return np.nan
    c = np.mean(null)
    return (1 + np.sum(np.abs(null - c) >= np.abs(obs - c))) / (1 + len(null))


def perm_null(gid, x, y, lc, n, rng):
    """Shuffle x and labels among each star's own sectors (gid must be sorted)."""
    bs, cs = np.empty(n), np.empty(n)
    for k in range(n):
        o = np.lexsort((rng.random(len(gid)), gid))
        bs[k], cs[k], _, _ = _stats(gid, x[o], y, lc[o])
    return bs, cs


def bootstrap(gid, xd, yd, n, rng):
    """Resample stars; demeaning is within star, so duplicates act as integer weights."""
    G = gid.max() + 1
    sxy = np.bincount(gid, weights=xd * yd, minlength=G)
    sxx = np.bincount(gid, weights=xd * xd, minlength=G)
    W = rng.multinomial(G, np.full(G, 1.0 / G), size=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        bs = (W @ sxy) / (W @ sxx)
    return np.nanpercentile(bs, [16, 84])


def build_swap_table(t, disks, evaluate):
    """X[s, d], L[s, d] (label code), U[s, d] for every sector s under every disk d."""
    S, D = len(t), len(disks)
    X = np.empty((S, D), np.float32)
    L = np.empty((S, D), np.int8)
    U = np.empty((S, D), bool)
    for j, d in enumerate(disks):
        xx, ll, uu = evaluate(d, t)
        X[:, j], L[:, j], U[:, j] = xx, code(ll), uu
    return X, L, U


def swap_null(sel, star_ids, y, table, disk_of_star, n, rng):
    """Give each star the WISE disk curve of another random star, at its own sector times.
    sel: sector indices (all finite-y sectors of the stars in the test)."""
    X, L, U = table
    D = X.shape[1]
    sid = star_ids[sel]
    stars, sidx = np.unique(sid, return_inverse=True)
    own = np.array([disk_of_star[s] for s in stars])
    o = np.argsort(sidx, kind="stable")
    sel, sidx = sel[o], sidx[o]
    ys = y[sel]
    G = len(stars)
    bs, cs = np.empty(n), np.empty(n)
    for k in range(n):
        donors = rng.integers(0, D - 1, size=G)
        donors += donors >= own
        dj = donors[sidx]
        xs = X[sel, dj].astype(float)
        w = U[sel, dj].astype(float)
        bs[k], cs[k], _, _ = _stats(sidx, xs, ys, L[sel, dj], w=w, G=G)
    return bs, cs


def run_test(star_ids, x, y, lab, usable, t, disks, disk_of_star, evaluate,
             n_null=10000, n_boot=2000, seed=20261104, do_swap=True, keep_nulls=False,
             swap_table=None):
    rng = np.random.default_rng(seed)
    lc_all = code(lab)
    ids, gid = compact(star_ids, usable & np.isfinite(y) & np.isfinite(x))
    out = dict(n_stars=int(gid.max() + 1) if len(gid) else 0, n_sectors=len(ids))
    if out["n_stars"] < 3:
        return out
    xs, ys, ls = x[ids], y[ids], lc_all[ids]
    beta, contrast, xd, yd = _stats(gid, xs, ys, ls)
    pb, pc = perm_null(gid, xs, ys, ls, n_null, rng)
    lo, hi = bootstrap(gid, xd, yd, n_boot, rng)
    out.update(beta=beta, beta_lo=lo, beta_hi=hi, rho=spearmanr(xd, yd).correlation,
               contrast=contrast, n_build=int(np.sum(ls == A_CODE)),
               n_diss=int(np.sum(ls == B_CODE)),
               p_perm=pval(beta, pb), p_perm_contrast=pval(contrast, pc))
    if keep_nulls:
        out["_perm_null"], out["_xd"], out["_yd"] = pb, xd, yd
    if do_swap:
        if swap_table is None:
            swap_table = build_swap_table(t, disks, evaluate)
        stars = np.unique(star_ids[ids])
        sel = np.where(np.isin(star_ids, stars) & np.isfinite(y))[0]
        sb, sc = swap_null(sel, star_ids, y, swap_table, disk_of_star, n_null, rng)
        out.update(p_swap=pval(beta, sb), p_swap_contrast=pval(contrast, sc),
                   swap_null_mean=float(np.nanmean(sb)))
        if keep_nulls:
            out["_swap_null"] = sb
    return out
