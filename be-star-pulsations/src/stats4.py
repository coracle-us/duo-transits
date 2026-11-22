"""Within-star test of a metric against disk excess (stage 4 helpers).

Per-sector arrays: star id, x = disk_excess, y = ln(metric), lab = phase label, t = mid MJD.
The cross-star null re-evaluates a donor star's DiskCurve at this star's own sector times;
`evaluate(disk, t)` must return (x, lab, usable) arrays.
"""
import numpy as np
from scipy.stats import spearmanr

PHASE_A, PHASE_B = "building", "dissipating"


def compact(star_ids, keep):
    """Indices of kept sectors of stars with >= 2 kept sectors, sorted by star; group ids."""
    ids = np.where(keep)[0]
    sid = star_ids[ids]
    u, inv, cnt = np.unique(sid, return_inverse=True, return_counts=True)
    ok = cnt[inv] >= 2
    ids, sid = ids[ok], sid[ok]
    o = np.argsort(sid, kind="stable")
    ids, sid = ids[o], sid[o]
    _, gid = np.unique(sid, return_inverse=True)
    return ids, gid


def demean(gid, v):
    n = np.bincount(gid)
    return v - (np.bincount(gid, weights=v) / n)[gid]


def stats(gid, x, y, lab):
    xd, yd = demean(gid, x), demean(gid, y)
    sxx = np.sum(xd * xd)
    beta = np.sum(xd * yd) / sxx if sxx > 0 else np.nan
    a, b = yd[lab == PHASE_A], yd[lab == PHASE_B]
    contrast = a.mean() - b.mean() if len(a) and len(b) else np.nan
    return beta, contrast, xd, yd


def pval(obs, null):
    """Two-sided, about the null's own centre (so a shifted null is handled)."""
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or len(null) == 0:
        return np.nan
    c = np.mean(null)
    return (1 + np.sum(np.abs(null - c) >= np.abs(obs - c))) / (1 + len(null))


def perm_null(gid, x, y, lab, n, rng):
    """Shuffle x and labels among each star's own sectors (gid must be sorted)."""
    bs, cs = np.empty(n), np.empty(n)
    for k in range(n):
        o = np.lexsort((rng.random(len(gid)), gid))
        bs[k], cs[k], _, _ = stats(gid, x[o], y, lab[o])
    return bs, cs


def bootstrap(gid, x, y, lab, n, rng):
    G = gid.max() + 1
    members = np.split(np.arange(len(gid)), np.where(np.diff(gid) != 0)[0] + 1)
    bs = np.empty(n)
    for k in range(n):
        pick = rng.integers(0, G, size=G)
        idx = np.concatenate([members[p] for p in pick])
        g2 = np.repeat(np.arange(G), [len(members[p]) for p in pick])
        bs[k], _, _, _ = stats(g2, x[idx], y[idx], lab[idx])
    return np.nanpercentile(bs, [16, 84])


def swap_null(star_ids, t, y, stars, disks, disk_of_star, evaluate, n, rng):
    """Give each star the WISE disk curve of another random star, at its own sector times."""
    rows = {s: np.where(star_ids == s)[0] for s in stars}
    tab = {}
    for s in stars:
        for j, d in enumerate(disks):
            tab[(s, j)] = evaluate(d, t[rows[s]])
    own = np.array([disk_of_star[s] for s in stars])
    N = len(disks)
    bs, cs = np.empty(n), np.empty(n)
    X = np.full(len(y), np.nan); L = np.empty(len(y), dtype=object); U = np.zeros(len(y), bool)
    for k in range(n):
        donors = rng.integers(0, N - 1, size=len(stars))
        donors = donors + (donors >= own)
        for s, j in zip(stars, donors):
            xx, ll, uu = tab[(s, j)]
            r = rows[s]
            X[r], L[r], U[r] = xx, ll, uu
        ids, gid = compact(star_ids, U & np.isfinite(y))
        if len(ids) == 0:
            bs[k] = cs[k] = np.nan
            continue
        bs[k], cs[k], _, _ = stats(gid, X[ids], y[ids], L[ids])
    return bs, cs


def run_test(star_ids, x, y, lab, usable, t, disks, disk_of_star, evaluate,
             n_null=10000, n_boot=2000, seed=20261104, do_swap=True, keep_nulls=False):
    rng = np.random.default_rng(seed)
    lab = np.asarray(lab, dtype=object)
    ids, gid = compact(star_ids, usable & np.isfinite(y) & np.isfinite(x))
    out = dict(n_stars=int(gid.max() + 1) if len(gid) else 0, n_sectors=len(ids))
    if out["n_stars"] < 3:
        return out
    xs, ys, ls = x[ids], y[ids], lab[ids]
    beta, contrast, xd, yd = stats(gid, xs, ys, ls)
    pb, pc = perm_null(gid, xs, ys, ls, n_null, rng)
    lo, hi = bootstrap(gid, xs, ys, ls, n_boot, rng)
    out.update(beta=beta, beta_lo=lo, beta_hi=hi, rho=spearmanr(xd, yd).correlation,
               contrast=contrast, n_build=int(np.sum(ls == PHASE_A)),
               n_diss=int(np.sum(ls == PHASE_B)),
               p_perm=pval(beta, pb), p_perm_contrast=pval(contrast, pc))
    if keep_nulls:
        out["_perm_null"], out["_xd"], out["_yd"] = pb, xd, yd
    if do_swap:
        stars = np.unique(star_ids[ids])
        sel = np.where(np.isin(star_ids, stars) & np.isfinite(y))[0]
        sb, sc = swap_null(star_ids[sel], t[sel], y[sel], stars, disks, disk_of_star,
                           evaluate, n_null, rng)
        out.update(p_swap=pval(beta, sb), p_swap_contrast=pval(contrast, sc),
                   swap_null_mean=float(np.nanmean(sb)))
        if keep_nulls:
            out["_swap_null"] = sb
    return out
