"""
Window extraction and features for the blend-ceiling classifier (RULES.md section 2 and amendment A).

Preprocessing follows the description in RULES.md section 2 (our own implementation; no code is copied from
helemysm/singletransiters): quality == 0, NaNs removed, 30 min bins, a 99.5 h window of 200 points centred on the
predicted time, linear interpolation with the light-curve median beyond the data, divided by the window median,
minus 1. Grid points with no binned data within 15 min are flagged 'filled'; features ignore them.
"""
import numpy as np

BIN_H = 0.5
WIN_H = 99.5
L = int(np.floor(WIN_H / BIN_H)) + 1  # 200
T_REL_H = -WIN_H / 2 + BIN_H * np.arange(L)
DURS_H = np.array([1.5, 2, 3, 4, 6, 8, 10, 12, 15])
CENTRES_H = np.arange(-12, 12.01, 0.5)
FEATURES = ["snr", "depth", "log_depth", "dur_h", "abs_c_h", "sig_out", "sig_p2p", "red_ratio", "v_u", "asym",
            "snr2", "snr2_ratio", "skew", "max_up", "frac_filled"]
DEPTH_SCALE = ["snr", "depth", "log_depth"]  # dropped in the shape-only model (amendment A3)


def mad_std(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan


def prep(lc):
    """lightkurve LightCurve -> binned (t [BTJD], f) with quality == 0 and finite flux."""
    if "quality" in lc.colnames:
        lc = lc[np.asarray(lc.quality) == 0]
    lc = lc.remove_nans()
    if len(lc) == 0:
        return np.array([]), np.array([])
    b = lc.bin(time_bin_size=BIN_H / 24).remove_nans()
    t = np.asarray(b.time.value, float)
    f = np.asarray(b.flux.value, float)
    ok = np.isfinite(t) & np.isfinite(f)
    return t[ok], f[ok]


def window(t, f, T):
    """Return (grid flux minus 1, filled mask) for a window centred on T (BTJD)."""
    tg = T + T_REL_H / 24
    med = np.nanmedian(f)
    fg = np.interp(tg, t, f, left=med, right=med)
    m2 = np.nanmedian(fg)
    if not np.isfinite(m2) or m2 == 0:
        return None, None
    idx = np.searchsorted(t, tg)
    near = np.full(L, np.inf)
    for j in (idx - 1, idx):
        jj = np.clip(j, 0, len(t) - 1)
        near = np.minimum(near, np.abs(t[jj] - tg) * 24)
    return fg / m2 - 1.0, near > BIN_H / 2


def valid(t, T, t14_h):
    """RULES.md section 2: >= 50% of bins within +-T14/2 and >= 50% of the 200 window bins hold data."""
    if len(t) == 0 or not np.isfinite(t14_h) or t14_h <= 0:
        return False
    dt_h = (t - T) * 24
    n_in = np.sum(np.abs(dt_h) < t14_h / 2)
    n_win = np.sum(np.abs(dt_h) <= WIN_H / 2)
    return n_in >= 0.5 * max(t14_h / BIN_H, 1) and n_win >= 0.5 * L


def _box(x, ok, c, D):
    d = T_REL_H - c
    inn = ok & (np.abs(d) < D / 2)
    out = ok & (np.abs(d) > D)
    if inn.sum() < 2 or out.sum() < 10:
        return np.nan, np.nan, np.nan
    so = mad_std(x[out])
    depth = np.median(x[out]) - np.mean(x[inn])
    return depth / (so / np.sqrt(inn.sum())) if so > 0 else np.nan, depth, so


def features(x, filled):
    ok = ~filled
    best = (-np.inf, None, None, None, None)
    for D in DURS_H:
        for c in CENTRES_H:
            snr, dep, so = _box(x, ok, c, D)
            if np.isfinite(snr) and snr > best[0]:
                best = (snr, dep, so, D, c)
    snr, dep, so, D, c = best
    if D is None:
        return None
    d = T_REL_H - c
    inn = ok & (np.abs(d) < D / 2)
    inner = ok & (np.abs(d) < D / 4)
    base = np.median(x[ok & (np.abs(d) > D)])
    full_mean = np.mean(x[inn]) - base
    v_u = (np.mean(x[inner]) - base) / full_mean if inner.sum() and full_mean != 0 else np.nan
    left, right = inn & (d < 0), inn & (d >= 0)
    asym = ((np.mean(x[left]) - np.mean(x[right])) / dep if left.sum() and right.sum() and dep != 0 else np.nan)
    # second dip: best box at least D away from the first (centres over the whole window)
    snr2 = -np.inf
    for D2 in DURS_H:
        for c2 in np.arange(-WIN_H / 2 + D2, WIN_H / 2 - D2 + 0.01, 0.5):
            if abs(c2 - c) < D + D2 / 2:
                continue
            s2, _, _ = _box(x, ok & (np.abs(d) > D), c2, D2)
            if np.isfinite(s2):
                snr2 = max(snr2, s2)
    xo = x[ok]
    p2p = mad_std(np.diff(xo)) / np.sqrt(2) if ok.sum() > 3 else np.nan
    sk = np.mean(((xo - xo.mean()) / xo.std()) ** 3) if xo.std() > 0 else np.nan
    return dict(snr=snr, depth=dep, log_depth=np.log10(max(dep, 1e-7)), dur_h=D, abs_c_h=abs(c), sig_out=so,
                sig_p2p=p2p, red_ratio=so / p2p if p2p and p2p > 0 else np.nan, v_u=v_u, asym=asym,
                snr2=snr2 if np.isfinite(snr2) else np.nan,
                snr2_ratio=snr2 / snr if np.isfinite(snr2) and snr > 0 else np.nan,
                skew=sk, max_up=(np.max(xo) - np.median(xo)) / so if so > 0 else np.nan,
                frac_filled=float(filled.mean()))


def sector_box_search(t, f, exclude_T, min_sep_d=3.0, step_h=0.5):
    """TCE-like windows (amendment A2): highest-SNR box in one light curve, centre >= min_sep_d from every
    excluded time. Uses the same in/out definition as features(), on a window centred at each trial centre."""
    if len(t) < L // 2:
        return None
    cs = np.arange(t.min() + 1.0, t.max() - 1.0, step_h / 24)
    if len(exclude_T):
        ex = np.asarray(exclude_T)
        cs = cs[np.min(np.abs(cs[:, None] - ex[None, :]), axis=1) >= min_sep_d]
    best = (-np.inf, None, None)
    med = np.nanmedian(f)
    for c in cs:
        dt_h = (t - c) * 24
        w = np.abs(dt_h) <= WIN_H / 2
        if w.sum() < 0.5 * L:
            continue
        x, d = f[w] / med - 1, dt_h[w]
        for D in DURS_H:
            inn, out = np.abs(d) < D / 2, np.abs(d) > D
            if inn.sum() < 0.5 * D / BIN_H or out.sum() < 10:
                continue
            so = mad_std(x[out])
            if so <= 0:
                continue
            snr = (np.median(x[out]) - np.mean(x[inn])) / (so / np.sqrt(inn.sum()))
            if snr > best[0]:
                best = (snr, c, D)
    return None if best[1] is None else dict(snr=best[0], t_c=best[1], dur_h=best[2])
