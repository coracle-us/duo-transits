"""Pure-numpy pieces of the second-transit search (no network, easy to test).

Times are BTJD (BJD - 2457000) in days. Flux is detrended and normalised to 1.
"""
import numpy as np

G = 6.674e-11
RSUN = 6.957e8
MSUN = 1.989e30
DAY = 86400.0


def robust_std(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 3:
        return np.nan
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def cadence_days(t):
    dt = np.diff(np.sort(t))
    dt = dt[dt > 0]
    return float(np.median(dt)) if dt.size else np.nan


def duration_scale_noise(t, f, dur_d):
    """Scatter of non-overlapping bin means at the transit duration.

    Captures red noise, which a per-point sigma/sqrt(N) misses.
    """
    cad = cadence_days(t)
    n_exp = dur_d / cad
    edges = np.arange(t.min(), t.max() + dur_d, dur_d)
    idx = np.digitize(t, edges)
    means = []
    for k in np.unique(idx):
        m = idx == k
        if m.sum() >= 0.5 * n_exp:
            means.append(f[m].mean())
    sig_bins = robust_std(np.array(means)) if len(means) >= 5 else np.nan
    sig_white = robust_std(f) / np.sqrt(max(n_exp, 1.0))
    if not np.isfinite(sig_bins):
        return sig_white
    return max(sig_bins, sig_white)


def box_stats(t, f, centers, dur_d):
    """Mean in-box depth and fill fraction for each trial centre."""
    order = np.argsort(t)
    t, f = t[order], f[order]
    cs = np.concatenate([[0.0], np.cumsum(f)])
    lo = np.searchsorted(t, centers - dur_d / 2)
    hi = np.searchsorted(t, centers + dur_d / 2)
    n = hi - lo
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = (cs[hi] - cs[lo]) / n
    depth = 1.0 - mean
    frac = n / (dur_d / cadence_days(t))
    return depth, frac, n


def box_search(t, f, dur_d, step_frac=0.125, min_frac=0.5):
    """Slide a box of fixed width across the data.

    Returns centres, depth, snr (depth / duration-scale noise, scaled for
    partial coverage), fill fraction and the noise used.
    """
    sig = duration_scale_noise(t, f, dur_d)
    centers = np.arange(t.min() + dur_d / 2, t.max() - dur_d / 2, dur_d * step_frac)
    depth, frac, n = box_stats(t, f, centers, dur_d)
    with np.errstate(invalid="ignore", divide="ignore"):
        snr = depth / (sig / np.sqrt(np.clip(frac, 1e-3, None)))
    snr[(frac < min_frac) | ~np.isfinite(snr)] = np.nan
    return centers, depth, snr, frac, sig


def pick_events(centers, snr, dur_d, n_max=5, sign=1):
    """Greedy peak picking: best SNR, mask +/- one duration, repeat."""
    s = sign * snr.copy()
    out = []
    for _ in range(n_max):
        if not np.any(np.isfinite(s)):
            break
        i = int(np.nanargmax(s))
        out.append(i)
        s[np.abs(centers - centers[i]) < dur_d] = np.nan
    return out


def refine_event(t, f, tc, dur0_d, sig_fn=duration_scale_noise):
    """Grid over centre (+/- dur0/2) and duration (0.5-2 x dur0) for best box SNR."""
    best = None
    for dur in dur0_d * np.geomspace(0.5, 2.0, 13):
        sig = sig_fn(t, f, dur)
        cs = tc + np.linspace(-dur0_d / 2, dur0_d / 2, 21)
        depth, frac, n = box_stats(t, f, cs, dur)
        ok = frac >= 0.5
        if not ok.any():
            continue
        snr = np.where(ok, depth / (sig / np.sqrt(np.clip(frac, 1e-3, None))), -np.inf)
        j = int(np.argmax(snr))
        if best is None or snr[j] > best["snr"]:
            best = dict(tc=float(cs[j]), dur_d=float(dur), depth=float(depth[j]),
                        depth_err=float(sig / np.sqrt(max(frac[j], 1e-3))),
                        snr=float(snr[j]), frac=float(frac[j]))
    return best


def gap_edge_distance(t, tc, gap_d=0.5):
    """Distance (days) from tc to the nearest data edge or gap longer than gap_d."""
    ts = np.sort(t)
    dt = np.diff(ts)
    edges = [ts[0], ts[-1]]
    for i in np.where(dt > gap_d)[0]:
        edges += [ts[i], ts[i + 1]]
    return float(np.min(np.abs(np.array(edges) - tc)))


def a_over_rstar(P_d, mstar, rstar):
    a = (G * mstar * MSUN * (P_d * DAY) ** 2 / (4 * np.pi ** 2)) ** (1 / 3)
    return a / (rstar * RSUN)


def max_duration_d(P_d, mstar, rstar, k):
    """Longest possible duration (b=0, circular) for period P."""
    x = (1 + k) / a_over_rstar(P_d, mstar, rstar)
    return P_d / np.pi * np.arcsin(np.clip(x, 0, 1))


def pmin_from_duration(dur_d, mstar, rstar, depth):
    """Shortest circular-orbit period that can host a transit this long.

    Eccentric orbits can beat this, so treat it as a soft limit.
    """
    if not (np.isfinite(mstar) and np.isfinite(rstar) and mstar > 0 and rstar > 0):
        return np.nan
    k = np.sqrt(max(depth, 0))
    lo, hi = 0.05, 1e5
    if max_duration_d(hi, mstar, rstar, k) < dur_d:
        return np.nan
    for _ in range(80):
        mid = np.sqrt(lo * hi)
        if max_duration_d(mid, mstar, rstar, k) >= dur_d:
            hi = mid
        else:
            lo = mid
    return hi


def alias_table(t1, t2, dur_d, depth0, depth0_err, t, f, p_floor=1.0,
                pmin_density=np.nan, n_max=3000, min_frac=0.5):
    """Status of every alias P_n = (t2 - t1)/n against all available data.

    For each alias, look at each predicted epoch (other than the two known
    events) that falls in the data. An epoch with coverage >= min_frac whose
    box depth is below max(depth0 - 3 sigma, 0.5 depth0) excludes the alias;
    one whose depth agrees with depth0 counts as a supporting (third) transit.
    """
    dT = t2 - t1
    rows = []
    tmin, tmax = np.nanmin(t), np.nanmax(t)
    sig = duration_scale_noise(t, f, dur_d)
    for n in range(1, n_max + 1):
        P = dT / n
        if P < p_floor:
            break
        kmin = int(np.floor((tmin - t1) / P)) - 1
        kmax = int(np.ceil((tmax - t1) / P)) + 1
        ks = np.array([k for k in range(kmin, kmax + 1) if k not in (0, n)])
        status, n_cov, n_supp, worst = "allowed", 0, 0, np.nan
        if ks.size:
            ep = t1 + ks * P
            depth, frac, _ = box_stats(t, f, ep, dur_d)
            cov = frac >= min_frac
            n_cov = int(cov.sum())
            if n_cov:
                d = depth[cov]
                err = sig / np.sqrt(frac[cov])
                thresh = np.maximum(depth0 - 3 * np.hypot(err, depth0_err), 0.5 * depth0)
                excl = d < thresh
                supp = np.abs(d - depth0) < 3 * np.hypot(err, depth0_err)
                n_supp = int(supp.sum())
                worst = float(np.min(d))
                if excl.any():
                    status = "excluded"
                elif n_supp:
                    status = "supported"
        if np.isfinite(pmin_density) and P < pmin_density and status != "excluded":
            status = "allowed_below_pmin" if status == "allowed" else status
        rows.append(dict(n=n, period_d=P, status=status, n_epochs_covered=n_cov,
                         n_epochs_supporting=n_supp, min_depth_at_epochs=worst))
    return rows
