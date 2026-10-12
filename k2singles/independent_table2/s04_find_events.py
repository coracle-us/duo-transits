"""Blind box search for dips in every source of every star.
Detrend: iterative running median (window W days) with in-dip points masked.
Box search: for durations D, slide box over time grid, SNR = depth / (sigma/sqrt(n_in)).
Report local SNR maxima >= 7 (merged within 1 d), with depth and number of points.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib2

CFG = {  # epic: (catalogue duration h, detrend window d, box durations h)
    '246194159': (8.83, 3.0, [3, 6, 9, 14, 20, 28]),
    '212012030': (4.90, 2.0, [2, 3.5, 5, 8, 12]),
    '212694013': (40.7, 6.0, [12, 20, 30, 40, 55]),
}
SKIP = ('SAPallq', 'K2-c06-K2SC', 'K2-c05-K2SC')


def robust_detrend(t, f, W):
    mask = np.zeros_like(f, bool)
    for it in range(4):
        fd, tr = lib2.detrend(t, f, W, mask)
        r = fd - 1
        s = 1.4826 * np.median(np.abs(r[~mask] - np.median(r[~mask])))
        new = r < -3 * s
        # grow mask by a few points
        new = np.convolve(new.astype(float), np.ones(5), 'same') > 0
        if np.array_equal(new, mask):
            break
        mask = new
    return fd, s


def box_search(t, fd, sig, Dh):
    D = Dh / 24
    grid = np.arange(t.min(), t.max(), min(D / 6, 0.05))
    lo = np.searchsorted(t, grid - D / 2)
    hi = np.searchsorted(t, grid + D / 2)
    cs = np.concatenate([[0], np.cumsum(fd - 1)])
    n = hi - lo
    dep = np.where(n > 0, -(cs[hi] - cs[lo]) / np.maximum(n, 1), 0)
    # expected number of points if fully covered
    cad = np.median(np.diff(t))
    full = n >= 0.5 * D / cad
    snr = np.where(full, dep / (sig / np.sqrt(np.maximum(n, 1))), 0)
    return grid, snr, dep, n


def run(epic):
    cat_d, W, Ds = CFG[epic]
    segs = lib2.load_all(epic)
    rows = []
    for s in segs:
        if any(k in s['src'] for k in SKIP):
            continue
        t, f = s['t'], s['f']
        # split at gaps > 1 d so detrending does not bridge orbits
        fd, sig = robust_detrend(t, f, W)
        # use binned-scale sigma: estimate red noise by comparing scatter of 1-h bins
        best = {}
        for Dh in Ds:
            g, snr, dep, n = box_search(t, fd, sig, Dh)
            # white-noise SNR is optimistic; compute red-noise factor at this scale
            D = Dh / 24
            nb = int((t.max() - t.min()) / D)
            if nb > 10:
                idx = np.floor((t - t.min()) / D).astype(int)
                means = np.array([np.mean(fd[idx == k] - 1) for k in range(nb) if np.sum(idx == k) > 3])
                npts = np.median(np.bincount(idx)[np.bincount(idx) > 3])
                beta = max(1.0, 1.4826 * np.median(np.abs(means - np.median(means))) / (sig / np.sqrt(npts)))
            else:
                beta = 1.0
            snr = snr / beta
            for i in np.where(snr >= 7)[0]:
                # local max over +-0.5 d
                win = np.abs(g - g[i]) < 0.5
                if snr[i] < snr[win].max():
                    continue
                key = round(g[i], 0)
                k2 = [k for k in best if abs(k - g[i]) < 1.0]
                if k2 and best[k2[0]][1] >= snr[i]:
                    continue
                for k in k2:
                    best.pop(k)
                best[g[i]] = (Dh, snr[i], dep[i], n[i], beta)
        for tc, (Dh, sn, dp, nn, beta) in sorted(best.items()):
            rows.append((s['src'], tc, Dh, sn, dp, nn, beta))
            print(f"{epic} {s['src']:26s} t={tc:9.3f} D={Dh:5.1f}h SNR={sn:7.1f} depth={dp*1e6:9.0f}ppm n={nn:4d} beta={beta:.2f}")
    return rows


if __name__ == '__main__':
    for epic in (sys.argv[1:] or CFG):
        run(epic)
