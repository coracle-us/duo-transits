#!/usr/bin/env python3
"""Saved record of the window 4 result in the RNAAS note (rnaas/rnaas.tex, Section 3, "Window 4").

Uses kepler_alias_scan.py's data and rule unchanged (same Kepler light curves from MAST's public AWS bucket, same
quiet level, same 0.8% additive-deficit test, same T0 and core). Two parts:

1. Scan 1025-1375 d in 0.01 d steps, the range of Madurga-Favieres et al. (2026)'s five windows. Surviving ranges go
   to window_scan_1025_1375.csv (kepler_alias_survivors.csv, the 770-845 d result, is left alone). The note says
   windows 1, 2, 3 and 5 are recovered and window 4 is not.
2. Detail for window 4 (P = 1345.648 +- 0.008 d, MF26): for P = 1345.640, 1345.648, 1345.656 and each predicted
   Kepler epoch, the status, the core's cadences, how many of them are dimmed by less than 0.8%, the median core
   deficit and the deepest deficit within +-1 d.

Usage:  python window4_check.py      (run from boyajian-companion/; writes window_scan_1025_1375.csv)
"""
import numpy as np

import kepler_alias_scan as kas

MF26_WINDOWS = {1: "~1030 d", 2: "~1035 d", 3: "1301.0-1303.6 d", 4: "1345.648 +- 0.008 d", 5: "~1368 d"}


def main():
    t, f = kas.kepler_lc()
    q = kas.quiet_level(t, f)
    deficit = q - f
    visible = deficit < 0.008

    def status(tc):
        i, j = np.searchsorted(t, [tc - kas.T23 / 2, tc + kas.T23 / 2])
        if visible[i:j].sum() >= 3:
            return "excluded"
        return "hidden" if j - i >= 3 else "gap"

    def epochs(P):
        return [n for n in range(-12, 0) if 120.0 < kas.T0 + n * P < 1591.5]

    print(f"Kepler cadences: {t.size}, BKJD {t.min():.2f}-{t.max():.2f}; T0 = BKJD {kas.T0:.3f}; core half-width "
          f"{kas.T23 / 2:.2f} d\n")

    # part 1: the MF26 range
    step = 0.01
    Ps = np.round(np.arange(1025.0, 1375.0 + 1e-9, step), 4)
    surv = [(P, [(n, status(kas.T0 + n * P)) for n in epochs(P)]) for P in Ps]
    surv = [(P, st) for P, st in surv if all(s != "excluded" for _, s in st)]
    rows, i = [], 0
    while i < len(surv):
        j = i
        while j + 1 < len(surv) and abs(surv[j + 1][0] - surv[j][0] - step) < 1e-6:
            j += 1
        rows.append((surv[i][0], surv[j][0], surv[i][1]))
        i = j + 1
    print("Part 1: surviving ranges, 1025-1375 d (MF26 windows: "
          + "; ".join(f"{k}: {v}" for k, v in MF26_WINDOWS.items()) + ")")
    with open("window_scan_1025_1375.csv", "w") as fh:
        fh.write("p_min,p_max,epochs\n")
        for lo, hi, st in rows:
            print(f"  {lo:.2f}-{hi:.2f}  {st}")
            fh.write(f"{lo:.2f},{hi:.2f},\"{st}\"\n")
    w4 = [r for r in rows if r[1] >= 1345.0 and r[0] <= 1346.3]
    print(f"  ranges overlapping 1345.0-1346.3 d: {len(w4)}\n")

    # part 2: window 4 detail
    print("Part 2: window 4 epochs")
    for P in (1345.640, 1345.648, 1345.656):
        for n in epochs(P):
            tc = kas.T0 + n * P
            i, j = np.searchsorted(t, [tc - kas.T23 / 2, tc + kas.T23 / 2])
            a, b = np.searchsorted(t, [tc - 1.0, tc + 1.0])
            core = deficit[i:j]
            print(f"  P={P:.3f} n={n:+d} BKJD {tc:8.2f}: {status(tc):8s} core cadences {j - i:3d}, "
                  f"dimmed <0.8%: {int(visible[i:j].sum()):3d}, median core deficit "
                  f"{(np.median(core) * 100 if core.size else np.nan):5.2f}%, deepest within +-1 d "
                  f"{(deficit[a:b].max() * 100 if b > a else np.nan):5.2f}%")


if __name__ == "__main__":
    main()
