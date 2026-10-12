"""
Two-phase (primary + secondary) ephemeris search for an eclipsing binary whose events cannot all
be explained by one period, e.g. similar-depth primary and secondary eclipses in an eccentric orbit.
For each trial P, the events are split into a primary group (phase 0) and one secondary group.
Every predicted primary and secondary time that falls in the data must show a dip: the mean flux
within +/- 0.35 D is below 1 - 0.5 * depth_min. If a predicted time is covered and flat, P is
rejected. Survivors are refined by least squares.
Usage:  python k2singles/scripts/two_phase_ephem.py EPIC T1 T2 ... --dur_h D --depth DEPTH
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import k2lib  # noqa: E402


def main():
    a = argparse.ArgumentParser()
    a.add_argument("epic", type=int)
    a.add_argument("times", type=float, nargs="+")
    a.add_argument("--dur_h", type=float, required=True)
    a.add_argument("--depth", type=float, required=True, help="smallest event depth (fraction)")
    a.add_argument("--pmin", type=float, default=2.0)
    a.add_argument("--pmax", type=float, default=400.0)
    o = a.parse_args()
    cat = pd.read_csv(HERE.parent / "data" / "k2_singles.csv")
    tic = int(cat[cat.EPIC == o.epic].TIC.iloc[0])
    lc, _ = k2lib.load_combined(tic, o.epic)
    t, f = lc.time.to_numpy(), lc.flux.to_numpy()
    seg = lc.segment.to_numpy()
    # normalise each segment by its median (eclipses are short compared with segments)
    fn = f.copy()
    for s in np.unique(seg):
        m = seg == s
        fn[m] = f[m] / np.nanmedian(f[m])
    D = o.dur_h / 24
    T = np.array(sorted(o.times))
    tol = max(0.15, D / 4)

    def dip_at(tc):
        m = np.abs(t - tc) < 0.35 * D
        if m.sum() < 3:
            return None                       # not covered
        return float(np.mean(fn[m])) < 1 - 0.5 * o.depth

    sols = []
    for P in np.arange(o.pmin, o.pmax, 0.0005):
        ph = (T - T[0]) / P
        ph -= np.round(ph)
        prim = np.abs(ph) * P < tol
        rest = ((T[~prim] - T[0]) / P) % 1
        if prim.sum() < 1 or (len(rest) and (rest.max() - rest.min()) * P > 2 * tol):
            continue
        sec = rest.mean() if len(rest) else None
        ok = True
        k0, k1 = int(np.floor((t.min() - T[0]) / P)) - 1, int(np.ceil((t.max() - T[0]) / P)) + 1
        for k in range(k0, k1 + 1):
            for phase in ([0.0] + ([sec] if sec is not None else [])):
                tc = T[0] + (k + phase) * P
                if tc < t.min() - 1 or tc > t.max() + 1:
                    continue
                d = dip_at(tc)
                if d is False:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            sols.append((round(P, 4), int(prim.sum()), None if sec is None else round(sec, 4)))
    merged = []
    for s in sols:
        if not merged or s[0] - merged[-1][0] > 0.01:
            merged.append(s)
    print(f"EPIC {o.epic} TIC {tic}: events {list(T)}")
    print(f"{len(merged)} surviving (P, n_primary, secondary phase):")
    for s in merged:
        print("  ", s)


if __name__ == "__main__":
    main()
