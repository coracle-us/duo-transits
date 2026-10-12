"""Compact per-epoch summary for a trial ephemeris: best (smallest-error) exclusion-grade
product measurement at each predicted epoch that has any data."""
import sys, os, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from epochs import epochs
GOOD = {"K2pdc", "K2ev", "K2sff", "spoc", "tspoc", "qlp", "qlpk", "tglc", "own", "eleanor"}

def summary(e, T0, P, phase=0.0, show=True):
    out = epochs(e, T0, P, phase, quiet=True)
    rows = {}
    for k, tc, src, d, s, cov, ex, det in out:
        rows.setdefault(k, []).append((tc, src, d, s, cov, ex, det))
    res = []
    for k, L in sorted(rows.items()):
        tc = L[0][0]
        Lg = [x for x in L if np.isfinite(x[3]) and any(g in x[1] for g in ["PDCSAP", "EVEREST", "K2SFF", "QLP", "own", "TGLC", "eleanor", "SPOC"])]
        if not Lg:
            res.append((k, tc, "gap", np.nan, np.nan, 0, False, False)); 
            if show: print(f"  E={k:5d} t={tc:10.3f}  no usable data (gap)")
            continue
        b = min(Lg, key=lambda x: x[3])
        anyex = any(x[5] for x in L); anydet = any(x[6] for x in L)
        res.append((k, tc, b[1], b[2], b[3], b[4], anyex, anydet))
        if show:
            print(f"  E={k:5d} t={tc:10.3f}  best={b[1]:34s} depth={b[2]*1e6:8.0f}+-{b[3]*1e6:6.0f} cov={b[4]:.2f}  anyExcl={anyex} anyDet={anydet}")
    return res

if __name__ == "__main__":
    e, T0, P = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    ph = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    summary(e, T0, P, ph)
