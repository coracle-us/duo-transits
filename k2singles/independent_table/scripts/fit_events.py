"""Fit trapezoids to listed events in every product that covers them; print table."""
import sys, os, pickle, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, lib
SCAN = os.path.join(lib.ROOT, "data", "scan")
OK = {"K2pdc", "K2ev", "K2sff", "spoc", "tspoc", "qlp", "qlpk", "tglc", "own", "eleanor", "K2sap", "tars"}

def fit_all(e, tc, durg=None, depg=None, kinds=OK, quiet=False):
    P = pickle.load(open(os.path.join(SCAN, f"{e}.pkl"), "rb")); r = P["row"]
    dur = durg or r["dur"] / 24; dep = depg or r["depth"] * 1e-6
    rows = []
    for s in P["segs"]:
        if s["kind"] not in kinds: continue
        m = np.abs(s["t"] - tc) < dur / 2 + max(1.5 * dur, 0.75)
        if m.sum() < 15 or np.sum(np.abs(s["t"] - tc) < dur / 2) < 3: continue
        f = lib.fit_event(s["t"], s["f"], tc, dep, dur)
        if f is None: continue
        rows.append((s["src"], f))
        if not quiet:
            print(f"  {s['src']:36s} t0={f['t0']:.5f}+-{f['e_t0']:.5f}  depth={f['depth']*1e6:7.0f}+-{f['e_depth']*1e6:5.0f}  T14={f['T14']*24:6.2f}h T12={f['T12']*24:5.2f}h beta={f['beta']:.2f} n={f['n']}")
    return rows

if __name__ == "__main__":
    e = int(sys.argv[1])
    for tc in sys.argv[2:]:
        print(e, tc); fit_all(e, float(tc))
