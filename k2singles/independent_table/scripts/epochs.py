"""For a trial ephemeris, list every predicted epoch inside any light curve and what each
product measures there (depth, error, coverage, excluded/detected flags)."""
import sys, os, pickle, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, lib
SCAN = os.path.join(lib.ROOT, "data", "scan")

def epochs(e, T0, P, phase=0.0, quiet=False, kinds=None):
    D = pickle.load(open(os.path.join(SCAN, f"{e}.pkl"), "rb")); r = D["row"]
    dur = r["dur"] / 24; dep = r["depth"] * 1e-6
    tmin = min(s["t"][0] for s in D["segs"]); tmax = max(s["t"][-1] for s in D["segs"])
    k0 = int(np.floor((tmin - T0) / P)) - 1; k1 = int(np.ceil((tmax - T0) / P)) + 1
    out = []
    for k in range(k0, k1 + 1):
        tc = T0 + (k + phase) * P
        for x in D["scan"]:
            if kinds and x["kind"] not in kinds: continue
            if not (x["tmin"] - dur < tc < x["tmax"] + dur): continue
            i = np.argmin(np.abs(x["t0"] - tc))
            if abs(x["t0"][i] - tc) > dur / 8: 
                out.append((k, tc, x["src"], np.nan, np.nan, 0.0, False, False)); continue
            out.append((k, tc, x["src"], x["depth"][i], x["err"][i], x["cov"][i], bool(x["excl"][i]), bool(x["det"][i])))
    if not quiet:
        for o in out:
            print(f"  E={o[0]:5d} t={o[1]:10.3f} {o[2]:36s} depth={o[3]*1e6:9.0f} +- {o[4]*1e6:7.0f} cov={o[5]:.2f} excl={o[6]} det={o[7]}")
    return out

if __name__ == "__main__":
    e, T0, P = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    ph = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    epochs(e, T0, P, ph)
