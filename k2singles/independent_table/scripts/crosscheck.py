"""Common-mode check: measure the same box (time, duration) in the light curves of the other
C5/C16/C18-field stars (same TESS sectors). A real event should appear only in its own star."""
import sys, os, pickle, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, lib
STARS = [211892898, 211503363, 211390677, 211411112, 211489484, 212152316]
def check(tc, dur, kinds=("qlp", "own", "tspoc", "spoc")):
    for e in STARS:
        D = pickle.load(open(os.path.join(lib.ROOT, "data", "scan", f"{e}.pkl"), "rb"))
        out = []
        for s in D["segs"]:
            if s["kind"] not in kinds or not (s["t"][0] < tc < s["t"][-1]): continue
            r = lib.local_fit(s["t"], s["f"], tc, dur)
            if np.isfinite(r["depth"]): out.append(f"{s['src'].split(' ',2)[2]}: {r['depth']*1e6:.0f}+-{r['err']*1e6:.0f}")
        print(f"   {e}: " + "; ".join(out))
for tc, dur, lab in [(2532.72, 0.96, "211503363 S45"), (3270.14, 0.96, "211503363 S72"), (2547.82, 1.62, "211390677 S45"), (3284.6, 1.0, "211390677 S72 partial (box placed on covered part)")]:
    print(lab, tc); check(tc, dur)
