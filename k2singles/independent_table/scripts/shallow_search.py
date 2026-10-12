"""Look for shallower dips (possible secondary eclipses) at any depth: grid points with
SNR > 7 in >= 2 independent-ish products (K2: PDCSAP/EVEREST/K2SFF; TESS: SPOC/TESS-SPOC/
QLP SAP/own), excluding +-1 duration around known primary events. Reports phase too."""
import sys, os, pickle, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, lib
SCAN = os.path.join(lib.ROOT, "data", "scan")
OK = {"K2pdc", "K2ev", "K2sff", "spoc", "tspoc", "qlp", "own"}
D = json.load(open(os.path.join(lib.ROOT, "data", "final_details.json")))

def run(e, P=None, T0=None, snrmin=7):
    S = pickle.load(open(os.path.join(SCAN, f"{e}.pkl"), "rb")); dur = S["row"]["dur"] / 24
    ev = [x["t0"] for x in D[str(e)]["events"] if "t0" in x]
    hits = []
    for x in S["scan"]:
        if x["kind"] not in OK: continue
        snr = x["depth"] / x["err"]
        m = (snr > snrmin) & (x["cov"] >= 0.7) & (x["nb"] >= 3) & (x["na"] >= 3)
        for i in np.where(m)[0]:
            t = x["t0"][i]
            if min(abs(t - np.array(ev))) < 1.5 * dur: continue
            hits.append((t, x["src"], x["depth"][i], snr[i]))
    hits.sort()
    # cluster within dur and require >= 2 products
    out = []; i = 0
    while i < len(hits):
        j = i
        while j + 1 < len(hits) and hits[j + 1][0] - hits[i][0] < dur: j += 1
        grp = hits[i:j + 1]; srcs = set(g[1] for g in grp)
        if len(srcs) >= 2:
            b = max(grp, key=lambda g: g[3])
            ph = ((b[0] - T0) / P) % 1 if P else np.nan
            out.append((round(b[0], 3), round(ph, 3), len(srcs), round(b[2] * 1e6), round(b[3], 1), sorted(srcs)))
        i = j + 1
    return out

if __name__ == "__main__":
    for e, v in D.items():
        P = max(a[0] for a in v["allowed"]) if v.get("allowed") else None
        T0 = v["events"][0]["t0"]
        o = run(int(e), P, T0)
        print(e, f"P={P}", len(o), "multi-product shallow dips")
        for r in o: print("   ", r)
