"""Per star: load every light curve, scan each with a local-polynomial box search at the
catalogue duration, record detections and per-time 'covered & flat' exclusion flags."""
import sys, os, pickle, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from targets import T
import lib

OUT = os.path.join(lib.ROOT, "data", "scan")
os.makedirs(OUT, exist_ok=True)
EXCL_KINDS = {"K2pdc", "K2ev", "K2sff", "spoc", "tspoc", "qlp", "qlpk", "tglc", "own", "eleanor"}
OTHER_IDS = {211411112: [200233362]}


def load_all(r):
    segs = lib.load_k2(int(r.epic), OTHER_IDS.get(int(r.epic), []))
    segs += lib.load_tess_pipeline(int(r.tic))
    for fn in lib.tesscut_files(r.ra, r.dec):
        try:
            s = lib.own_tesscut_phot(fn)
            if s: segs.append(s)
        except Exception as ex:
            print("  own phot fail", os.path.basename(fn), ex)
    return segs


def split_gaps(seg, gap=1.0):
    """Split a segment at gaps > gap days (orbits), so the local fit never spans a gap."""
    t = seg["t"]; br = np.where(np.diff(t) > gap)[0]
    edges = [0] + list(br + 1) + [len(t)]
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a < 30: continue
        s = dict(seg); s["t"] = seg["t"][a:b]; s["f"] = seg["f"][a:b]
        out.append(s)
    return out


def run(epic_only=None):
    for _, r in T.iterrows():
        e = int(r.epic)
        if epic_only and e not in epic_only: continue
        dur = r.dur / 24.0; dep = r.depth * 1e-6
        segs = load_all(r)
        print(e, "segments:", len(segs), flush=True)
        res = []
        for s in segs:
            for piece in split_gaps(s, gap=max(0.5, dur)):
                rows = lib.scan(piece, dur, step=dur / 8)
                if not rows: continue
                A = {k: np.array([x.get(k, np.nan) for x in rows], float) for k in ["t0", "depth", "err", "cov", "cov_full", "nb", "na"]}
                is_tess = not s["kind"].startswith("K2")
                deff = dep * (0.5 if is_tess else 1.0)
                snr = A["depth"] / A["err"]
                A["det"] = (snr > 5) & (A["depth"] > 0.3 * deff) & (A["cov"] >= 0.5)
                A["excl"] = (A["cov"] >= 0.5) & ((deff - A["depth"]) / A["err"] > 3) & (A["depth"] < 0.5 * deff) & (s["kind"] in EXCL_KINDS)
                res.append(dict(src=s["src"], kind=s["kind"], camp=s["camp"], tmin=piece["t"][0], tmax=piece["t"][-1],
                                cad=lib.cadence(piece["t"]), rms=lib.robust_std(np.diff(piece["f"])) / np.sqrt(2), **A))
        pickle.dump(dict(row=r.to_dict(), segs=segs, scan=res), open(os.path.join(OUT, f"{e}.pkl"), "wb"))
        # summary of detections
        for x in res:
            m = x["det"]
            if not m.any(): continue
            t0 = x["t0"][m]; snr = (x["depth"] / x["err"])[m]; d = x["depth"][m]
            # cluster
            groups = np.split(np.arange(len(t0)), np.where(np.diff(t0) > dur)[0] + 1)
            for g in groups:
                i = g[np.argmax(snr[g])]
                print(f"   DET {x['src']:38s} t0={t0[i]:10.3f} depth={d[i]*1e6:8.0f} ppm snr={snr[i]:6.1f}")
        print(flush=True)


if __name__ == "__main__":
    run([int(a) for a in sys.argv[1:]] or None)
