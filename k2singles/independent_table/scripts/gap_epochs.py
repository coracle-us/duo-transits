"""List predicted epochs that fall inside a TESScut sector's time span (or a K2 campaign
span) but where no pipeline/own product has usable data -> candidates for raw inspection."""
import sys, os, glob, re, pickle, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from astropy.io import fits
import lib
from targets import T
from epoch_summary import summary

def spans(e):
    r = T[T.epic == e].iloc[0]; out = []
    for fn in lib.tesscut_files(r.ra, r.dec):
        t = fits.open(fn)[1].data["TIME"]; sec = int(re.search(r"tess-s(\d{4})", fn).group(1))
        out.append((f"S{sec}", np.nanmin(t), np.nanmax(t)))
    D = pickle.load(open(os.path.join(lib.ROOT, "data", "scan", f"{e}.pkl"), "rb"))
    for s in D["segs"]:
        if s["kind"] == "K2sap": out.append((f"K2 C{s['camp']}", s["t"][0], s["t"][-1]))
    return out

if __name__ == "__main__":
    e, T0, P = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    ph = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    sp = spans(e)
    for k, tc, src, *_ in summary(e, T0, P, ph, show=False):
        if src != "gap": continue
        inn = [n for n, a, b in sp if a - 0.2 < tc < b + 0.2]
        if inn: print(f"  E={k} t={tc:.3f} in span of {inn} but no usable product data")
