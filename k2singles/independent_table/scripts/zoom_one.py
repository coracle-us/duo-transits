"""Zoom plot of all products around a time: python zoom_one.py EPIC TIME [halfwin] [tag]"""
import sys, os, pickle
sys.path.insert(0, os.path.dirname(__file__))
from s6_plots import zoom
import lib
e = int(sys.argv[1]); tc = float(sys.argv[2])
hw = float(sys.argv[3]) if len(sys.argv) > 3 else None
tag = sys.argv[4] if len(sys.argv) > 4 else "check"
P = pickle.load(open(os.path.join(lib.ROOT, "data", "scan", f"{e}.pkl"), "rb")); r = P["row"]
segs = [s for s in P["segs"] if s["kind"] not in ("t16", "cdips")]
zoom(e, segs, tc, r["dur"] / 24, r["depth"] * 1e-6, f"{e}_{tag}_{tc:.2f}.png", halfwin=hw)
print("ok")
