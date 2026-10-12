"""Download all light curves for EPIC 211924561 / TIC 238606535 into ./data."""
import os, sys, warnings
import lightkurve as lk
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data")
os.makedirs(OUT, exist_ok=True)

def grab(target, **kw):
    sr = lk.search_lightcurve(target, **kw)
    print(f"--- {target} {kw}: {len(sr)} products")
    print(sr.table["mission", "author", "exptime", "target_name", "distance"] if len(sr) else "")
    return sr

k2 = grab("EPIC 211924561", mission="K2")
tess = grab("TIC 238606535", mission="TESS")
tess_rad = grab("131.0849 20.3497", radius=20, mission="TESS")

for name, sr in [("k2", k2), ("tess", tess)]:
    for i in range(len(sr)):
        row = sr[i]
        try:
            row.download(download_dir=OUT, quality_bitmask="none")
            print("ok", name, row.table["mission"][0], row.table["author"][0], row.table["exptime"][0])
        except Exception as e:
            print("FAIL", name, i, row.table["mission"][0], row.table["author"][0], e)
