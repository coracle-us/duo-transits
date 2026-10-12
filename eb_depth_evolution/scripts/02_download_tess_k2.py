"""Download TESScut cutouts (all sectors) and K2 C0 TPF + light curves for EPIC 202126877 / TIC 86491361."""
import os, warnings; warnings.filterwarnings("ignore")
import lightkurve as lk
out = os.path.join(os.path.dirname(__file__), "..", "data")
d = os.path.join(out, "tesscut"); os.makedirs(d, exist_ok=True)
sr = lk.search_tesscut("TIC 86491361")
for i in range(len(sr)):
    try:
        sr[i].download(cutout_size=21, download_dir=d)
        print("ok", sr[i].mission[0])
    except Exception as e: print("fail", sr[i].mission[0], e)
d2 = os.path.join(out, "k2"); os.makedirs(d2, exist_ok=True)
for kw in [dict(author="K2"), ]:
    lk.search_targetpixelfile("EPIC 202126877", **kw).download_all(download_dir=d2)
lk.search_lightcurve("EPIC 202126877").download_all(download_dir=d2)
d3 = os.path.join(out, "tess_lc"); os.makedirs(d3, exist_ok=True)
lk.search_lightcurve("TIC 86491361").download_all(download_dir=d3)
print("done")
