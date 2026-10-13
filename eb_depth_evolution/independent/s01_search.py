"""Step 1: list available data products for EPIC 202126877 / TIC 86491361."""
import warnings
warnings.filterwarnings("ignore")
import lightkurve as lk

out = []
for target, kw in [("EPIC 202126877", dict(mission="K2")),
                   ("TIC 86491361", dict(mission="TESS"))]:
    sr = lk.search_lightcurve(target, **kw)
    print(target, "light curves:")
    print(sr)
    out.append(sr)
sr = lk.search_targetpixelfile("EPIC 202126877", mission="K2")
print("K2 TPF:"); print(sr)
sr = lk.search_tesscut("TIC 86491361")
print("TESScut:"); print(sr)
