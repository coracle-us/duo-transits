"""Step 2: download K2 + TESS light curves, K2 TPF, and TESScut cutouts into data/."""
import os, warnings
warnings.filterwarnings("ignore")
import lightkurve as lk

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(D, exist_ok=True)

sr = lk.search_lightcurve("EPIC 202126877", mission="K2")
for i in range(len(sr)):
    try:
        sr[i].download(download_dir=D, flux_column=None) if False else sr[i].download(download_dir=D)
        print("ok K2", sr[i].author[0])
    except Exception as e:
        print("fail K2", sr[i].author[0], e)

lk.search_targetpixelfile("EPIC 202126877", mission="K2").download(download_dir=D)
print("ok K2 TPF")

sr = lk.search_lightcurve("TIC 86491361", mission="TESS")
for i in range(len(sr)):
    try:
        sr[i].download(download_dir=D)
        print("ok TESS", sr[i].mission[0], sr[i].author[0])
    except Exception as e:
        print("fail TESS", sr[i].mission[0], sr[i].author[0], e)

# TESScut via astroquery (handles odd sector labels)
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
c = SkyCoord(97.98012, 25.48207, unit="deg")
sec = Tesscut.get_sectors(coordinates=c)
print(sec)
paths = Tesscut.download_cutouts(coordinates=c, size=25, path=os.path.join(D, "tesscut"))
print(paths)
