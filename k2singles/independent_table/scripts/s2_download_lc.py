"""Download all K2 (K2 mission, EVEREST, K2SFF) and TESS pipeline light curves."""
import sys, os, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
from targets import T
import lightkurve as lk
from astropy.coordinates import SkyCoord
D = os.path.join(os.path.dirname(__file__), "..", "data", "mast")
for _, r in T.iterrows():
    e, tic = int(r.epic), int(r.tic)
    for q, kw in [(f"EPIC {e}", dict(mission="K2")), (SkyCoord(r.ra, r.dec, unit="deg"), dict(mission="K2", radius=3)), (f"TIC {tic}", dict(mission="TESS"))]:
        try:
            sr = lk.search_lightcurve(q, **kw)
        except Exception as ex:
            print(e, "search fail", ex); continue
        for i in range(len(sr)):
            try:
                sr[i].download(download_dir=D, quality_bitmask="none")
            except Exception as ex:
                print(e, "dl fail", sr.table["mission"][i], sr.table["author"][i], ex, flush=True)
        print(e, "done", q if isinstance(q, str) else "coord", len(sr), flush=True)
