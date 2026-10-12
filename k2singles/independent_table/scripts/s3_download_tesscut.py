"""Download TESScut 11x11 cutouts for every sector covering each target."""
import sys, os, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
from targets import T
import lightkurve as lk
from astropy.coordinates import SkyCoord
D = os.path.join(os.path.dirname(__file__), "..", "data", "tesscut")
os.makedirs(D, exist_ok=True)
for _, r in T.iterrows():
    e = int(r.epic)
    sr = lk.search_tesscut(SkyCoord(r.ra, r.dec, unit="deg"))
    for i in range(len(sr)):
        sec = str(sr.table["mission"][i]).split()[-1]
        import glob
        if glob.glob(os.path.join(D, "tesscut", f"tess-s{int(sec):04d}-*_{r.ra:.6f}_{r.dec:.6f}_11x11_astrocut.fits")): 
            print(e, sec, "exists"); continue
        try:
            tpf = sr[i].download(cutout_size=11, download_dir=D)
            print(e, sec, "ok", flush=True)
        except Exception as ex:
            print(e, sec, "fail", ex, flush=True)
