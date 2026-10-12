"""Inventory: K2 LCs (all campaigns) and TESS products + TESScut sectors for each target."""
import sys, os, warnings, json
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
from targets import T
import lightkurve as lk
from astropy.coordinates import SkyCoord
import astropy.units as u
out = {}
for _, r in T.iterrows():
    e = int(r.epic); d = {}
    try:
        sr = lk.search_lightcurve(f"EPIC {e}", mission="K2")
        d["k2"] = [(str(x["mission"]), str(x["author"]), float(x["exptime"]), str(x["target_name"])) for x in sr.table] if len(sr) else []
    except Exception as ex:
        d["k2"] = "ERR " + str(ex)
    try:
        sr = lk.search_lightcurve(f"TIC {int(r.tic)}")
        d["tess"] = [(str(x["mission"]), str(x["author"]), float(x["exptime"])) for x in sr.table] if len(sr) else []
    except Exception as ex:
        d["tess"] = "ERR " + str(ex)
    try:
        sc = lk.search_tesscut(SkyCoord(r.ra, r.dec, unit="deg"))
        d["tesscut"] = [str(m) for m in sc.table["mission"]]
    except Exception as ex:
        d["tesscut"] = "ERR " + str(ex)
    out[e] = d
    print(e, "K2:", d["k2"] if isinstance(d["k2"], str) else len(d["k2"]), "TESS:", d["tess"] if isinstance(d["tess"], str) else len(d["tess"]), "cut:", d["tesscut"], flush=True)
json.dump(out, open(os.path.join(os.path.dirname(__file__), "..", "data", "inventory.json"), "w"), indent=1)
