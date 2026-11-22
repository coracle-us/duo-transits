"""Download TESS-SPOC and QLP light curves, and TESScut cutouts for sectors that
have no LC or are >= S94 (QLP clock-offset caveat). Also fetch TIC stellar params."""
import warnings, os, json; warnings.filterwarnings("ignore")
import lightkurve as lk
from astroquery.mast import Catalogs
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data")
TICS = [275647095, 241550976, 107113048, 230133616]
stell = {}
for t in TICS:
    sr = lk.search_lightcurve(f"TIC {t}", mission="TESS", author=["TESS-SPOC", "QLP"])
    sr.download_all(download_dir=os.path.join(D, "lk"), quality_bitmask="none")
    lcsec = set(int(m.split()[-1]) for m in sr.table["mission"])
    sc = lk.search_tesscut(f"TIC {t}")
    want = [i for i, m in enumerate(sc.table["mission"])
            if int(m.split()[-1]) not in lcsec or int(m.split()[-1]) >= 94]
    if want:
        print(t, "TESScut sectors:", [sc.table["mission"][i] for i in want])
        sc[want].download_all(cutout_size=15, download_dir=os.path.join(D, "tesscut"))
    r = Catalogs.query_object(f"TIC {t}", catalog="TIC", radius=0.001)[0]
    stell[t] = {k: (float(r[k]) if r[k] is not None and str(r[k]) != "--" else None)
                for k in ["ra", "dec", "Tmag", "Teff", "rad", "mass", "logg", "rho", "contratio", "d"]}
    stell[t]["GAIA"] = str(r["GAIA"])
    print(t, stell[t])
json.dump(stell, open(os.path.join(ROOT, "out", "stellar.json"), "w"), indent=1)
