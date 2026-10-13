"""Step 4: TIC parameters of target and neighbours within 2.5 arcmin (for dilution)."""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np
from astroquery.mast import Catalogs
from astropy.coordinates import SkyCoord
import astropy.units as u

H = os.path.dirname(os.path.abspath(__file__))
c = SkyCoord(97.98012, 25.48207, unit="deg")
t = Catalogs.query_region(c, radius=150 * u.arcsec, catalog="TIC")
sep = SkyCoord(t["ra"], t["dec"], unit="deg").separation(c).arcsec
t["sep"] = sep
t.sort("sep")
cols = ["ID", "ra", "dec", "sep", "Tmag", "GAIAmag", "Vmag", "Kepmag", "Teff", "logg", "rad", "mass", "d", "plx", "ruwe" if "ruwe" in t.colnames else "objType"]
t = t[[cc for cc in cols if cc in t.colnames]]
t.write(os.path.join(H, "data", "tic_neighbours.csv"), overwrite=True)
print(t[t["sep"] < 150])
T0 = t["Tmag"][0]
for row in t[1:]:
    if row["sep"] < 100:
        print("TIC", row["ID"], "sep %.1f\"  dT=%.2f  flux ratio %.4f" % (row["sep"], row["Tmag"] - T0, 10 ** (-0.4 * (row["Tmag"] - T0))))
