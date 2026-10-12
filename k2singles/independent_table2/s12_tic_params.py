"""TIC stellar parameters and nearby sources (for dilution / plausibility only)."""
from astroquery.mast import Catalogs
import numpy as np
for tic in [301328232, 14394716, 178941241]:
    r = Catalogs.query_criteria(catalog='TIC', ID=tic)
    c = r[0]
    print(f"TIC {tic}: RA {c['ra']:.5f} Dec {c['dec']:.5f} Tmag {c['Tmag']:.2f} Teff {c['Teff']} logg {c['logg']} "
          f"R {c['rad']} M {c['mass']} rho {c['rho']} lumclass {c['lumclass']} contratio {c['contratio']} dist {c['d']} GAIA {c['GAIA']}")
    nb = Catalogs.query_region(f"{c['ra']} {c['dec']}", radius=0.02, catalog='TIC')
    nb = nb[np.argsort(nb['dstArcSec'])]
    for n in nb[:8]:
        print(f"    nb TIC {n['ID']} sep {n['dstArcSec']:.1f}\" Tmag {n['Tmag']:.2f}")
