"""Stage 1b: crowding proxy for QLP / ELEANOR stars whose TIC contratio is empty.

TIC v8 gives contratio only for CTL stars; 309 of 338 pilot QLP/ELEANOR sectors had none.
Proxy (flagged crowd_source = 'tic_proxy'): flux of TIC neighbours within 63" (3 TESS
pixels), each weighted by a Gaussian PSF of sigma 16" at its separation, over the target
flux, both from Tmag. Comparable in spirit to contratio, not identical (no aperture model).
"""
import numpy as np, pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import Catalogs
from common import PILOT, RESULTS

SIG = 16.0  # arcsec


def proxy(tic, ra, dec):
    t = Catalogs.query_region(SkyCoord(ra, dec, unit="deg"), radius=63 * u.arcsec,
                              catalog="TIC").to_pandas()
    me = t[t.ID.astype(str) == str(tic)]
    if me.empty or not np.isfinite(me.Tmag.iloc[0]):
        return np.nan, np.nan
    tm = me.Tmag.iloc[0]
    nb = t[(t.ID.astype(str) != str(tic)) & np.isfinite(t.Tmag)]
    nb = nb[nb.disposition.fillna("") != "ARTIFACT"] if "disposition" in nb else nb
    w = np.exp(-0.5 * (nb.dstArcSec / SIG) ** 2)
    return float(np.sum(w * 10 ** (-0.4 * (nb.Tmag - tm)))), float(tm)


def main():
    sec = pd.read_csv(PILOT / "sectors.csv")
    cat = pd.read_parquet(RESULTS / "catalogue.parquet").set_index("tic")
    rows = []
    for tic in sec.tic.unique():
        c = cat.loc[tic]
        p, tm = proxy(tic, c.ra, c.dec)
        rows.append(dict(tic=tic, contratio_proxy=p, tmag=tm))
        print(tic, p, tm, flush=True)
    pr = pd.DataFrame(rows)
    pr.to_csv(PILOT / "crowding_proxy.csv", index=False)
    print(pr.describe())


if __name__ == "__main__":
    main()
