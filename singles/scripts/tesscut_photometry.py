"""
Own aperture photometry from a TESScut FFI cutout, for sectors that have no SPOC/QLP light
curve yet (e.g. a sector released days ago).

Aperture: pixels within 1.6 px of the target's WCS position whose median flux is above the
local background; background: median of the faintest 40% of pixels per frame. quality == 0.
Output CSV in the same format as the fetched light curves (time, flux, flux_err, quality), so
alias_pipeline.load_lightcurves can use it: data/raw/lightcurves/<TIC>/s<NNN>_OWN_<cad>s.csv.gz
(and 'OWN' is added to the product ranking only by the caller).

Usage:  python singles/scripts/tesscut_photometry.py TIC SECTOR [RA DEC]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import diff_image as di  # noqa: E402


def photometry(tic, sector, ra, dec):
    path = di.get_cutout(tic, ra, dec, sector)
    if path is None:
        raise SystemExit(f"no TESScut cutout for TIC {tic} sector {sector}")
    with fits.open(path) as h:
        tab, wcs = h[1].data, WCS(h[2].header)
    t, flux, q = tab["TIME"], tab["FLUX"].astype(float), tab["QUALITY"]
    tx, ty = wcs.all_world2pix(ra, dec, 0)
    yy, xx = np.mgrid[0:flux.shape[1], 0:flux.shape[2]]
    med = np.nanmedian(flux, axis=0)
    bkg_pix = med <= np.nanpercentile(med, 40)
    bkg = np.nanmedian(flux[:, bkg_pix], axis=1)
    ap = (np.hypot(xx - tx, yy - ty) <= 1.6) & (med > np.nanpercentile(med, 60))
    f = np.nansum(flux[:, ap] - bkg[:, None], axis=1)
    good = np.isfinite(t) & np.isfinite(f) & (q == 0) & (f > 0)
    cad = int(round(np.nanmedian(np.diff(t[good])) * 86400))
    fe = np.sqrt(np.abs(f)) / np.nanmedian(f[good])  # rough Poisson-like error, relative
    df = pd.DataFrame({"time": t[good], "flux": f[good], "flux_err": fe[good] * np.nanmedian(f[good]),
                       "quality": 0})
    out = ROOT / "data" / "raw" / "lightcurves" / str(tic) / f"s{sector:03d}_OWN_{cad}s.csv.gz"
    df.to_csv(out, index=False, compression="gzip")
    print(f"TIC {tic} S{sector}: {len(df)} cadences ({cad} s), aperture {ap.sum()} px -> {out.name}")
    return out


if __name__ == "__main__":
    tic, sector = int(sys.argv[1]), int(sys.argv[2])
    if len(sys.argv) > 4:
        ra, dec = float(sys.argv[3]), float(sys.argv[4])
    else:
        c = pd.read_csv(ROOT / "singles" / "data" / "ctoi_singles.csv").drop_duplicates("TIC").set_index("TIC").loc[tic]
        ra, dec = c.RA, c.Dec
    photometry(tic, sector, ra, dec)
