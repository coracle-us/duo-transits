"""Inventory every downloaded FITS: header info, time span, columns, cadence counts."""
import glob, os
import numpy as np
from astropy.io import fits
HERE = os.path.dirname(os.path.abspath(__file__))
for f in sorted(glob.glob(os.path.join(HERE, "data", "**", "*.fits"), recursive=True)):
    with fits.open(f) as h:
        h0 = h[0].header
        d = h[1].data
        cols = d.columns.names
        tcol = "TIME" if "TIME" in cols else ("T" if "T" in cols else cols[0])
        t = np.asarray(d[tcol], float)
        good = np.isfinite(t)
        print(os.path.basename(f))
        print("   obj=%s sector=%s camp=%s cam=%s ccd=%s  n=%d  t=[%.3f, %.3f]  ext=%d" % (
            h0.get("OBJECT"), h0.get("SECTOR"), h0.get("CAMPAIGN"), h0.get("CAMERA"), h0.get("CCD"),
            good.sum(), np.nanmin(t), np.nanmax(t), len(h)))
        print("   cols:", cols[:25])
        if "TIMEREF" in h[1].header or "BJDREFI" in h[1].header:
            print("   BJDREFI=", h[1].header.get("BJDREFI"), "BJDREFF=", h[1].header.get("BJDREFF"), "TIMESYS=", h[1].header.get("TIMESYS"))
