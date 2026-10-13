"""Step 3: load all light curves, report coverage, plot raw curves."""
import os, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits

H = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(H, "data")
files = sorted(glob.glob(os.path.join(D, "mastDownload", "**", "*.fits"), recursive=True))
fig, axes = plt.subplots(len(files), 1, figsize=(14, 2.2 * len(files)))
for ax, f in zip(axes, files):
    h = fits.open(f)
    hd = h[1].header
    cols = h[1].columns.names
    d = h[1].data
    t = d["TIME"] if "TIME" in cols else d["T"]
    name = os.path.basename(f)
    print(name, "cols:", cols[:12])
    y = None
    for c in ["PDCSAP_FLUX", "SAP_FLUX", "KSPSAP_FLUX", "FLUX", "FCOR", "DETFLUX"]:
        if c in cols:
            y = np.array(d[c], float)
            break
    good = np.isfinite(t) & np.isfinite(y)
    print("   using", c, "N=", good.sum(), "t:", np.nanmin(t[good]), np.nanmax(t[good]))
    ax.plot(t[good], y[good] / np.nanmedian(y[good]), ".", ms=1)
    ax.set_title(name[:60] + " " + c, fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(H, "fig_raw_lcs.png"), dpi=80)

# TESScut headers
for f in sorted(glob.glob(os.path.join(D, "tesscut", "*.fits"))):
    h = fits.open(f)
    t = h[1].data["TIME"]
    q = h[1].data["QUALITY"]
    print(os.path.basename(f), h[1].data["FLUX"].shape, np.nanmin(t), np.nanmax(t),
          "q!=0:", (q != 0).sum(), "ccd", h[0].header.get("CAMERA"), h[0].header.get("CCD"),
          "sector", h[0].header.get("SECTOR"))
