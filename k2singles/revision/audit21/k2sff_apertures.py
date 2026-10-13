"""
Aperture dependence of the catalogued event in the K2SFF HLSP, for the 21 non-supported events of audit21.
For each K2SFF aperture (CIRC_APER0-9, PRF_APER0-9, and BESTAPER) run k2lib.test_window_k2 on FCOR at the
catalogued epoch, depth and duration. Number of pixels per aperture is taken from CIRC_APER_TBL / PRF_APER_TBL.
A dip that appears in K2SFF but not in PDCSAP/EVEREST and whose depth changes strongly with aperture size points to
a contaminating source or a background problem rather than the target; a dip present at similar depth in every
aperture but absent from FRAW points to the K2SFF arclength correction.
Output: k2singles/revision/audit21/k2sff_apertures.csv
Usage:  python k2singles/revision/audit21/k2sff_apertures.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import audit21 as A  # noqa: E402

warnings.filterwarnings("ignore")


def main():
    from astropy.io import fits
    rows = []
    for _, r in A.targets().iterrows():
        epic, camp, T = int(r.EPIC), int(r.campaign), float(r.epoch_btjd)
        D, delta = r.dur_h / 24, r.depth_ppm * 1e-6
        ps = list((A.DATA / "mastDownload" / "HLSP").glob(f"hlsp_k2sff_k2_lightcurve_{epic}-c{camp:02d}*/*.fits"))
        if not ps:
            continue
        with fits.open(ps[0]) as h:
            npx = {}
            for kind in ("CIRC", "PRF"):
                tbl = h[f"{kind}_APER_TBL"].data
                for i in range(tbl.shape[0]):
                    npx[f"{kind}_APER{i}"] = int(np.sum(tbl[i] > 0))
            for ext in ["BESTAPER"] + [f"{k}_APER{i}" for k in ("CIRC", "PRF") for i in range(10)]:
                a = h[ext].data
                t = np.asarray(a["T"], float) + A.BK
                for col in ("FCOR", "FRAW"):
                    f = np.asarray(a[col], float)
                    t2, f2, _ = A._norm(t, f, np.zeros(len(t), bool), T, D)
                    w = A.window(t2, f2, T, D, delta)
                    rows.append(dict(EPIC=epic, aperture=ext, column=col, npix=npx.get(ext, np.nan),
                                     verdict=w["verdict"], ratio=w.get("depth", np.nan) / delta, snr=w.get("snr")))
        print(epic, "done", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(HERE / "k2sff_apertures.csv", index=False)
    piv = d[d.column == "FCOR"].pivot_table(index="EPIC", columns="aperture", values="ratio")
    pd.set_option("display.width", 250)
    cols = ["BESTAPER"] + [f"CIRC_APER{i}" for i in range(10)]
    print(piv[cols].round(2).to_string())
    print(piv[[f"PRF_APER{i}" for i in range(10)]].round(2).to_string())
    print(d[d.column == "FCOR"].drop_duplicates("aperture")[["aperture", "npix"]].to_string(index=False))


if __name__ == "__main__":
    main()
