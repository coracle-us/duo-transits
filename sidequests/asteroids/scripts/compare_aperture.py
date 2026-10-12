"""
Validation-time method check: hard pixel-centre aperture at the Horizons position (version 1)
vs soft-edged aperture re-centred on the stacked source (version 2). For each object, writes the
light curve with each variant, runs period.analyse, and records the residual rms (30-min bins)
and the recovered period. Output: data/aperture_comparison.csv
Usage: python compare_aperture.py NUM:SECTOR:P_known [NUM:SECTOR:P_known ...]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mtphot  # noqa: E402
import period  # noqa: E402

HERE = Path(__file__).resolve().parents[1]


def main():
    rows = []
    for arg in sys.argv[1:]:
        num, sec, pk = arg.split(":")
        obj = f"{num};"
        tab = mtphot.load_cutout(mtphot.fetch_cutout(obj, int(sec)))
        eph = mtphot.fetch_ephem(obj, tab["t"][np.isfinite(tab["t"])])
        ph1 = mtphot.photometry(tab, soft=False)
        c = mtphot.centroid_check(tab, ph1)
        ph2 = mtphot.photometry(tab, c["cen_dx"], c["cen_dy"], soft=True)
        for tag, ph in [("v1_hard_horizons", ph1), ("v2_soft_recentred", ph2)]:
            t = tab["t"]
            ok = np.isfinite(t)
            e = {k: np.interp(np.where(ok, t, np.nanmedian(t)), t[ok], v) for k, v in
                 mtphot.fetch_ephem(obj, t[ok]).items()}
            df = pd.DataFrame(dict(time=t, t_lt=t - e["delta"] / mtphot.C_AU_PER_DAY, flux=ph["F"],
                                   flux_err=ph["err"], flag=ph["flag"], star_ap=ph["star_ap"],
                                   r=e["r"], delta=e["delta"], alpha=e["alpha"], V=e["V"]))
            sid = f"cmp_{num}_{tag}"
            df.to_csv(mtphot.LCDIR / f"{sid}_s{int(sec):03d}.csv.gz", index=False, compression="gzip")
            o = period.analyse(sid, int(sec), float(pk), plot=False)
            rows.append(dict(number=int(num), sector=int(sec), variant=tag, cen_dx=round(c["cen_dx"], 3),
                             cen_dy=round(c["cen_dy"], 3), resid_rms_mag=o.get("resid_rms_mag"),
                             P2_h=o.get("P2_h"), P_known_h=float(pk), amp2=o.get("amp2")))
            (mtphot.LCDIR / f"{sid}_s{int(sec):03d}.csv.gz").unlink()
            print(rows[-1], flush=True)
        del eph
    pd.DataFrame(rows).to_csv(HERE / "data" / "aperture_comparison.csv", index=False)


if __name__ == "__main__":
    main()
