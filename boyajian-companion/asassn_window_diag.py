import sys; sys.path.insert(0, ".")
import numpy as np, pandas as pd
from astropy.time import Time
import ground_check as gc
d = pd.read_csv("ground_photometry_merged.csv")
cov = pd.read_csv("epoch_coverage.csv")
print("groups:", d.groupby("group").size().to_dict(), "last good point", Time(d.jd.max(), format="jd").iso)
for name, n, pmin, pmax in gc.WINDOWS:
    a, b = gc.T0 + n*pmin, gc.T0 + n*pmax
    x = d[(d.jd > a - gc.CORE) & (d.jd < b + gc.CORE)]
    r = cov[cov.window == name]
    used = (r.groups.notna() & (r.groups != "")).mean()
    nb = d[(d.jd > a - 15) & (d.jd < b + 15)].groupby("group").size().to_dict()
    print(f"{name:22s} {Time(a,format='jd').iso[:10]}..{Time(b,format='jd').iso[:10]} pts {len(x):2d} nights {np.floor(x.jd).nunique():2d} | "
          f"tc with a usable group {used:4.0%} | |dchi2| max {r.dchi2.abs().max():6.1f} pp {r.dchi2_per_point.abs().max():7.1f} | +-15 d: {nb}")
# scatter of nightly means over 30-d stretches: what one night can do against a 0.012 mag dip
for g, x in d.groupby("group"):
    nm = x.groupby(np.floor(x.jd)).mag.mean().sort_index()
    print(g, "nights", len(nm), "overall nightly MAD", round(gc.mad(nm.values), 4), "first-diff scatter/sqrt2", round(gc.mad(np.diff(nm.values))/np.sqrt(2), 4),
          "5-sigma nights needed (from first-diff)", round((5*gc.mad(np.diff(nm.values))/np.sqrt(2)/gc.DEPTH_MAG)**2, 1))
