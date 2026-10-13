"""Shared data loading: eclipse windows from data/windows.csv, cleaned and centred."""
import os
import numpy as np
import pandas as pd

H = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(H, "data")

# refined rough ephemeris (from visual minima; final values come from the fits)
P_R, T0_R, DSEC_R = 21.7750, 2488.600, 7.185
HALF = {"P": 0.60, "S": 0.70}            # half-width of fit window (d)
OOT = {"P": 0.32, "S": 0.40}             # |dt| beyond this is out of eclipse (for noise/baseline guess)
EXCLUDE = {(-1, "S"), (2, "S"), (-127, "S")}   # sector-edge ramp / scattered-light / coarse pointing


def load(src="own", col="flux", seasons=None):
    W = pd.read_csv(os.path.join(D, "windows.csv"))
    W = W[W.src == src].copy()
    if col != "flux":
        W["flux"] = W[col]
    out = []
    for (n, typ), d in W.groupby(["epoch", "type"]):
        if (n, typ) in EXCLUDE:
            continue
        if seasons is not None and d.season.iloc[0] not in seasons:
            continue
        tc = T0_R + n * P_R + (DSEC_R if typ == "S" else 0.0)
        d = d[np.abs(d.time - tc) < HALF[typ]].copy()
        d = d[np.isfinite(d.flux) & (d.flux > 0)]
        if len(d) < 20:
            continue
        oot = np.abs(d.time - tc) > OOT[typ]
        if oot.sum() >= 10:
            c = np.polyfit(d.time[oot] - tc, d.flux[oot], 1)
            base = np.polyval(c, d.time - tc)
        else:
            base = np.nanmax(d.flux) * np.ones(len(d))
        d["fn"] = d.flux / base
        # outlier rejection: 5 sigma from a running median (window 7) of normalised flux
        rm = d.fn.rolling(7, center=True, min_periods=3).median()
        r = d.fn - rm
        s = 1.4826 * np.nanmedian(np.abs(r - np.nanmedian(r)))
        d = d[np.abs(r) < 5 * s].copy()
        oot = np.abs(d.time - tc) > OOT[typ]
        sig = np.std(d.fn[oot]) if oot.sum() >= 10 else 1.4826 * np.nanmedian(np.abs(np.diff(d.fn))) / np.sqrt(2)
        d["err"] = sig
        d["tc_rough"] = tc
        d["ecl"] = "%s%d" % (typ, n)
        out.append(d)
    return pd.concat(out, ignore_index=True)
