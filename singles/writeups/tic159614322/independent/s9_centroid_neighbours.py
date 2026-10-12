"""In-transit minus out-of-transit centroid (own photometry) for the S21 dip vs the known S14/S41 events,
(Gaia neighbour query not performed.)"""
import json, os
import numpy as np
from phot import photometry
from lcio import OUT
E = json.load(open(os.path.join(OUT, "epochs.json")))
T14 = E["T14_d"]
tc21 = json.load(open(os.path.join(OUT, "highn_density_ephem.json")))["tc_S21"]["value"]
res = {}
for s, tc, hi in [(14, E["TC1"], None), (41, E["TC2"], None), (21, tc21, 1883.90)]:
    ph = photometry(s)
    t, q = ph["time"], ph["quality"]
    ok = (q == 0) & np.isfinite(ph["cx"])
    if hi:
        ok &= t < hi
    intr = ok & (np.abs(t - tc) < T14 / 2 - 0.03)
    oot = ok & (np.abs(t - tc) > T14 / 2 + 0.03) & (np.abs(t - tc) < 0.8)
    o = {}
    for c in ["cx", "cy"]:
        # remove a linear drift fitted to out-of-transit points
        p = np.polyfit(t[oot] - tc, ph[c][oot], 1)
        r_in = ph[c][intr] - np.polyval(p, t[intr] - tc)
        r_out = ph[c][oot] - np.polyval(p, t[oot] - tc)
        o[c] = [float(np.mean(r_in) * 1e3), float(np.std(r_out) / np.sqrt(intr.sum()) * 1e3)]
    res[f"S{s}"] = dict(tc=tc, n_in=int(intr.sum()), n_out=int(oot.sum()), shift_mpix=o)
    print(f"S{s}: in-transit centroid shift col {o['cx'][0]:+.2f}±{o['cx'][1]:.2f} mpix, row {o['cy'][0]:+.2f}±{o['cy'][1]:.2f} mpix (n_in={intr.sum()})")
# (Gaia DR3 neighbour query removed: the archive call hung from this machine; not done.)
json.dump(res, open(os.path.join(OUT, "centroid_neighbours.json"), "w"), indent=1)
