#!/usr/bin/env python3
"""Check of ind_common.lnp_excluded_fraction against a literal loop over grid periods.

Purpose: the fast ln P fraction maps runs of excluding grid times onto period-grid index ranges;
this compares it with the direct period-by-period lookup on a thinned period grid.
Inputs: work/raw/<tic>.npz (a star already loaded by ind_arm_a.py). Outputs: prints both numbers.
Usage: python test_lnp.py <tic> <t_cat> <depth_ppm> <dur_h>
"""
import sys
import numpy as np
import ind_common as ic
from ind_common import core, pl

tic, t_cat, depth, dur = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]) * 1e-6, float(sys.argv[4]) / 24
secs, _ = ic.load_star_raw(tic)
T, F = ic.concat([pl.detrend(s["t"], s["f"], dur) for s in secs])
sig, edges = core.duration_scale_noise(T, F, dur), core.gap_edges(T)
a = ic.lnp_excluded_fraction(t_cat, T, F, dur, depth, 0.1 * depth, sig, edges)
b = ic.lnp_excluded_fraction_brute(t_cat, T, F, dur, depth, 0.1 * depth, sig, edges)
print("fast", a, "brute (thinned grid)", b)
