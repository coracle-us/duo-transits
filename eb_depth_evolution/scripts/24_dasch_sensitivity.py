"""Fraction of random epochs at which EPIC 202126877 would be >0.3 (>0.4) mag below maximum, as a function of inclination, from the
base physical model (blue-ish plates: LD 0.5/0.25, J=0.85). Combined with the DASCH faint-tail counts gives the sensitivity."""
import numpy as np, json, sys
from scipy.stats import binom
sys.path.insert(0, "scripts")
from ebmodel import flux, ecosw_for_phase
p = json.load(open("data/physical_fit_base.json"))["params"]
ec = ecosw_for_phase(0.33221, p["esinw"]); e = np.hypot(ec, p["esinw"]); w = np.arctan2(p["esinw"], ec)
Pa = 21.7748; tt = np.linspace(-0.7, 0.7, 2801); dt = tt[1] - tt[0]
n_t, k_t = 1206, 8; base_rate = 0.0064
print(" inc   f(>0.3mag)  f(>0.4mag)  expected N(>0.3) of 1206  P(N<=8)")
for inc in [85.0, 86.0, 87.26, 87.91, 88.77, 90.0]:
    fr = []
    for thr in (0.3, 0.4):
        tot = 0
        for kind in "PS":
            f = flux(tt, 0, kind, Pa, e, w, np.radians(inc), p["r1"], p["r2"], 0.85, u=(0.5, 0.25))
            tot += np.sum(-2.5*np.log10(f) > thr)*dt
        fr.append(tot/Pa)
    lam = n_t*(base_rate + fr[0])
    print(f"{inc:6.2f}  {fr[0]:.4f}      {fr[1]:.4f}      {lam:6.1f}                {binom.cdf(k_t, n_t, base_rate + fr[0]):.2e}")
