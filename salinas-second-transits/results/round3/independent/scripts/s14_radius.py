"""Stellar radius/density for a target (usage: s14_radius.py TIC) from Gaia DR3 parallax + G, since DR3 gives no
FLAME/GSP-Phot parameters for it. L from M_G with the DR3 BC_G(Teff) polynomial (Creevey+2022,
4000-8000 K); extinction A_G = 2.74 E(B-V) with E(B-V) from the TIC; Teff from the TIC.
Mass for the density from the TIC (and a +-10% range). Compared with the transit-implied
circular-orbit density from out/shape_<tic>.json."""
import sys, os, json, numpy as np
from astroquery.mast import Catalogs
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TIC = int(sys.argv[1])
_lines = []
_print = print
def print(*a):
    _lines.append(" ".join(str(x) for x in a)); _print(*a)
g = json.load(open(os.path.join(ROOT, "out", f"gaia_{TIC}.json")))["target"]
sh = json.load(open(os.path.join(ROOT, "out", f"shape_{TIC}.json")))
r = Catalogs.query_criteria(catalog="TIC", ID=TIC)[0]
teff, eteff = float(r["Teff"]), float(r["e_Teff"])
ebv = float(r["ebv"]) if str(r["ebv"]) != "--" else 0.0
eebv = float(r["e_ebv"]) if str(r["e_ebv"]) != "--" else 0.5 * ebv
mass = float(r["mass"])
print(f"TIC: Teff={teff:.0f}±{eteff:.0f}  E(B-V)={ebv:.3f}±{eebv:.3f}  M={mass:.2f}  R_TIC={float(r['rad']):.3f}  logg={float(r['logg']):.2f}")
a = [6.000e-02, 6.731e-05, -6.647e-08, 2.859e-11, -7.197e-15]


def radius(teff, ag, plx):
    bc = sum(c * (teff - 5772) ** i for i, c in enumerate(a))
    mg = g["G"] - ag + 5 * np.log10(plx / 1000) + 5
    mbol = mg + bc
    L = 10 ** (-0.4 * (mbol - 4.74))
    return np.sqrt(L) * (5772 / teff) ** 2, L


rng = np.random.default_rng(1)
N = 20000
T = rng.normal(teff, eteff, N)
AG = np.clip(rng.normal(2.74 * ebv, 2.74 * eebv, N), 0, None)
PL = rng.normal(g["plx"] + 0.02, 0.015, N)  # +0.02 mas: approx. DR3 parallax zero-point
R, L = radius(T, AG, PL)
M = mass * rng.normal(1, 0.1, N)
rho = M / R**3
P = sh["P0"]
rho_tr = 3 * np.pi * sh["a"] ** 3 / (6.674e-8 * (P * 86400) ** 2) / 1.41
q = lambda x: np.percentile(x, [16, 50, 84])
print(f"Gaia-based: L = {q(L)[1]:.2f} (+{q(L)[2]-q(L)[1]:.2f}/-{q(L)[1]-q(L)[0]:.2f}) Lsun, "
      f"R = {q(R)[1]:.3f} (+{q(R)[2]-q(R)[1]:.3f}/-{q(R)[1]-q(R)[0]:.3f}) Rsun, "
      f"rho = {q(rho)[1]:.2f} (+{q(rho)[2]-q(rho)[1]:.2f}/-{q(rho)[1]-q(rho)[0]:.2f}) rho_sun")
print(f"Transit-implied (circular, b={sh['b']:.2f}, a/R*={sh['a']:.2f}): rho = {rho_tr:.2f} rho_sun")
print(f"Fraction of Gaia-based posterior with rho <= transit value: {np.mean(rho <= rho_tr):.4f}")
Rneed = (np.median(M) / rho_tr) ** (1 / 3)
print(f"Radius needed for the circular-orbit density: {Rneed:.2f} Rsun")
# eccentricity needed: duration ratio = (rho_star/rho_tr)^(1/3) = (1+e sin w)/sqrt(1-e^2); minimum e at w=90
ratio = (np.median(rho) / rho_tr) ** (1 / 3)
es = np.linspace(0, 0.95, 9501)
if ratio > 1:
    emin = es[np.argmax((1 + es) / np.sqrt(1 - es**2) >= ratio)]
    print(f"Duration / circular-expectation = {ratio:.2f} -> minimum eccentricity (transit at apoapsis) e >= {emin:.2f}")
else:
    print(f"Duration / circular-expectation = {ratio:.2f} (<= 1: a circular orbit with b > fitted b, or e > 0 "
          f"near periapsis, fits; no eccentricity needed)")
open(os.path.join(ROOT, "out", f"radius_{TIC}.txt"), "w").write("\n".join(_lines) + "\n")
