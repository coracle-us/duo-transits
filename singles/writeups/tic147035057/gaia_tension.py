"""
Second-reviewer item 2 for TIC 147035057: does a large Gaia RV amplitude fit with a low RUWE?

Inputs: Gaia DR3 values queried 2026-10-11 (source 2947647054662200192):
  rv_amplitude_robust 49.76 km/s (peak-to-peak over 27 transits), rv_chisq_pvalue 0,
  rv_renormalised_gof 99.8, vbroad 25.9 km/s, ruwe 1.118, astrometric_excess_noise 0.118 mas
  (sig 28), ipd_frac_multi_peak 0, ipd_gof_harmonic_amplitude 0.012, duplicated_source True,
  non_single_star 0, parallax 6.023 mas, G 8.83, BP-RP 0.624, FLAME R 2.00 / M 1.44 / L 6.05,
  logg GSP-Phot 3.90, GSP-Spec 4.06, DSC binary-star probability 0.007; TIC: R 2.00, M 1.32.

Computes:
  1  K ~ amplitude/2 -> mass function -> minimum companion mass (e = 0, i ~ 90 deg; the eclipses
     make sin i ~ 1)
  2  photocentre semi-major axis for a dark companion, and for a luminous one
  3  the SB2-twin alternative: two equal stars sharing L = 6.05 Lsun, each star's R and log g
     and absolute G mag, against the single-star values
Output: gaia_tension.json
Run:    python singles/writeups/tic147035057/gaia_tension.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
G = 6.674e-11
MSUN, AU, RSUN = 1.989e30, 1.496e11, 6.957e8
P_D = 37.31497
M1, R1, L1 = 1.44, 2.00, 6.05
PLX = 6.023
AMP = 49.76


def main():
    res = {}
    K = AMP / 2 * 1e3
    P = P_D * 86400
    fm = P * K ** 3 / (2 * np.pi * G) / MSUN
    m2 = brentq(lambda m: m ** 3 / (M1 + m) ** 2 - fm, 1e-3, 5)
    res["K_kms_assumed"] = K / 1e3
    res["mass_function_Msun"] = fm
    res["M2_min_Msun"] = m2
    a = ((M1 + m2) * (P_D / 365.25) ** 2) ** (1 / 3)
    a1 = a * m2 / (M1 + m2)
    res["a_total_AU"] = a
    res["a1_AU"] = a1
    res["photocentre_mas_dark_companion"] = a1 * PLX
    # luminous companion: main-sequence mass-luminosity L ~ M^4 (rough, for the G band)
    beta = m2 ** 4 / (m2 ** 4 + L1)
    res["companion_light_fraction_rough"] = beta
    res["photocentre_mas_luminous_companion"] = a * abs(m2 / (M1 + m2) - beta) * PLX
    # companion radius for an M2 ~ 0.6 Msun dwarf, as a check on the grazing-fit k
    r2 = 0.6
    res["k_if_M2_dwarf_R0.6"] = r2 / R1
    # SB2 twins: two equal stars share L
    Rt = R1 / np.sqrt(2)
    Mt = 1.25          # a plausible mass for each twin (~F6V-F8V, slightly evolved)
    res["twin_each_R_Rsun"] = Rt
    res["twin_each_logg"] = 4.438 + np.log10(Mt) - 2 * np.log10(Rt)
    res["single_logg_FLAME"] = 4.438 + np.log10(M1) - 2 * np.log10(R1)
    res["logg_gspspec"] = 4.06
    res["logg_gspphot"] = 3.90
    res["M_G_total"] = 8.83 - 5 * np.log10(1000 / PLX / 10)
    res["M_G_each_twin"] = res["M_G_total"] + 2.5 * np.log10(2)
    # twin orbit: relative velocity and each star's K for e = 0
    at = ((2 * Mt) * (P_D / 365.25) ** 2) ** (1 / 3)
    vrel = 2 * np.pi * at * AU / P / 1e3
    res["twin_K_each_kms"] = vrel / 2
    res["twin_photocentre_mas"] = 0.0
    res["twin_eclipse_true_depth_on_one_star_ppt"] = 2 * 3.2
    (HERE / "gaia_tension.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
