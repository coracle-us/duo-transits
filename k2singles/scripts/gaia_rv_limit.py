"""
Gaia DR3 radial-velocity limit on the companion of TOI-5560.01 (EPIC 211924561, TIC 238606535).

DR3 gives only summary RV statistics. Two conservative limits on the semi-amplitude K:
  (a) K <= rv_amplitude_robust / 2 (half the robust peak-to-peak of the epoch RVs);
  (b) per-epoch scatter sigma_ep = (radial_velocity_error - 0.11) * sqrt(2N/pi), following the DR3
      definition of radial_velocity_error for bright stars (Katz et al. 2023). A circular orbit adds
      variance K^2/2, so the 99.7% upper limit on the excess variance gives K_max.
M2 sin i then follows from the mass function for the host mass M1 (Gaia FLAME), circular orbit.
Usage:  python k2singles/scripts/gaia_rv_limit.py
"""
import numpy as np
from scipy import optimize, stats

G = 6.674e-11
MSUN, MJUP, DAY = 1.989e30, 1.898e27, 86400.0

# Gaia DR3 source 660640844292494336, queried 2026-10-12 (k2singles/data/lead_status.csv; FLAME mass
# from gaiadr3.astrophysical_parameters)
RV_ERR, N, AMP_ROBUST = 0.601236, 22, 7.698715       # km/s, transits, km/s
M1 = 1.424                                           # Msun, FLAME


def m2sini(K_kms, P_d, m1=M1):
    """Companion mass (Msun) for semi-amplitude K, circular orbit (solve the full mass function)."""
    f = (K_kms * 1e3) ** 3 * P_d * DAY / (2 * np.pi * G) / MSUN
    return optimize.brentq(lambda m2: m2 ** 3 / (m1 + m2) ** 2 - f, 1e-6, 10)


def main():
    sig_ep = (RV_ERR - 0.11) * np.sqrt(2 * N / np.pi)
    # 99.7% upper limit on the sample variance s^2, given (N-1) s^2 / sigma^2 ~ chi2(N-1)
    s2_max = sig_ep ** 2 * stats.chi2.ppf(0.997, N - 1) / (N - 1)
    K_b = np.sqrt(2 * max(s2_max - sig_ep ** 2, 0))
    K_a = AMP_ROBUST / 2
    print(f"per-epoch scatter {sig_ep:.2f} km/s; K limits: (a) {K_a:.2f}  (b) {K_b:.2f} km/s")
    for P in (140.49392, 70.246961):
        print(f"P = {P:9.5f} d: M2 sin i < {m2sini(K_a, P):.3f} Msun (a), < {m2sini(K_b, P):.3f} Msun (b); "
              f"1 MJ gives K = {1e-3 * (2 * np.pi * G / (P * DAY)) ** (1 / 3) * MJUP / (M1 * MSUN) ** (2 / 3):.3f} km/s")


if __name__ == "__main__":
    main()
