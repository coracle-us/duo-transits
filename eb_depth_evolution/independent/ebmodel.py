"""Eclipsing-binary light-curve model: two limb-darkened spheres on an eccentric orbit.

Geometry (relative orbit of star 2 about star 1, observer along +z):
  r      = a (1 - e^2) / (1 + e cos f)
  z      = r sin(w + f) sin i        (> 0: star 2 in front of star 1)
  d_sky  = r sqrt(1 - sin^2 i sin^2(w + f))
Primary eclipse (star 1 occulted by star 2) occurs at w + f = pi/2; T0 = that conjunction time.
Occulted fraction for a quadratically limb-darkened star uses the Mandel & Agol (2002)
routine from batman (_quadratic_ld), which handles radius ratios > 1.
Flux = [ (1 - a1) + Lr (1 - a2) ] / (1 + Lr) * (1 - l3) + l3,  Lr = J * (r2/r1)^2,
J = ratio of disk-averaged surface brightness (star 2 / star 1) in the bandpass.
"""
import numpy as np
from batman import _quadratic_ld


def kepler_E(M, e):
    M = np.mod(M, 2 * np.pi)
    E = M + e * np.sin(M) if e < 0.8 else np.pi * np.ones_like(M)
    for _ in range(50):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E = E - dE
        if np.max(np.abs(dE)) < 1e-12:
            break
    return E


def true_anom(t, P, Tperi, e):
    M = 2 * np.pi * (t - Tperi) / P
    E = kepler_E(M, e)
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))


def tperi_from_t0(T0, P, e, w):
    """Time of periastron given time of primary conjunction (w + f = pi/2)."""
    fc = np.pi / 2 - w
    Ec = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(fc / 2), np.sqrt(1 + e) * np.cos(fc / 2))
    Mc = Ec - e * np.sin(Ec)
    return T0 - Mc * P / (2 * np.pi)


def secondary_offset(P, e, w):
    """Time from primary conjunction to secondary conjunction (w + f = 3 pi / 2)."""
    def tconj(fc):
        Ec = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(fc / 2), np.sqrt(1 + e) * np.cos(fc / 2))
        return (Ec - e * np.sin(Ec)) * P / (2 * np.pi)
    dt = tconj(3 * np.pi / 2 - w) - tconj(np.pi / 2 - w)
    return np.mod(dt, P)


def occ(z, p, ua, ub):
    """Unocculted flux fraction; nudges batman's degenerate points (z = |1-p|, z = p)."""
    z = np.array(z, float)
    bad = (np.abs(z - abs(1 - p)) < 1e-6) | (np.abs(z - p) < 1e-6)
    z[bad] += 2e-6
    return _quadratic_ld._quadratic_ld(np.ascontiguousarray(z), p, ua, ub, 1)


def lc(t, P, T0, e, w, inc_deg, r1, r2, J, u1a, u1b, u2a, u2b, l3=0.0, texp=0.0, nsub=1):
    t = np.asarray(t, float)
    if nsub > 1 and texp > 0:
        off = (np.arange(nsub) + 0.5) / nsub - 0.5
        tt = (t[:, None] + off[None, :] * texp).ravel()
    else:
        tt = t
    Tp = tperi_from_t0(T0, P, e, w)
    f = true_anom(tt, P, Tp, e)
    inc = np.radians(inc_deg)
    r = (1 - e ** 2) / (1 + e * np.cos(f))
    s = np.sin(w + f)
    z = r * s * np.sin(inc)
    d = r * np.sqrt(np.clip(1 - np.sin(inc) ** 2 * s ** 2, 0, None))
    k = r2 / r1
    F1 = np.ones_like(tt)
    F2 = np.ones_like(tt)
    m1 = (z > 0) & (d < r1 + r2)
    m2 = (z < 0) & (d < r1 + r2)
    if m1.any():
        F1[m1] = occ(d[m1] / r1, k, u1a, u1b)
    if m2.any():
        F2[m2] = occ(d[m2] / r2, 1.0 / k, u2a, u2b)
    Lr = J * k ** 2
    F = (F1 + Lr * F2) / (1 + Lr)
    F = F * (1 - l3) + l3
    if nsub > 1 and texp > 0:
        F = F.reshape(-1, nsub).mean(axis=1)
    return F


def _numeric_occult(z, p, ua, ub, n=4000):
    """Brute-force check of the occulted flux fraction (for testing)."""
    r = (np.arange(n) + 0.5) / n
    mu = np.sqrt(1 - r ** 2)
    I = 1 - ua * (1 - mu) - ub * (1 - mu) ** 2
    out = []
    for zz in np.atleast_1d(z):
        c = (r ** 2 + zz ** 2 - p ** 2) / (2 * r * max(zz, 1e-12))
        th = np.where(c <= -1, np.pi, np.where(c >= 1, 0.0, np.arccos(np.clip(c, -1, 1))))
        if zz < 1e-12:
            th = np.where(r < p, np.pi, 0.0)
        blocked = np.sum(I * 2 * th * r) / np.sum(I * 2 * np.pi * r)
        out.append(1 - blocked)
    return np.array(out)


if __name__ == "__main__":
    z = np.array([0.1, 0.5, 0.9, 1.2, 1.5, 2.0, 2.4])
    for p in [0.4, 0.9, 1.3, 2.0]:
        a = _quadratic_ld._quadratic_ld(z, p, 0.3, 0.2, 1)
        b = _numeric_occult(z, p, 0.3, 0.2)
        print("p=%.1f max|batman-numeric|=%.2e" % (p, np.max(np.abs(a - b))))
    # sanity: secondary offset for e cos w ~ -0.27
    e, w = 0.3, np.radians(150)
    print("sec offset", secondary_offset(21.78, e, w), "approx",
          21.78 * (0.5 + 2 / np.pi * e * np.cos(w)))
