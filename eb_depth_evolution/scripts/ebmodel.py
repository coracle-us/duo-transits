"""Minimal detached-EB eclipse model: two limb-darkened spheres on a Keplerian orbit (no reflection/ellipsoidal; fine for a
21.8-d detached A+A binary). Units: radii in units of semi-major axis a."""
import numpy as np
from scipy.optimize import brentq
NR = 400
_rr = (np.arange(NR) + 0.5)/NR            # radial nodes on the occulted disc (midpoint rule)
def ld_I(r, u1, u2):
    mu = np.sqrt(1 - r**2); return 1 - u1*(1 - mu) - u2*(1 - mu)**2
def blocked_fraction(z, p, u1, u2):
    """fraction of a quadratic-LD star's flux (radius 1) blocked by an opaque disc of radius p at centre distance z (arrays)."""
    z = np.atleast_1d(z)[:, None]; r = _rr[None, :]
    I = ld_I(_rr, u1, u2)[None, :]
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (r**2 + z**2 - p**2)/(2*r*z)
        frac = np.where(c <= -1, 1.0, np.where(c >= 1, 0.0, np.arccos(np.clip(c, -1, 1))/np.pi))
    frac = np.where(z[:, 0:1] + r <= p, 1.0, frac)            # annulus fully inside occultor
    num = np.sum(I*2*np.pi*r*frac, axis=1)/NR
    tot = np.sum(ld_I(_rr, u1, u2)*2*np.pi*_rr)/NR
    return num/tot
def kepler_E(M, e):
    E = M + e*np.sin(M)
    for _ in range(30):
        E = E - (E - e*np.sin(E) - M)/(1 - e*np.cos(E))
    return E
def M_from_nu(nu, e):
    E = 2*np.arctan2(np.sqrt(1 - e)*np.sin(nu/2), np.sqrt(1 + e)*np.cos(nu/2))
    return E - e*np.sin(E)
def nu_from_M(M, e):
    E = kepler_E(M, e); return 2*np.arctan2(np.sqrt(1 + e)*np.sin(E/2), np.sqrt(1 - e)*np.cos(E/2))
def sec_phase(e, w, i=np.pi/2):
    """phase of secondary conjunction after primary conjunction (conjunction = nu+w = pi/2 [primary], 3pi/2 [secondary])"""
    Mp = M_from_nu(np.pi/2 - w, e); Ms = M_from_nu(3*np.pi/2 - w, e)
    return ((Ms - Mp)/(2*np.pi)) % 1
def ecosw_for_phase(phase, esinw):
    f = lambda ec: sec_phase(np.hypot(ec, esinw), np.arctan2(esinw, ec)) - phase
    return brentq(f, -0.8, 0.8)
def flux(t, tc, kind, P, e, w, inc, r1, r2, J, u=(0.2, 0.27), l3=0.0, texp=None, nsub=1):
    """normalised flux around one eclipse. tc = time of the relevant conjunction (kind 'P' or 'S')."""
    t = np.asarray(t, float)
    if texp and nsub > 1:
        off = (np.arange(nsub) + 0.5)/nsub - 0.5
        tt = (t[:, None] + off[None, :]*texp).ravel()
    else: tt = t
    nu_c = np.pi/2 - w if kind == "P" else 3*np.pi/2 - w
    Mc = M_from_nu(nu_c, e)
    M = Mc + 2*np.pi*(tt - tc)/P
    nu = nu_from_M(M, e)
    rsep = (1 - e**2)/(1 + e*np.cos(nu))
    X = -rsep*np.cos(w + nu); Y = -rsep*np.sin(w + nu)*np.cos(inc); Z = rsep*np.sin(w + nu)*np.sin(inc)
    d = np.hypot(X, Y)
    u1, u2 = u
    Lfac = 1 - u1/3 - u2/6
    L1 = r1**2*Lfac; L2 = r2**2*J*Lfac
    f = np.full(tt.shape, L1 + L2)
    # star 2 in front (Z>0): star 1 eclipsed
    m1 = (Z > 0) & (d < r1 + r2)
    if m1.any(): f[m1] -= L1*blocked_fraction(d[m1]/r1, r2/r1, u1, u2)
    m2 = (Z < 0) & (d < r1 + r2)
    if m2.any(): f[m2] -= L2*blocked_fraction(d[m2]/r2, r1/r2, u1, u2)
    f = f/(L1 + L2)
    f = (f + l3)/(1 + l3)
    if texp and nsub > 1: f = f.reshape(len(t), nsub).mean(axis=1)
    return f
