import numpy as np
from scipy.optimize import least_squares
P0, T00 = 21.77465, 2488.59020     # TESS-only ephemeris (k2singles/results/ephemerides/202126877.json)

def gg(t, tc, D, w, p):
    """generalized Gaussian eclipse profile: depth D (fractional), half-width scale w (d), shape p"""
    return 1 - D*np.exp(-np.abs((t-tc)/w)**p)

def fit_event(t, f, e, tguess, win=0.75, inner=0.30, poly=1, fix_p=None):
    """Detrend with polynomial baseline fitted jointly with a generalized-Gaussian eclipse.
    Returns dict with tc, D, w, p, errors (cov scaled by reduced chi2 and red-noise beta), n, rms."""
    m = np.abs(t - tguess) < win
    t, f, e = t[m], f[m], e[m]
    if len(t) < 20: return None
    x0 = [tguess, 0.2, 0.12, 2.0] + [1.0] + [0.0]*poly
    def model(par):
        tc, D, w, p = par[:4]
        if fix_p is not None: p = fix_p
        base = np.polyval(par[4:][::-1], t - tguess)
        return base*gg(t, tc, D, w, p)
    def res(par): return (f - model(par))/e
    lb = [tguess-0.3, 0, 0.02, 0.8] + [-np.inf]*(poly+1)
    ub = [tguess+0.3, 1, 0.5, 8] + [np.inf]*(poly+1)
    # crude start: deepest point
    i0 = np.argmin(f); x0[0] = float(np.clip(t[i0], tguess-0.25, tguess+0.25)); x0[1] = min(0.9, max(0.01, 1 - f[i0]/np.median(f)))
    r = least_squares(res, x0, bounds=(lb, ub))
    J = r.jac; chi2 = np.sum(r.fun**2); dof = max(1, len(t) - len(x0))
    try: cov = np.linalg.inv(J.T @ J) * chi2/dof
    except Exception: cov = np.full((len(x0), len(x0)), np.nan)
    resid = f - model(r.x)
    # time-averaged beta (red noise) on ~1h bins
    dt = np.median(np.diff(t)); nb = max(2, int(round(1/24/dt)))
    sig1 = np.std(resid); bins = [resid[i:i+nb].mean() for i in range(0, len(resid)-nb, nb)]
    beta = max(1.0, np.std(bins)/(sig1/np.sqrt(nb)*np.sqrt(len(bins)/(len(bins)-1)))) if len(bins) > 3 else 1.0
    err = np.sqrt(np.diag(cov))*beta
    covered_min = np.any(np.abs(t - r.x[0]) < 0.05)
    cov_in = np.sum(np.abs(t - r.x[0]) < 2*r.x[2]); 
    return dict(tc=r.x[0], e_tc=err[0], D=r.x[1], e_D=err[1], w=r.x[2], p=r.x[3], n=len(t), rms=sig1, beta=beta,
                n_in=int(cov_in), covers_min=bool(covered_min), t=t, f=f, model=model(r.x), base=np.polyval(r.x[4:][::-1], t - tguess))
