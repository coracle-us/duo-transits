"""Injection model shared by s3_search.py and s5_cnn.py (PLAN.md section 7): a batman quadratic limb-darkened
transit with Rp/R* = sqrt(depth_draw), impact parameter b, and a/R* set so that T14 equals the drawn duration at
P = 100 d. Multiplies SAP_FLUX before any processing."""
import numpy as np

U_LD = [0.40, 0.25]   # quadratic limb darkening (solar-type star, TESS band)
P_INJ = 100.0         # d; only sets a/R* together with the duration (one transit per sector)
EXP_S = 200.0         # S104 FFI cadence


def model(t, inj):
    import batman
    p = batman.TransitParams()
    k, b, dur = np.sqrt(inj["depth_draw"]), inj["b"], inj["dur_h"] / 24
    a = P_INJ * np.sqrt((1 + k) ** 2 - b ** 2) / (np.pi * dur)
    p.t0, p.per, p.rp, p.a, p.inc = inj["T"], P_INJ, k, a, np.degrees(np.arccos(b / a))
    p.ecc, p.w, p.u, p.limb_dark = 0.0, 90.0, U_LD, "quadratic"
    good = np.isfinite(t)
    m = np.ones_like(t, dtype=float)
    m[good] = batman.TransitModel(p, t[good], supersample_factor=3, exp_time=EXP_S / 86400).light_curve(p)
    return m
