"""Fixed transit shape from the joint fit of the two known events (normalised to unit depth)."""
import json, os
import numpy as np, batman

OUT = os.path.dirname(os.path.abspath(__file__))
FK = json.load(open(os.path.join(OUT, "fit_known.json")))["primary"]["best"]
P_REF = 50.0


def shape(t, tc, cad):
    t = np.ascontiguousarray(t, dtype=np.float64)
    pm = batman.TransitParams()
    a = np.exp(FK["ln_aR"])
    pm.t0, pm.per, pm.rp, pm.a = tc, P_REF, FK["rp"], a
    pm.inc = np.degrees(np.arccos(FK["b"] / a))
    pm.ecc, pm.w, pm.u, pm.limb_dark = 0.0, 90.0, [0.36, 0.23], "quadratic"
    if cad > 300:
        m = batman.TransitModel(pm, t, supersample_factor=11, exp_time=cad / 86400)
    else:
        m = batman.TransitModel(pm, t)
    lc = m.light_curve(pm)
    return (1 - lc) / (1 - lc.min()) if lc.min() < 1 else np.zeros_like(t)
