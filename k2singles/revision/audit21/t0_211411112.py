"""
EPIC 211411112 (Sanders 617, M67) is the one Table 1 target among the 21 audited events. Its Table 1 period uses
a C5 eclipse time fitted to K2 PDCSAP, in which the eclipse core is missing (ingress only). This script re-times the
C5 eclipse in the reductions where the full eclipse is present (EVEREST, K2SFF, SAP) with a simple supersampled
trapezoid + linear baseline, and times the C16 eclipse in PDCSAP with the same code as a control. It then compares
with the PDCSAP-ingress time used for Table 1 (k2singles/results/ephemerides/211411112.json) and with the
C16 + S46 ephemeris (no C5 data at all).
Output: k2singles/revision/audit21/t0_211411112.csv and stdout.
Usage:  python k2singles/revision/audit21/t0_211411112.py
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import audit21 as A  # noqa: E402

warnings.filterwarnings("ignore")
CAD = 29.4244 / 1440


def trap(t, t0, dep, t14, t12):
    x = np.abs(t - t0)
    t12 = np.clip(t12, 1e-3, t14 / 2)
    return np.where(x >= t14 / 2, 0.0, np.where(x <= t14 / 2 - t12, dep, dep * (t14 / 2 - x) / t12))


def model(p, t):
    t0, dep, t14, t12, a, b = p
    off = (np.arange(10) + 0.5) / 10 * CAD - CAD / 2
    m = np.mean([trap(t + o, t0, dep, t14, t12) for o in off], axis=0)
    return (a + b * (t - t0)) * (1 - m)


def fit(t, f, guess):
    w = np.abs(t - guess) < 0.5
    t, f = t[w], f[w]
    p0 = [guess, 0.025, 0.17, 0.05, np.median(f), 0.0]
    lo, hi = [guess - 0.05, 0, 0.08, 0.005, 0.5, -1], [guess + 0.05, 0.2, 0.3, 0.15, 1.5, 1]
    r = least_squares(lambda p: model(p, t) - f, p0, bounds=(lo, hi))
    res = r.fun
    J = r.jac
    s2 = np.sum(res ** 2) / max(1, len(t) - len(p0))
    cov = np.linalg.pinv(J.T @ J) * s2
    return dict(t0=r.x[0], sig_t0=np.sqrt(cov[0, 0]), depth=r.x[1], t14_h=r.x[2] * 24, n=len(t),
                rms=np.sqrt(s2))


def main():
    eph = json.loads((A.ROOT / "k2singles" / "results" / "ephemerides" / "211411112.json").read_text())
    ev = {round(e["T"]): e for e in eph["events"]}
    t16, t46 = ev[1118]["T"], ev[2558]["T"]
    P_2 = (t46 - t16) / 23
    pred_c5 = t16 - 15 * P_2
    P_lit = 62.59484
    rows = []
    prods, _ = A.load_products(211411112, 5, 178.9526, 4.4112 / 24)
    for name in ("EVEREST", "K2SFF", "SAP_q", "PDC_q"):
        t, f, flag = prods[name]
        use = ~flag
        try:
            r = fit(t[use], f[use], 178.955)
        except Exception as e:  # noqa: BLE001
            r = dict(error=repr(e))
        rows.append(dict(event="C5", product=name, **r))
    # control: C16 in PDCSAP with the same code
    p16, _ = A.load_products(211411112, 16, 1117.876, 4.4112 / 24)
    t, f, _ = p16["PDC_q"]
    rows.append(dict(event="C16 (control)", product="PDC_q", **fit(t, f, 1117.876)))
    d = pd.DataFrame(rows)
    d["table1_fit_t0_C5"] = ev[179]["T"]
    d["pred_C5_from_C16_S46"] = pred_c5
    d["pred_C5_from_C16_and_Plit"] = t16 - 15 * P_lit
    d["dt_vs_table1_min"] = (d.t0 - ev[179]["T"]) * 1440
    d["dt_vs_C16S46_pred_min"] = (d.t0 - pred_c5) * 1440
    d.to_csv(HERE / "t0_211411112.csv", index=False)
    pd.set_option("display.width", 250)
    print(d.round(5).to_string(index=False))
    print(f"P from C16+S46 only (n=23): {P_2:.6f} d; Table 1: {eph['P']:.5f} +- {eph['sigma_P']:.5f}; "
          f"Sandquist+2021: {P_lit}")
    # period from C5(EVEREST) + C16 + S46
    c5 = d[(d.event == "C5") & (d["product"] == "EVEREST")].iloc[0]
    T = np.array([c5.t0, t16, t46])
    s = np.array([c5.sig_t0, ev[1118]["sigma_T"], ev[2558]["sigma_T"]])
    n = np.array([0, 15, 38])
    W = np.diag(1 / s ** 2)
    X = np.column_stack([np.ones(3), n])
    cov = np.linalg.inv(X.T @ W @ X)
    beta = cov @ X.T @ W @ T
    chi2 = float(np.sum(((T - X @ beta) / s) ** 2))
    print(f"P with C5 re-timed in EVEREST: {beta[1]:.6f} +- {np.sqrt(cov[1, 1]):.6f} d, chi2 = {chi2:.2f} (dof 1)")


if __name__ == "__main__":
    main()
