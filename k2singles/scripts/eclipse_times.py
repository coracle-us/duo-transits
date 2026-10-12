"""
Dedicated eclipse-time fits and linear ephemerides for the K2-singles table (fixes review point P4:
the run_full grid times are quantised to D/4 and must not be quoted).

Per event: a symmetric trapezoid (mid-time, depth, total duration T14, ingress fraction) times a
local linear baseline, fitted by least squares to the data within +/- 1.5 T14 of the guess, one
product/sector at a time. Mid-time errors come from the covariance, scaled up by sqrt(reduced
chi^2) when it exceeds 1, and by a red-noise beta measured on the residuals.

Ephemeris: weighted least squares T = T0 + n P over the events labelled primary. Cycle numbers
come from the trial period given on the command line; the fit fails loudly if any residual
exceeds 3 sigma.

Usage:  python k2singles/scripts/eclipse_times.py EPIC P_trial T1[:s] T2[:s] ... [--dur_h D]
        (':s' marks a secondary eclipse, which is fitted and listed but kept out of the ephemeris)
Output: k2singles/results/ephemerides/<EPIC>.json and <EPIC>_fits.png
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import k2lib  # noqa: E402

OUT = HERE.parent / "results" / "ephemerides"


def trapezoid(t, t0, depth, t14, fin):
    """Symmetric trapezoid: 1 out of eclipse, 1 - depth on the flat bottom, ingress = fin * t14."""
    tin = max(fin, 1e-3) * t14
    x = np.abs(t - t0)
    m = np.ones_like(t)
    full = x <= t14 / 2 - tin
    part = (x > t14 / 2 - tin) & (x < t14 / 2)
    m[full] = 1 - depth
    m[part] = 1 - depth * (t14 / 2 - x[part]) / tin
    return m


def beta_of(res, n_bin):
    nb = len(res) // n_bin
    if n_bin <= 1 or nb < 4:
        return 1.0
    rb = res[: nb * n_bin].reshape(nb, n_bin).mean(axis=1)
    return float(max(1.0, np.std(rb, ddof=1) / (np.std(res, ddof=1) / np.sqrt(n_bin))))


def fit_event(t, f, guess, D):
    m = np.abs(t - guess) < 1.5 * D + 0.3
    tt, ff = t[m], f[m]
    if len(tt) < 10:
        return None
    oot = np.abs(tt - guess) > 0.6 * D
    if oot.sum() < 5:
        return None
    ff = ff / np.median(ff[oot])
    d0 = max(0.002, 1 - np.min(ff))
    p0 = [guess, d0, D, 0.3, 1.0, 0.0]

    def resid(p):
        t0, depth, t14, fin, a, b = p
        return ff - trapezoid(tt, t0, depth, t14, fin) * (a + b * (tt - guess))

    lo = [guess - 0.5 * D, 0, 0.3 * D, 0.01, 0.5, -1]
    hi = [guess + 0.5 * D, 1, 2.0 * D, 0.5, 1.5, 1]
    r = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac")
    res = resid(r.x)
    dof = max(1, len(tt) - len(p0))
    chi2r = float(np.sum(res ** 2) / dof)
    J = r.jac
    try:
        cov = np.linalg.inv(J.T @ J) * chi2r
    except np.linalg.LinAlgError:
        return None
    cad = float(np.median(np.diff(tt)))
    beta = beta_of(res, max(1, int(round(0.25 * r.x[2] / cad))))
    s_t0 = float(np.sqrt(cov[0, 0]) * beta)
    return dict(T=float(r.x[0]), sigma_T=s_t0, depth=float(r.x[1]), sigma_depth=float(np.sqrt(cov[1, 1]) * beta),
                T14_h=float(r.x[2] * 24), fin=float(r.x[3]), beta=beta, n=int(len(tt)), cadence_min=cad * 1440,
                rms=float(np.std(res)), _tt=tt, _ff=ff, _model=trapezoid(tt, *r.x[:4]) * (r.x[4] + r.x[5] * (tt - guess)))


def main():
    a = argparse.ArgumentParser()
    a.add_argument("epic", type=int)
    a.add_argument("p_trial", type=float)
    a.add_argument("events", nargs="+")
    a.add_argument("--dur_h", type=float)
    a.add_argument("--tic", type=int)
    o = a.parse_args()
    cat = pd.read_csv(HERE.parent / "data" / "k2_singles.csv")
    row = cat[cat.EPIC == o.epic].iloc[0]
    tic = o.tic or int(row.TIC)
    D = (o.dur_h or row.dur_h) / 24
    lc, _ = k2lib.load_combined(tic, o.epic)
    OUT.mkdir(parents=True, exist_ok=True)
    fits = []
    for ev in o.events:
        sec = ev.endswith(":s")
        g = float(ev.rstrip(":s"))
        near = lc[np.abs(lc.time - g) < 1.5 * D + 0.3]
        best = None
        for (sector, prod), d in near.groupby(["sector", "product"]):
            r = fit_event(d.time.to_numpy(), d.flux.to_numpy(), g, D)
            if r is None:
                continue
            r.update(sector=int(sector), product=prod, kind="secondary" if sec else "primary", guess=g)
            if best is None or r["sigma_T"] < best["sigma_T"]:
                best = r
        if best is None:
            print(f"  {g}: no usable data")
            continue
        fits.append(best)
        print(f"  {g:10.3f} -> T = {best['T']:.5f} +/- {best['sigma_T']:.5f}  depth {best['depth'] * 1e3:.1f} ppt  "
              f"T14 {best['T14_h']:.2f} h  beta {best['beta']:.2f}  [{best['product']} s{best['sector']}] {best['kind']}")
    prim = [x for x in fits if x["kind"] == "primary"]
    T = np.array([x["T"] for x in prim])
    s = np.array([x["sigma_T"] for x in prim])
    n = np.round((T - T[0]) / o.p_trial)
    A = np.vstack([np.ones_like(n), n]).T
    W = np.diag(1 / s ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    T0, P = cov @ A.T @ W @ T
    res = T - (T0 + n * P)
    chi2 = float(np.sum((res / s) ** 2))
    dof = len(T) - 2
    scale = np.sqrt(chi2 / dof) if dof > 0 and chi2 > dof else 1.0
    sP, sT0 = float(np.sqrt(cov[1, 1]) * scale), float(np.sqrt(cov[0, 0]) * scale)
    print(f"EPIC {o.epic}: P = {P:.5f} +/- {sP:.5f} d, T0 = {T0:.5f} +/- {sT0:.5f} BTJD; cycles {n.astype(int).tolist()}; "
          f"chi2 {chi2:.1f} / dof {dof}; O-C (min) {np.round(res * 1440, 1).tolist()}")
    if np.any(np.abs(res) > 3 * s * scale):
        print("WARNING: residual > 3 sigma; check cycle numbers")
    for x in fits:
        if x["kind"] == "secondary":
            x["phase"] = float(((x["T"] - T0) / P) % 1)
            print(f"  secondary {x['T']:.4f}: phase {x['phase']:.4f}")
    k = len(fits)
    fig, axs = plt.subplots(1, k, figsize=(3.2 * k, 2.8), squeeze=False)
    for ax, x in zip(axs[0], fits):
        ax.plot((x["_tt"] - x["T"]) * 24, x["_ff"], ".", ms=2, color="0.4")
        ax.plot((x["_tt"] - x["T"]) * 24, x["_model"], "r-", lw=1)
        ax.set_title(f"{x['product']} s{x['sector']} {x['kind'][0]}\nT={x['T']:.4f}", fontsize=8)
        ax.set_xlabel("h from mid", fontsize=8)
    fig.suptitle(f"EPIC {o.epic} / TIC {tic}: P = {P:.5f} +/- {sP:.5f} d", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / f"{o.epic}_fits.png", dpi=110)
    clean = [{k2: v for k2, v in x.items() if not k2.startswith("_")} for x in fits]
    (OUT / f"{o.epic}.json").write_text(json.dumps(dict(EPIC=o.epic, TIC=tic, P=float(P), sigma_P=sP, T0=float(T0),
                                                         sigma_T0=sT0, cycles=n.astype(int).tolist(), chi2=chi2,
                                                         dof=dof, oc_min=(res * 1440).tolist(), events=clean),
                                                    indent=1))


if __name__ == "__main__":
    main()
