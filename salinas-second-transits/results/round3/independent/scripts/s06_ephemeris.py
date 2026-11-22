"""Linear ephemeris fit to the timing-source mid-times in out/times_<tic>.csv.

Rows used: timing_source == 1 and (by default) at least one of ingress/egress covered and
in-transit coverage >= 0.4.  Weighted least squares Tc = T0 + E*P.  T0 is quoted at the
epoch that minimises the P-T0 covariance (closest integer epoch) and also at epoch 0 (the
paper's event).  If reduced chi2 > 1 the errors are also quoted scaled by sqrt(chi2_red).
Usage: s06_ephemeris.py TIC [--exclude sector:epoch ...] [--full-only]
"""
import sys, os, csv, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import OUT


def load(tic):
    rows = list(csv.DictReader(open(os.path.join(OUT, f"times_{tic}.csv"))))
    return [r for r in rows if r["tc"]]


def fit(E, T, S):
    w = 1 / S**2
    A = np.vstack([np.ones_like(E), E]).T
    C = np.linalg.inv(A.T @ (A * w[:, None]))
    x = C @ (A.T @ (w * T))
    res = T - A @ x
    chi2 = np.sum(res**2 * w)
    return x, C, res, chi2


def run(tic, exclude=(), full_only=False, quiet=False):
    rows = load(tic)
    sel = []
    for r in rows:
        if r["timing_source"] != "1":
            continue
        cov, ing, egr = float(r["cov_in"]), float(r["ingress_cov"]), float(r["egress_cov"])
        if f"{r['sector']}:{r['epoch']}" in exclude:
            continue
        if cov < 0.4 or max(ing, egr) < 0.5:
            continue
        if full_only and min(ing, egr) < 0.5:
            continue
        sel.append(r)
    E = np.array([int(r["epoch"]) for r in sel], float)
    T = np.array([float(r["tc"]) for r in sel])
    S = np.array([float(r["tc_err"]) for r in sel])
    x, C, res, chi2 = fit(E, T, S)
    dof = len(T) - 2
    rchi2 = chi2 / dof if dof > 0 else np.nan
    # pivot epoch
    Ep = int(round(-C[0, 1] / C[1, 1]))
    T0p = x[0] + Ep * x[1]
    sT0p = np.sqrt(C[0, 0] + 2 * Ep * C[0, 1] + Ep**2 * C[1, 1])
    out = dict(tic=tic, n=len(T), P=x[1], sP=np.sqrt(C[1, 1]), T0_E0=x[0], sT0_E0=np.sqrt(C[0, 0]),
               Epiv=Ep, T0_piv=T0p, sT0_piv=sT0p, chi2=chi2, dof=dof, rchi2=rchi2,
               scale=max(1.0, np.sqrt(rchi2)) if dof > 0 else 1.0,
               used=[(r["sector"], r["source"], int(r["epoch"]), float(r["tc"]), float(r["tc_err"]),
                      float(r["cov_in"]), float(r["ingress_cov"]), float(r["egress_cov"])) for r in sel])
    if not quiet:
        print(f"TIC {tic}: N={len(T)}  P = {x[1]:.6f} ± {out['sP']:.6f} d   "
              f"T0(E=0) = {x[0]:.5f} ± {out['sT0_E0']:.5f}   T0(E={Ep}) = {T0p:.5f} ± {sT0p:.5f}")
        print(f"   chi2 = {chi2:.2f} / dof {dof}  -> reduced chi2 = {rchi2:.2f}")
        print("   sector source     E    Tc            err[min]  O-C[min]  O-C/err  cov ing egr")
        for r, e_, t_, s_, o in zip(sel, E, T, S, res):
            flag = "  <-- >3 sigma" if abs(o / s_) > 3 else ""
            print(f"   {r['sector']:>4} {r['source']:<9} {int(e_):4d}  {t_:.5f}  {s_*1440:6.2f}  {o*1440:+7.2f}  {o/s_:+6.2f}"
                  f"   {float(r['cov_in']):.2f} {float(r['ingress_cov']):.2f} {float(r['egress_cov']):.2f}{flag}")
        # also: non-timing rows vs ephemeris (comparison only)
        print("   comparison-only rows (not in fit):")
        for r in rows:
            if r in sel:
                continue
            e_ = int(r["epoch"]); t_ = float(r["tc"]); s_ = float(r["tc_err"])
            o = t_ - (x[0] + e_ * x[1])
            print(f"   {r['sector']:>4} {r['source']:<9} {e_:4d}  {t_:.5f}  {s_*1440:6.2f}  {o*1440:+7.2f}  {o/s_:+6.2f}"
                  f"   {float(r['cov_in']):.2f} {float(r['ingress_cov']):.2f} {float(r['egress_cov']):.2f}"
                  f"  timing_src={r['timing_source']}")
    return out


if __name__ == "__main__":
    tic = int(sys.argv[1])
    ex = [a for a in sys.argv[2:] if ":" in a]
    full = "--full-only" in sys.argv
    o = run(tic, ex, full)
    tag = "_full" if full else ""
    json.dump(o, open(os.path.join(OUT, f"ephem_{tic}{tag}.json"), "w"), indent=1)
