"""Step 2: time every primary and secondary eclipse with the S58 template.

Per eclipse: weighted least squares (shift, depth scale, offset, slope) over |phase|<0.06,
moving-block bootstrap of residuals (0.05 d blocks, 500 resamples) for the error,
coverage cuts (RULES 3 + change log 2), then per-sector/type error scaling.
Writes independent/eclipse_times.csv and independent/out/dropped_eclipses.csv.
"""
import os, sys, pickle
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(__file__))
import hiplib as H

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
NBOOT = 500
BLOCK_D = 0.05
rng = np.random.default_rng(7666)


def fit_one(dt, y, w, model, p0):
    def res(p):
        return (y - model(dt, *p)) * w
    r = least_squares(res, p0, x_scale=[1e-3, 0.1, 1e-3, 1e-3], method="lm")
    return r.x, r


def time_sector(sector, tmpl, log, data=None, texp=None, nboot=None):
    texp = H.MAIN[sector]["texp"] if texp is None else texp
    NB = NBOOT if nboot is None else nboot
    if data is None:
        d = np.load(os.path.join(OUT, "prepared_s%02d.npz" % sector))
        t, y, e = d["t"], d["flux"], d["err"]
    else:
        t, y, e = data
    cad = texp / 86400.0
    L = max(1, int(round(BLOCK_D / cad)))
    Emin = int(np.floor((t.min() - H.T0_LEE) / H.P_LEE)) - 1
    Emax = int(np.ceil((t.max() - H.T0_LEE) / H.P_LEE)) + 1
    rows, dropped = [], []
    for kind, off in (("P", 0.0), ("S", 0.5)):
        model = H.eclipse_model_factory(tmpl, kind, texp)
        for E in range(Emin, Emax + 1):
            tp = H.T0_LEE + H.P_LEE * (E + off)
            dph = (t - tp) / H.P_LEE
            if np.all(np.abs(dph) > H.FIT_HALF):
                continue
            n_exp = 2 * H.ECL_HALF * H.P_LEE / cad
            n_in = np.sum(np.abs(dph) < H.ECL_HALF)
            side_exp = (H.ECL_HALF - 0.015) * H.P_LEE / cad
            n_ing = np.sum((dph > -H.ECL_HALF) & (dph < -0.015))
            n_egr = np.sum((dph > 0.015) & (dph < H.ECL_HALF))
            frac, fi, fe = n_in / n_exp, n_ing / side_exp, n_egr / side_exp
            cyc = E + off
            if frac < 0.8 or fi < 0.5 or fe < 0.5:
                if n_in > 0 or np.any(np.abs(dph) < H.FIT_HALF):
                    dropped.append(dict(sector=sector, type=kind, cycle=cyc, t_pred=tp,
                                        frac_in=frac, frac_ingress=fi, frac_egress=fe))
                continue
            m = np.abs(dph) < H.FIT_HALF
            dt, yy, ww = t[m] - tp, y[m], 1.0 / e[m]
            p0 = [0.0, 1.0, 1.0, 0.0]
            p, r = fit_one(dt, yy, ww, model, p0)
            best = model(dt, *p)
            resid = yy - best
            n = len(dt)
            shifts = np.empty(NB)
            nb = int(np.ceil(n / L))
            for b in range(NB):
                starts = rng.integers(0, max(1, n - L + 1), nb)
                idx = (starts[:, None] + np.arange(L)[None, :]).ravel()[:n]
                yb = best + resid[idx]
                pb, _ = fit_one(dt, yb, ww, model, p)
                shifts[b] = pb[0]
            err = np.std(shifts, ddof=1)
            rows.append(dict(sector=sector, type=kind, cycle=cyc, t_btjd=tp + p[0], err_raw=err,
                             shift_s=p[0] * 86400, scale=p[1], n=n, frac_in=frac,
                             rms_ppt=np.std(resid) * 1e3))
            log.write("S%d %s E=%.1f shift=%+.1f s err=%.1f s scale=%.3f n=%d rms=%.2f ppt\n" %
                      (sector, kind, cyc, p[0] * 86400, err * 86400, p[1], n, np.std(resid) * 1e3))
    return rows, dropped


def scale_errors(rows, log):
    for s in sorted(set(r["sector"] for r in rows)):
        for k in ("P", "S"):
            sel = [r for r in rows if r["sector"] == s and r["type"] == k]
            if len(sel) < 3:
                f = 1.0
                chi2nu = np.nan
            else:
                x = np.array([r["cycle"] for r in sel]); y = np.array([r["t_btjd"] for r in sel])
                sg = np.array([r["err_raw"] for r in sel])
                A = np.vstack([np.ones_like(x), x - x.mean()]).T / sg[:, None]
                c, *_ = np.linalg.lstsq(A, y / sg, rcond=None)
                chi2 = np.sum((y / sg - A @ c) ** 2)
                chi2nu = chi2 / (len(sel) - 2)
                f = np.sqrt(chi2nu) if chi2nu > 1 else 1.0
            for r in sel:
                r["err_d"] = r["err_raw"] * f
                r["scale_factor"] = f
            log.write("scaling S%d %s: N=%d chi2nu=%.3f factor=%.3f\n" % (s, k, len(sel), chi2nu, f))
            print("scaling S%d %s: N=%d chi2nu=%.3f factor=%.3f" % (s, k, len(sel), chi2nu, f))
    return rows


if __name__ == "__main__":
    os.chdir(os.path.join(HERE, ".."))
    with open(os.path.join(OUT, "template_s58.pkl"), "rb") as fh:
        tmpl = H.Template(pickle.load(fh))
    log = open(os.path.join(OUT, "times_log.txt"), "w")
    allrows, alldrop = [], []
    for s in (18, 58, 85):
        rows, dropped = time_sector(s, tmpl, log)
        allrows += rows; alldrop += dropped
        print("S%d: %d primaries, %d secondaries kept; %d dropped" %
              (s, sum(r["type"] == "P" for r in rows), sum(r["type"] == "S" for r in rows), len(dropped)))
    allrows = scale_errors(allrows, log)
    allrows.sort(key=lambda r: r["t_btjd"])
    with open(os.path.join(HERE, "eclipse_times.csv"), "w") as fh:
        fh.write("sector,type,cycle,t_btjd,err_d,err_d_unscaled,scale_factor,shift_vs_lee_s,depth_scale,n_points,rms_ppt\n")
        for r in allrows:
            fh.write("%d,%s,%.1f,%.7f,%.7f,%.7f,%.4f,%.2f,%.4f,%d,%.3f\n" %
                     (r["sector"], "primary" if r["type"] == "P" else "secondary", r["cycle"],
                      r["t_btjd"], r["err_d"], r["err_raw"], r["scale_factor"], r["shift_s"],
                      r["scale"], r["n"], r["rms_ppt"]))
    with open(os.path.join(OUT, "dropped_eclipses.csv"), "w") as fh:
        fh.write("sector,type,cycle,t_pred_btjd,frac_in,frac_ingress,frac_egress\n")
        for r in alldrop:
            fh.write("%d,%s,%.1f,%.5f,%.3f,%.3f,%.3f\n" % (r["sector"], r["type"], r["cycle"], r["t_pred"],
                                                          r["frac_in"], r["frac_ingress"], r["frac_egress"]))
