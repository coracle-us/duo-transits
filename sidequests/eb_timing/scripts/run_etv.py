"""
Per-system eclipse timing -> O-C -> model comparison -> plot.

Usage:  python sidequests/eb_timing/scripts/run_etv.py TIC [TIC ...]
        (ephemeris seed from data/prsa2022_tess_ebs.csv; light curves from cache/lc)
Outputs (sidequests/eb_timing/):
  results/oc/<TIC>_oc.csv        every timed eclipse (BTJD), formal + scaled errors, flags
  results/oc/<TIC>_fit.json      ephemeris, model comparison, LTTE solution
  figs/<TIC>_oc.png              template fold + O-C diagnostics (vetted by eye)

Steps
 1. Seed template from data within 400 d of the catalogue epoch (2-min data if possible),
    phase 0 moved to the template's primary minimum.
 2. Time all primaries with a wide search, fit a linear ephemeris (iterative 5-sigma clip).
 3. Rebuild the template from ALL data (2-min only if >= 3 sectors of it, else everything),
    folded with the improved ephemeris AND with each point shifted by the running-median
    O-C (so a large ETV does not smear the template), then time primaries + secondaries.
 4. Clean: points > 5 robust-sigma from a running median of their own kind are flagged.
 5. Error scaling: z = consecutive-eclipse differences within a sector / combined sigma;
    errors *= max(1, 1.4826 MAD(z)) per kind. This measures white noise only (signals
    with P3 >> sector length do not inflate it).
 6. O-C models fitted jointly to primaries + secondaries (separate offsets):
      lin   c0 + c0s*[sec] + c1*t
      quad  + c2*t^2
      aps   quad + c1s*t*[sec]       (secondary drifts relative to primary: apsidal motion)
      ltte  quad + LTTE(P3, T3, e3, w3, a12 sin i / c), common to both eclipse types
    LTTE: circular-orbit periodogram over P3 in [20 d, 2 x span], best 4 peaks refined with
    eccentric fits from several starts. BIC compared with the scaled errors.
    Also: at the best P3, separate sinusoids for primaries and secondaries -> correlation
    rho = cos(angle between their amplitude vectors): +1 = shift together (LTTE / dynamical),
    -1 = anticorrelated (spots or apsidal motion).
 7. Spot diagnostic: O'Connell difference per sector (mean template-subtracted flux at
    phase 0.20-0.30 minus 0.70-0.80) and its scatter.
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import etvlib as el  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
OCDIR = HERE / "results" / "oc"
FIG = HERE / "figs"
CAT = HERE / "data" / "prsa2022_tess_ebs.csv"
MIN_P3 = 20.0
MAX_P3_SPAN = 1.5   # P3 upper bound in units of the data span (changed from 3 during validation)


def robust_sigma(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan


def fit_linear(t, oc, sig, n_iter=5):
    keep = np.ones(len(t), bool)
    for _ in range(n_iter):
        A = np.column_stack([np.ones(keep.sum()), t[keep]])
        coef, *_ = np.linalg.lstsq(A / sig[keep, None], oc[keep] / sig[keep], rcond=None)
        r = oc - (coef[0] + coef[1] * t)
        s = robust_sigma(r[keep])
        new = np.abs(r) < 5 * s
        if (new == keep).all():
            break
        keep = new
    return coef, keep


def smooth_fn(t, y, n=15):
    """Callable: running median of y(t) (n eclipses), linearly interpolated, constant
    beyond the ends (used as a prior O-C shift)."""
    o = np.argsort(t)
    tt, rm = np.asarray(t)[o], running_median(np.asarray(t)[o], np.asarray(y)[o], n)
    ok = np.isfinite(rm)
    return lambda x: np.interp(x, tt[ok], rm[ok])


def running_median(t, y, n=15):
    o = np.argsort(t)
    s = pd.Series(y[o]).rolling(n, center=True, min_periods=3).median().values
    out = np.empty_like(s)
    out[o] = s
    return out


# ------------------------------------------------------------------ O-C model fitting
class OCData:
    def __init__(self, df, tref):
        self.t = df["T"].values - tref
        self.y = df.oc.values
        self.s = df.sig_s.values
        self.sec = (df.kind == "s").values.astype(float)

    def design(self, kind):
        cols = [np.ones_like(self.t), self.sec, self.t]
        if kind in ("quad", "aps", "ltte"):
            cols.append(self.t ** 2)
        if kind == "aps":
            cols.append(self.t * self.sec)
        if not self.sec.any():
            cols.pop(1)
            if kind == "aps":
                cols.pop(-1)
        return np.column_stack(cols)

    def linfit(self, A, y=None):
        y = self.y if y is None else y
        w = 1 / self.s
        coef, *_ = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)
        r = (y - A @ coef) / self.s
        return coef, float(r @ r)


def bic(chi2, k, n):
    return chi2 + k * np.log(n)


def ltte_search(d, span, base="quad"):
    A0 = d.design(base)
    P_grid = np.exp(np.linspace(np.log(MIN_P3), np.log(MAX_P3_SPAN * span), 4000))
    chi = np.empty(len(P_grid))
    for i, P in enumerate(P_grid):
        ph = 2 * np.pi * d.t / P
        A = np.column_stack([A0, np.sin(ph), np.cos(ph)])
        chi[i] = d.linfit(A)[1]
    # pick up to 4 well-separated minima
    order = np.argsort(chi)
    peaks = []
    for i in order:
        if all(abs(np.log(P_grid[i] / P_grid[j])) > 0.05 for j in peaks):
            peaks.append(i)
        if len(peaks) == 4:
            break
    best = None
    for i in peaks:
        P0 = P_grid[i]
        ph = 2 * np.pi * d.t / P0
        A = np.column_stack([A0, np.sin(ph), np.cos(ph)])
        coef, _ = d.linfit(A)
        a0 = max(np.hypot(coef[-2], coef[-1]), 1e-6)
        for e0 in (0.05, 0.3, 0.6):
            for w0 in (0.5, 2.0, 3.6, 5.2):
                for T30 in (0.0, 0.33, 0.66):
                    p0 = [P0, T30 * P0, e0, w0, a0]

                    def resid(p):
                        lt = el.ltte(d.t, *p)
                        coef, _ = d.linfit(A0, d.y - lt)
                        return (d.y - lt - A0 @ coef) / d.s

                    try:
                        r = least_squares(resid, p0, bounds=([MIN_P3, -np.inf, 0, -np.inf, 0],
                                                             [MAX_P3_SPAN * span, np.inf, 0.95, np.inf, 1.0]),
                                          x_scale=[P0 * 0.05, P0 * 0.1, 0.1, 0.5, a0], max_nfev=400)
                    except Exception:  # noqa: BLE001
                        continue
                    c2 = float(r.fun @ r.fun)
                    if best is None or c2 < best[1]:
                        best = (r.x, c2, r)
    # covariance for errors
    p, c2, r = best
    try:
        J = r.jac
        cov = np.linalg.inv(J.T @ J)
        err = np.sqrt(np.diag(cov) * max(1.0, c2 / max(len(d.t) - J.shape[1] - d.design(base).shape[1], 1)))
    except np.linalg.LinAlgError:
        err = np.full(5, np.nan)
    return P_grid, chi, p, err, c2


def corr_pri_sec(d, P3, base="quad"):
    """Separate sinusoids at P3 for primaries and secondaries; rho = cos(angle)."""
    if d.sec.sum() < 10 or (1 - d.sec).sum() < 10:
        return np.nan, np.nan, np.nan
    ph = 2 * np.pi * d.t / P3
    A0 = d.design(base)
    pri = 1 - d.sec
    A = np.column_stack([A0, np.sin(ph) * pri, np.cos(ph) * pri, np.sin(ph) * d.sec, np.cos(ph) * d.sec])
    coef, _ = d.linfit(A)
    vp, vs = coef[-4:-2], coef[-2:]
    rho = float(vp @ vs / (np.linalg.norm(vp) * np.linalg.norm(vs) + 1e-30))
    return rho, float(np.linalg.norm(vp)), float(np.linalg.norm(vs))


# ------------------------------------------------------------------ spots
def oconnell(lc, tmpl, P, T0):
    rows = []
    for s, g in lc.groupby("sector"):
        ph = ((g.time.values - T0) / P) % 1
        f = g.flux.values
        m1 = (ph > 0.2) & (ph < 0.3)
        m2 = (ph > 0.7) & (ph < 0.8)
        if m1.sum() < 5 or m2.sum() < 5:
            continue
        rows.append(dict(sector=s, t=g.time.median(), oc_diff=np.median(f[m1]) - np.median(f[m2])))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ main per system
def run(tic, cat):
    row = cat[cat.TIC == tic].iloc[0]
    P, T0 = float(row.Per), float(row.BJD0)
    lc = el.load(tic)
    if len(lc) < 1000:
        print(f"TIC {tic}: no data")
        return None
    t, f = lc.time.values, lc.flux.values
    P_cat, T0_cat = P, T0
    P_pdm = el.refine_period(t, f, P, T0)
    sc = (lc["product"] == "SPOC_120").values

    def pass0(P, T0):
        """Seed template near the catalogue epoch, wide-search primaries, linear ephemeris."""
        near = np.abs(t - T0) < 400
        m = near & sc if (near & sc).sum() > 5000 else (near if near.sum() > 3000 else np.ones(len(t), bool))
        tm0 = el.build_template(t[m], f[m], P, T0)
        ph1, *_ = el.find_minima(tm0)
        T0 = T0 + ph1 * P
        tm0 = el.build_template(t[m], f[m], P, T0)
        hw1 = el.eclipse_halfwidth(tm0, 0.0)
        pri0 = el.time_eclipses(lc, tm0, P, T0, 0.0, hw1, search=min(0.05 * P, 0.1), kind="p")
        if len(pri0) < 10:
            return None
        oc0 = pri0["T"].values - (T0 + pri0.epoch.values * P)
        coef, keep = fit_linear(pri0.epoch.values.astype(float), oc0, pri0.sig.values)
        T0, P = T0 + coef[0], P + coef[1]
        oc0 = pri0["T"].values - (T0 + pri0.epoch.values * P)
        spread = robust_sigma(oc0[keep]) if keep.sum() > 5 else 0.01
        return dict(P=P, T0=T0, pri0=pri0, oc0=oc0, keep=keep, spread=spread, n=int(keep.sum()))

    # 1-2. seed ephemeris: try the catalogue period and the PDM-refined period (rule added
    # during the pilot: PDM fixes catalogue periods that are off by ~1e-5 d, e.g. TIC 140661916,
    # but for systems with hour-scale ETVs, e.g. TIC 167692429, the catalogue period can lock
    # better). Keep the catalogue seed unless the PDM seed times >10% more primaries or has
    # less than half the O-C scatter.
    c_cat = pass0(P_cat, T0_cat)
    c_pdm = pass0(P_pdm, T0_cat) if abs(P_pdm - P_cat) > 1e-7 * P_cat else None
    cands = [c for c in (c_cat, c_pdm) if c is not None]
    if not cands:
        print(f"TIC {tic}: too few primaries in pass 0")
        return None
    best0 = c_cat if c_cat is not None else c_pdm
    if c_cat is not None and c_pdm is not None and (c_pdm["n"] > 1.1 * c_cat["n"] or c_pdm["spread"] < 0.5 * c_cat["spread"]):
        best0 = c_pdm
    P, T0, pri0, oc0, keep, spread = (best0[k] for k in ("P", "T0", "pri0", "oc0", "keep", "spread"))
    print(f"TIC {tic}: catalogue P={P_cat:.7f}, PDM P={P_pdm:.7f}; pass-0 candidates "
          f"{[(round(c['P'], 7), c['n'], round(c['spread'] * 86400, 1)) for c in cands]} -> P={P:.7f}", flush=True)

    # 3. iterate: template from points shifted by the running-median O-C of their own
    #    eclipse type (primary half / secondary half of the phase), then re-time with the
    #    running median as the prior so the search can stay narrow.
    n_sc_sectors = lc[sc].sector.nunique()
    mt = sc if n_sc_sectors >= 3 else np.ones(len(t), bool)
    rm_fn = {"p": smooth_fn(pri0["T"].values[keep], oc0[keep])}
    rm_fn["s"] = rm_fn["p"]
    resid_spread = spread
    ph2 = np.nan
    for it in range(3):
        ph_pt = ((t - T0) / P) % 1
        if np.isfinite(ph2):
            near_s = np.abs(((ph_pt - ph2) + 0.5) % 1 - 0.5) < 0.25
        else:
            near_s = np.zeros(len(t), bool)
        shift = np.where(near_s, rm_fn["s"](t), rm_fn["p"](t))
        tmpl = el.build_template(t[mt] - shift[mt], f[mt], P, T0)
        ph1, ph2_new, d1, d2 = el.find_minima(tmpl)
        if abs(ph1) > 0.002:
            T0 = T0 + ph1 * P
            tmpl = el.build_template(t[mt] - shift[mt], f[mt], P, T0)
            ph1, ph2_new, d1, d2 = el.find_minima(tmpl)
        ph2 = ph2_new
        hw1 = el.eclipse_halfwidth(tmpl, 0.0)
        search = float(np.clip(8 * resid_spread + 0.002, 0.003, 0.05 * P))
        pri = el.time_eclipses(lc, tmpl, P, T0, 0.0, hw1, search=search, kind="p", shift_fn=rm_fn["p"])
        parts = [pri]
        hw2 = np.nan
        if np.isfinite(ph2):
            hw2 = el.eclipse_halfwidth(tmpl, ph2)
            sec = el.time_eclipses(lc, tmpl, P, T0, ph2, hw2, search=search, kind="s", shift_fn=rm_fn["s"])
            parts.append(sec)
        oc = pd.concat(parts, ignore_index=True)
        oc["phase_ecl"] = np.where(oc.kind == "p", 0.0, ph2 if np.isfinite(ph2) else 0)
        oc["oc"] = oc["T"] - (T0 + (oc.epoch + oc.phase_ecl) * P)
        # update running medians per kind and the residual spread around them
        res_all = []
        for k, g in oc.groupby("kind"):
            if len(g) < 10:
                continue
            rmv = running_median(g["T"].values, g.oc.values, 15)
            r = g.oc.values - rmv
            ok = np.abs(r) < 5 * max(robust_sigma(r), 1e-6)
            rm_fn[k] = smooth_fn(g["T"].values[ok], g.oc.values[ok])
            res_all.append(r[ok])
        if "s" not in oc.kind.values:
            rm_fn["s"] = rm_fn["p"]
        new_spread = robust_sigma(np.concatenate(res_all)) if res_all else resid_spread
        if it > 0 and abs(new_spread - resid_spread) < 0.1 * resid_spread:
            resid_spread = new_spread
            break
        resid_spread = new_spread

    # 4. clean + 5. error scaling
    oc["flag"] = ""
    oc["sig_s"] = oc.sig
    scale = {}
    for k, g in oc.groupby("kind"):
        rmk = running_median(g["T"].values, g.oc.values, 15)
        r = g.oc.values - rmk
        rs = robust_sigma(r)
        bad = np.abs(r) > 5 * max(rs, 1e-6)
        oc.loc[g.index[bad], "flag"] = "outlier"
        gg = g[~bad].sort_values("T")
        same = (np.diff(gg.sector.values) == 0) & (np.diff(gg.epoch.values) == 1)
        z = (np.diff(gg.oc.values) / np.sqrt(gg.sig.values[1:] ** 2 + gg.sig.values[:-1] ** 2))[same]
        sc_k = max(1.0, robust_sigma(z)) if len(z) > 10 else 1.0
        scale[k] = float(sc_k)
        oc.loc[g.index, "sig_s"] = g.sig * sc_k
    good = oc[oc.flag == ""].copy()
    tref = float(np.round(np.median(good["T"]), 1))
    span = float(good["T"].max() - good["T"].min())
    d = OCData(good, tref)
    n = len(good)
    res = dict(tic=int(tic), P=P, T0=T0, ph_sec=float(ph2) if np.isfinite(ph2) else None,
               hw_pri_phase=hw1, hw_sec_phase=float(hw2) if np.isfinite(hw2) else None,
               depth_pri=float(1 - d1 / np.median(tmpl.y)) if np.isfinite(d1) else None,
               n_pri=int((good.kind == "p").sum()), n_sec=int((good.kind == "s").sum()),
               n_outlier=int((oc.flag != "").sum()), n_sectors=int(lc.sector.nunique()),
               span_d=span, tref=tref, err_scale=scale,
               median_sig_pri_s=float(good[good.kind == "p"].sig_s.median() * 86400),
               median_sig_sec_s=float(good[good.kind == "s"].sig_s.median() * 86400) if (good.kind == "s").any() else None,
               products=lc.groupby("product").sector.nunique().to_dict())
    models = {}
    for k in ("lin", "quad", "aps"):
        A = d.design(k)
        coef, c2 = d.linfit(A)
        models[k] = dict(chi2=c2, k=A.shape[1], bic=bic(c2, A.shape[1], n), coef=coef.tolist())
    fits = {}
    for base in ("lin", "quad"):
        Pg_b, chi_b, p_b, perr_b, c2_b = ltte_search(d, span, base)
        kb = d.design(base).shape[1]
        models["ltte_" + base] = dict(chi2=c2_b, k=kb + 5, bic=bic(c2_b, kb + 5, n))
        fits[base] = (Pg_b, chi_b, p_b, perr_b, c2_b)
    # preferred LTTE variant by BIC
    ltbase = "lin" if models["ltte_lin"]["bic"] <= models["ltte_quad"]["bic"] else "quad"
    # rule added during validation (TIC 219900027): a solution with P3 at the upper bound is
    # a quadratic/LTTE degeneracy, not a measurement; prefer the other variant if it is not
    at_bound = {b: fits[b][2][0] > 0.98 * MAX_P3_SPAN * span for b in fits}
    other = "quad" if ltbase == "lin" else "lin"
    if at_bound[ltbase] and not at_bound[other]:
        ltbase = other
    variants = {b: dict(P3=float(fits[b][2][0]), e3=float(fits[b][2][2]),
                        A_s=float(el.amp_ltte(fits[b][2][4], fits[b][2][2], fits[b][2][3]) * 86400),
                        bic=models["ltte_" + b]["bic"], at_bound=bool(at_bound[b])) for b in fits}
    Pg, chi_g, p, perr, c2 = fits[ltbase]
    kq = d.design(ltbase).shape[1]
    models["ltte"] = dict(models["ltte_" + ltbase])
    # conservative comparison: errors inflated so the best model has reduced chi2 = 1
    best = min(models.values(), key=lambda m: m["bic"])
    c2r_best = max(1.0, best["chi2"] / max(n - best["k"], 1))
    for m_ in models.values():
        m_["bic_cons"] = bic(m_["chi2"] / c2r_best, m_["k"], n)
    P3, T3, e3, w3, asini = p
    A_ltte = el.amp_ltte(asini, e3, w3)
    rho, amp_p, amp_s = corr_pri_sec(d, P3, ltbase)
    res.update(models={k: {kk: vv for kk, vv in v.items() if kk != "coef"} for k, v in models.items()},
               ltte_base=ltbase, ltte_variants=variants, P_catalogue=P_cat, P_pdm=P_pdm, chi2r_best=c2r_best,
               dBIC_ltte_vs_quad=models["quad"]["bic"] - models["ltte"]["bic"],
               dBICc_ltte_vs_quad=models["quad"]["bic_cons"] - models["ltte"]["bic_cons"],
               dBICc_quad_vs_lin=models["lin"]["bic_cons"] - models["quad"]["bic_cons"],
               dBICc_aps_vs_quad=models["quad"]["bic_cons"] - models["aps"]["bic_cons"],
               dBICc_ltte_vs_aps=models["aps"]["bic_cons"] - models["ltte"]["bic_cons"],
               dBIC_quad_vs_lin=models["lin"]["bic"] - models["quad"]["bic"],
               dBIC_aps_vs_quad=models["quad"]["bic"] - models["aps"]["bic"],
               chi2r_quad=models["quad"]["chi2"] / max(n - models["quad"]["k"], 1),
               chi2r_ltte=c2 / max(n - kq - 5, 1), n_used=n,
               ltte=dict(P3=float(P3), P3_err=float(perr[0]), T3=float(T3 + tref), e3=float(e3),
                         e3_err=float(perr[2]), w3_deg=float(np.degrees(w3) % 360),
                         asini_s=float(asini * 86400), asini_err_s=float(perr[4] * 86400),
                         A_ltte_s=float(A_ltte * 86400),
                         f_m3=float(el.mass_function(asini, P3)), P3_over_span=float(P3 / span)),
               rho_pri_sec=rho, amp_pri_s=amp_p * 86400, amp_sec_s=amp_s * 86400)
    # 7. spots
    occ = oconnell(lc, tmpl, P, T0)
    res["oconnell_sd"] = float(occ.oc_diff.std()) if len(occ) > 2 else None
    res["oconnell_mean"] = float(occ.oc_diff.mean()) if len(occ) else None

    # residual from the preferred LTTE model, stored per eclipse (for product-offset checks)
    A_b = d.design(ltbase)
    lt_b = el.ltte(d.t, *p)
    c_b, _ = d.linfit(A_b, d.y - lt_b)
    oc["resid_ltte"] = np.nan
    oc.loc[good.index, "resid_ltte"] = d.y - A_b @ c_b - lt_b
    occ_b = d.design("quad")
    c_q, _ = d.linfit(occ_b)
    oc["resid_quad"] = np.nan
    oc.loc[good.index, "resid_quad"] = d.y - occ_b @ c_q
    OCDIR.mkdir(parents=True, exist_ok=True)
    oc.to_csv(OCDIR / f"{tic}_oc.csv", index=False, float_format="%.7f")
    with open(OCDIR / f"{tic}_fit.json", "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    plot(tic, lc, tmpl, P, T0, ph2, oc, good, d, models, p, tref, Pg, chi_g, occ, res, ltbase)
    print(f"TIC {tic}: P={P:.7f} n_p={res['n_pri']} n_s={res['n_sec']} sig_p={res['median_sig_pri_s']:.1f}s "
          f"dBICc(ltte_{ltbase}-quad)={res['dBICc_ltte_vs_quad']:.1f} P3={P3:.1f} A={A_ltte*86400:.1f}s e={e3:.2f} "
          f"rho={rho:.2f} dBICc(aps)={res['dBICc_aps_vs_quad']:.1f} chi2r_q={res['chi2r_quad']:.2f}", flush=True)
    return res


def plot(tic, lc, tmpl, P, T0, ph2, oc, good, d, models, p, tref, Pg, chi_g, occ, res, ltbase):
    fig = plt.figure(figsize=(13, 10))
    gs = fig.add_gridspec(4, 3, width_ratios=[1, 1, 2.2], hspace=0.35, wspace=0.3)
    # template fold
    ax = fig.add_subplot(gs[0, 0:2])
    ph = ((lc.time.values - T0) / P + 0.25) % 1 - 0.25
    sub = np.random.default_rng(0).choice(len(ph), min(len(ph), 30000), replace=False)
    ax.plot(ph[sub], lc.flux.values[sub], ",", color="0.6")
    xx = np.linspace(-0.25, 0.75, 2000)
    ax.plot(xx, tmpl(xx), "r-", lw=1)
    ax.set_title(f"TIC {tic}  P={P:.6f} d  fold + template", fontsize=9)
    ax.set_xlabel("phase")
    # periodogram
    ax = fig.add_subplot(gs[1, 0:2])
    ax.semilogx(Pg, models[ltbase]["chi2"] - chi_g, "k-", lw=0.6)
    ax.axvline(p[0], color="r", ls=":")
    ax.axvline(res["span_d"], color="b", ls="--", lw=0.8)
    ax.set_xlabel("P3 trial (d)  [blue: span]")
    ax.set_ylabel("Δχ² (sinusoid)")
    # text
    ax = fig.add_subplot(gs[2:, 0:2])
    ax.axis("off")
    L = res["ltte"]
    txt = (f"n_pri={res['n_pri']} n_sec={res['n_sec']} outliers={res['n_outlier']}\n"
           f"sectors={res['n_sectors']} span={res['span_d']:.0f} d\n"
           f"median σ pri={res['median_sig_pri_s']:.1f} s  sec={res['median_sig_sec_s']}\n"
           f"err scale={res['err_scale']}\n"
           f"products={res['products']}\n"
           f"χ²r quad={res['chi2r_quad']:.2f}  ltte={res['chi2r_ltte']:.2f}\n"
           f"ΔBIC, errors scaled so best model χ²r=1 (best χ²r={res['chi2r_best']:.2f}):\n"
           f"  quad vs lin={res['dBICc_quad_vs_lin']:.1f}  aps vs quad={res['dBICc_aps_vs_quad']:.1f}\n"
           f"  ltte({ltbase}) vs quad={res['dBICc_ltte_vs_quad']:.1f}  ltte vs aps={res['dBICc_ltte_vs_aps']:.1f}\n"
           f"LTTE: P3={L['P3']:.1f}±{L['P3_err']:.1f} d  A={L['A_ltte_s']:.1f} s\n"
           f"   e={L['e3']:.2f}±{L['e3_err']:.2f} ω={L['w3_deg']:.0f}°  P3/span={L['P3_over_span']:.2f}\n"
           f"   f(m3)={L['f_m3']:.3g} Msun\n"
           f"pri/sec corr at P3: ρ={res['rho_pri_sec']:.2f} "
           f"(A_p={res['amp_pri_s']:.1f}s, A_s={res['amp_sec_s']:.1f}s)\n"
           f"O'Connell sd={res['oconnell_sd']}")
    ax.text(0, 1, txt, va="top", family="monospace", fontsize=8)

    # O-C
    tt = good["T"].values
    ax1 = fig.add_subplot(gs[0, 2])
    Aq = d.design(ltbase)
    lt = el.ltte(d.t, *p)
    cq, _ = d.linfit(Aq, d.y - lt)
    base = Aq @ cq  # quad part of the ltte model
    off_sec = cq[1] if d.sec.any() else 0
    y_plot = (d.y - (base - off_sec * d.sec - cq[0])) * 1440  # remove linear+quad trend, keep LTTE
    ylin = d.y.copy()
    cl, _ = d.linfit(d.design("lin"))
    ylin = (d.y - d.design("lin") @ cl) * 1440
    for k, c in (("p", "C0"), ("s", "C3")):
        mk = (good.kind == k).values
        ax1.errorbar(tt[mk], ylin[mk], good.sig_s.values[mk] * 1440, fmt=".", ms=2, color=c, alpha=0.25, elinewidth=0.3)
    # sector means
    for k, c in (("p", "navy"), ("s", "darkred")):
        g = good[good.kind == k]
        if len(g) == 0:
            continue
        mk = (good.kind == k).values
        dfm = pd.DataFrame({"t": tt[mk], "y": ylin[mk], "w": 1 / good.sig_s.values[mk] ** 2, "s": g.sector.values})
        bm = dfm.groupby("s").apply(lambda x: pd.Series({"t": x.t.mean(), "y": np.average(x.y, weights=x.w),
                                                         "e": 1 / np.sqrt(x.w.sum())}), include_groups=False)
        ax1.errorbar(bm.t, bm.y, bm.e, fmt="o", ms=4, color=c, mfc="white")
    tg = np.linspace(tt.min(), tt.max(), 3000)
    cols_g = [np.ones_like(tg)] + ([np.zeros_like(tg)] if d.sec.any() else []) + [tg - tref]
    if ltbase == "quad":
        cols_g.append((tg - tref) ** 2)
    Aq_g = np.column_stack(cols_g)
    model_g = Aq_g @ cq + el.ltte(tg - tref, *p)
    lin_g = (np.column_stack([np.ones_like(tg), np.zeros_like(tg), tg - tref]) if d.sec.any()
             else np.column_stack([np.ones_like(tg), tg - tref])) @ cl
    ax1.plot(tg, (model_g - lin_g) * 1440, "k-", lw=0.8)
    lo, hi = np.nanpercentile(ylin, [1, 99])
    pad = 0.3 * (hi - lo)
    ax1.set_ylim(lo - pad, hi + pad)
    ax1.set_ylabel("O−C (min), lin removed")
    ax1.set_title("blue=primary red=secondary; open=sector means; line=quad+LTTE", fontsize=8)
    # residuals after ltte model
    ax2 = fig.add_subplot(gs[1, 2], sharex=ax1)
    full = Aq @ cq + lt
    r = (d.y - full) * 1440
    for k, c in (("p", "C0"), ("s", "C3")):
        mk = (good.kind == k).values
        ax2.plot(tt[mk], r[mk], ".", ms=2, color=c, alpha=0.3)
    lo, hi = np.nanpercentile(r, [1, 99])
    ax2.set_ylim(lo - 0.3 * (hi - lo), hi + 0.3 * (hi - lo))
    ax2.set_ylabel("resid (min)")
    # primary - secondary (sector means)
    ax3 = fig.add_subplot(gs[2, 2], sharex=ax1)
    if d.sec.any():
        mp = good.kind == "p"
        a = pd.DataFrame({"s": good.sector, "y": ylin, "k": good.kind}).groupby(["s", "k"]).y.mean().unstack()
        tm = good.groupby("sector")["T"].mean()
        if "s" in a and "p" in a:
            dd = (a["p"] - a["s"]).dropna()
            ax3.plot(tm.loc[dd.index], dd.values, "ko", ms=3)
    ax3.set_ylabel("pri−sec (min)")
    ax3.axhline(0, color="0.7", lw=0.5)
    ax4 = fig.add_subplot(gs[3, 2], sharex=ax1)
    if len(occ):
        ax4.plot(occ.t, occ.oc_diff * 1e3, "go", ms=3)
    ax4.set_ylabel("O'Connell (ppt)")
    ax4.set_xlabel("BTJD")
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{tic}_oc.png", dpi=90, bbox_inches="tight")
    plt.close(fig)


def main():
    """Options (before the TICs): --ephem=CSV (TIC, Per, BJD0 seed ephemeris instead of the
    Prsa catalogue), --max-sector=N (ignore sectors > N), --tag=NAME (write outputs to
    results/oc_NAME/ and figs/NAME/ instead of results/oc/ and figs/)."""
    global OCDIR, FIG
    args = sys.argv[1:]
    cat = pd.read_csv(CAT).drop_duplicates("TIC")
    while args and args[0].startswith("--"):
        k, v = args.pop(0)[2:].split("=", 1)
        if k == "ephem":
            cat = pd.read_csv(v, comment="#").drop_duplicates("TIC")
        elif k == "max-sector":
            el.MAX_SECTOR = int(v)
        elif k == "tag":
            OCDIR = HERE / "results" / f"oc_{v}"
            FIG = HERE / "figs" / v
    for tic in [int(x) for x in args]:
        try:
            run(tic, cat)
        except Exception as e:  # noqa: BLE001
            import traceback
            traceback.print_exc()
            print(f"TIC {tic}: FAILED {e}")


if __name__ == "__main__":
    main()
