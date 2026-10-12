"""
O-C noise model, LTTE fit with honest uncertainties, and automatic classification.

Works on any table of eclipse times with columns T (BTJD), kind ('p'/'s'), sector, sig (d),
oc (d, relative to any linear ephemeris), so external times of minima (e.g. Borkovits+2026
on request) can be fitted the same way (see compare_minima.py).

Noise model (per system), all fitted jointly with the deterministic O-C model:
  O-C_i = base(t_i) [+ LTTE(t_i)]                     deterministic (shared by p and s)
          + k_i * delta_s(i)                          spot term: per-sector ANTIsymmetric
                                                      offset (k=+1 primary, -1 secondary),
                                                      free, only for sectors with >= 3 of
                                                      each eclipse type. Starspots shift
                                                      primary and secondary minima in
                                                      opposite directions (Tran+2013); a
                                                      light-time effect shifts both equally,
                                                      so it is orthogonal to this term.
          + gamma_s(i)                                per-sector COMMON jitter, random effect
                                                      N(0, tau^2) (ridge rows in the LS)
          + eps_i,  eps_i ~ N(0, (s_w * sig_i)^2)     white, s_w >= 1 rescales the formal sig
  tau and s_w are estimated (method of moments) from the residuals of the best-fitting model
  and then held fixed for every model, so models are compared on the same noise.
Uncertainties: wild cluster bootstrap over sectors (Rademacher sign per sector applied to the
  residual block of that sector, which keeps any within-sector correlation and the sector
  jitter), B refits started from the best solution; sigma = half the 16-84% interval.
  Optionally multiplied by a global calibration factor CAL (from the validation pulls).
Classification (thresholds fixed on the southern pilot before the NCVZ false-alarm test):
  third-body flag if ALL of
    dBIC = BIC(best of lin/quad/aps) - BIC(LTTE) > DBIC_MIN
    A / sigma_A(raw bootstrap, before CAL) > SNR_MIN
    P3 not at the search bound
    if both eclipse types: rho(pri,sec) >= RHO_MIN and A_pri/A_sec in [1/AMP_RATIO, AMP_RATIO]
  class 'LTTE' if P3/span <= 0.7 else 'LTTE-long'; otherwise 'no-flag'  (reported as cls_v1).
  v2 (added after vetting the NCVZ false-alarm set; reported as cls): additionally
    P3 >= 5 P1, A >= 3 tau, and a common non-linear O-C (quad or LTTE vs linear) with
    dBIC > DBIC_MIN in every leave-one-season-out fit.
"""
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import etvlib as el

MIN_P3 = 20.0
MAX_P3_SPAN = 1.5
MIN_PER_KIND = 3
B_BOOT = 150
CAL = 1.0            # global sigma calibration factor (set from validation pulls)
DBIC_MIN = 30.0
SNR_MIN = 4.0         # tuned on the southern pilot (was 5): TIC 278861826 has A/sigma_A = 4.6 at P3~span
RHO_MIN = 0.8         # tuned on the southern pilot (was 0.5): spotted 25078924 (0.71), 308451526 (0.69)
AMP_RATIO = 1.5       # tuned on the southern pilot (was 2.0)
LONG_SPAN = 0.7
MIN_P3_P1 = 5.0       # v2 (after NCVZ false-alarm vetting)
A_TAU_MIN = 3.0       # v2


def robust_sigma(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan


class OC:
    """Design matrices for one system."""

    def __init__(self, df, spot=True):
        df = df.sort_values("T").reset_index(drop=True)
        self.df = df
        self.tref = float(np.round(np.median(df["T"]), 1))
        self.t = df["T"].values - self.tref
        self.y = df.oc.values.copy()
        self.sig0 = df.sig.values.copy()
        self.k = np.where(df.kind.values == "s", -1.0, 1.0)
        self.has_sec = bool((df.kind == "s").sum() >= 10 and (df.kind == "p").sum() >= 10)
        self.span = float(df["T"].max() - df["T"].min())
        sec = df.sector.values
        self.sectors = np.unique(sec)
        self.S = (sec[:, None] == self.sectors[None, :]).astype(float)    # sector dummies
        cnt_p = pd.Series(df.kind == "p").groupby(sec).sum()
        cnt_s = pd.Series(df.kind == "s").groupby(sec).sum()
        both = [s for s in self.sectors if cnt_p.get(s, 0) >= MIN_PER_KIND and cnt_s.get(s, 0) >= MIN_PER_KIND]
        self.spot_cols = (self.S[:, np.isin(self.sectors, both)] * self.k[:, None]) if (spot and self.has_sec) \
            else np.zeros((len(df), 0))
        self.tau = 0.0
        self.sw = 1.0

    def base(self, kind):
        cols = [np.ones_like(self.t)]
        if self.has_sec:
            cols.append((self.k < 0).astype(float))
        cols.append(self.t)
        if kind in ("quad", "aps"):
            cols.append(self.t ** 2)
        if kind == "aps" and self.has_sec:
            cols.append(self.t * (self.k < 0))
        return np.column_stack(cols)

    def solve(self, A, y):
        """GLS via ridge-augmented LS: fixed columns A + spot columns, random sector
        offsets with prior N(0, tau^2). Returns (coef_fixed, chi2, resid_unwhitened, ranef)."""
        sig = self.sig0 * self.sw
        X = np.column_stack([A, self.spot_cols])
        nfix = X.shape[1]
        if self.tau > 0:
            X = np.column_stack([X, self.S])
        Xw = X / sig[:, None]
        yw = y / sig
        if self.tau > 0:
            prior = np.zeros((len(self.sectors), X.shape[1]))
            prior[:, nfix:] = np.eye(len(self.sectors)) / self.tau
            Xw = np.vstack([Xw, prior])
            yw = np.concatenate([yw, np.zeros(len(self.sectors))])
        coef, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        rw = yw - Xw @ coef
        fixed = X[:, :nfix] @ coef[:nfix]
        ranef = coef[nfix:] if self.tau > 0 else np.zeros(len(self.sectors))
        return coef[:A.shape[1]], float(rw @ rw), y - fixed, ranef

    def nfix(self, A):
        return A.shape[1] + self.spot_cols.shape[1]

    def projector(self, A):
        """Q (orthonormal basis of the whitened, prior-augmented linear design) for fixed
        tau and s_w; the profiled whitened residual of any data vector y is
        yw - Q Q^T yw with yw = [y/sig, 0...]. Cached per design."""
        key = (id(A), A.shape, self.tau, self.sw)
        if getattr(self, "_pk", None) == key:
            return self._pq
        sig = self.sig0 * self.sw
        X = np.column_stack([A, self.spot_cols])
        nfix = X.shape[1]
        if self.tau > 0:
            X = np.column_stack([X, self.S])
        Xw = X / sig[:, None]
        if self.tau > 0:
            prior = np.zeros((len(self.sectors), X.shape[1]))
            prior[:, nfix:] = np.eye(len(self.sectors)) / self.tau
            Xw = np.vstack([Xw, prior])
        Q, _ = np.linalg.qr(Xw)
        self._pk, self._pq = key, Q
        return Q

    def whiten(self, y):
        yw = y / (self.sig0 * self.sw)
        if self.tau > 0:
            yw = np.concatenate([yw, np.zeros(len(self.sectors))])
        return yw

    def pres(self, Q, y):
        yw = self.whiten(y)
        return yw - Q @ (Q.T @ yw)


def bic(chi2, k, n):
    return chi2 + k * np.log(n)


def fit_ltte(d, base, y=None, p0=None, search=True):
    """LTTE + base under the current noise model. If search: sinusoid periodogram over P3
    (projected, cheap) + multi-start eccentric refinement of the best 4 peaks; else a single
    refinement from p0. Returns (p, chi2) or None."""
    y = d.y if y is None else y
    A0 = d.base(base)
    Q = d.projector(A0)
    hi = MAX_P3_SPAN * d.span

    def resid(p):
        return d.pres(Q, y - el.ltte(d.t, *p))

    starts = []
    if search:
        r0 = d.pres(Q, y)
        Pg = np.exp(np.linspace(np.log(MIN_P3), np.log(hi), 2000))
        gain = np.empty(len(Pg))
        for i, P in enumerate(Pg):
            ph = 2 * np.pi * d.t / P
            B = np.column_stack([d.pres(Q, np.sin(ph)), d.pres(Q, np.cos(ph))])
            c, *_ = np.linalg.lstsq(B, r0, rcond=None)
            gain[i] = r0 @ (B @ c)
        order = np.argsort(-gain)
        peaks = []
        for i in order:
            if all(abs(np.log(Pg[i] / Pg[j])) > 0.05 for j in peaks):
                peaks.append(i)
            if len(peaks) == 4:
                break
        for i in peaks:
            P0 = Pg[i]
            ph = 2 * np.pi * d.t / P0
            B = np.column_stack([d.pres(Q, np.sin(ph)), d.pres(Q, np.cos(ph))])
            c, *_ = np.linalg.lstsq(B, r0, rcond=None)
            a0 = max(np.hypot(c[0], c[1]), 1e-6)
            for e0 in (0.1, 0.45):
                for w0 in (0.5, 2.0, 3.6, 5.2):
                    starts.append([P0, 0.0, e0, w0, a0])
    else:
        starts = [p0]
    best = None
    for st in starts:
        st = list(st)
        st[0] = float(np.clip(st[0], MIN_P3 * 1.001, hi * 0.999))
        st[2] = float(np.clip(st[2], 0.0, 0.94))
        st[4] = float(np.clip(st[4], 1e-7, 0.99))
        try:
            r = least_squares(resid, st, bounds=([MIN_P3, -np.inf, 0, -np.inf, 0], [hi, np.inf, 0.95, np.inf, 1.0]),
                              x_scale=[st[0] * 0.05, st[0] * 0.1, 0.1, 0.5, st[4]], max_nfev=200 if search else 100)
        except Exception:  # noqa: BLE001
            continue
        c2 = float(r.fun @ r.fun)
        if best is None or c2 < best[1]:
            best = (r.x, c2)
    return best


def estimate_noise(d, model_fn):
    """Method-of-moments tau (sector common jitter) and s_w (white rescale) from the
    residuals of model_fn(d) -> residual vector (with tau=0, s_w=1)."""
    d.tau, d.sw = 0.0, 1.0
    r = model_fn(d)
    sig = d.sig0
    df = pd.DataFrame({"r": r, "w": 1 / sig ** 2, "s": d.df.sector.values})
    g = df.groupby("s")
    m = g.apply(lambda x: np.sum(x.w * x.r) / x.w.sum(), include_groups=False)
    vm = g.apply(lambda x: 1 / x.w.sum(), include_groups=False)
    # white rescale from within-sector scatter
    within = (df.r - df.s.map(m)) / sig
    sw = max(1.0, robust_sigma(within))
    tau2 = max(0.0, np.mean(m.values ** 2) - sw ** 2 * np.mean(vm.values))
    d.sw = float(sw)
    d.tau = float(np.sqrt(tau2)) if tau2 > 0 else 1e-7
    return d.tau, d.sw


def season_jackknife(d, base, p, gap=100.0):
    """Minimum over observing seasons (groups of eclipses separated by > gap days) of the
    dBIC of a common non-linear O-C (best of quadratic and LTTE) against the linear ephemeris,
    with that season removed (noise model fixed). Tests that the signal is not carried by one
    season; it does not ask the LTTE to beat a parabola, which a partial orbit cannot do."""
    t = d.df["T"].values
    o = np.argsort(t)
    season = np.zeros(len(t), int)
    season[o] = np.concatenate([[0], np.cumsum(np.diff(t[o]) > gap)])
    best = np.inf
    for s_ in np.unique(season):
        keep = season != s_
        if keep.sum() < 20:
            continue
        sub = d.df[keep].copy()
        dd = OC(sub, spot=d.spot_cols.shape[1] > 0)
        dd.tau, dd.sw = d.tau, d.sw
        n = len(dd.t)
        b0 = bic(dd.solve(dd.base("lin"), dd.y)[1], dd.nfix(dd.base("lin")), n)
        bq = bic(dd.solve(dd.base("quad"), dd.y)[1], dd.nfix(dd.base("quad")), n)
        pp = list(p)
        pp[1] = p[1] + d.tref - dd.tref           # T3 is relative to each table's tref
        f = fit_ltte(dd, base, p0=pp, search=False)
        if f is None:
            return -np.inf
        b1 = min(bic(f[1], dd.nfix(dd.base(base)) + 5, n), bq)
        best = min(best, b0 - b1)
    return float(best)


def analyse(df, spot=True, boot=B_BOOT, seed=0, P1=None):
    """Full analysis of one system's O-C table (P1 = eclipsing period, for the hierarchy
    check). Returns dict."""
    d = OC(df, spot=spot)
    n = len(d.t)
    # 1. preliminary LTTE (white noise only) to set the noise model from the best model
    res0 = {}
    for b in ("lin", "quad", "aps"):
        res0[b] = d.solve(d.base(b), d.y)[1]
    pl = {b: fit_ltte(d, b) for b in ("lin", "quad")}

    def resid_best(dd):
        cands = []
        for b in ("lin", "quad", "aps"):
            A = dd.base(b)
            _, c2, r, _ = dd.solve(A, dd.y)
            cands.append((bic(c2, dd.nfix(A), n), r))
        for b in ("lin", "quad"):
            A = dd.base(b)
            lt = el.ltte(dd.t, *pl[b][0])
            _, c2, r, _ = dd.solve(A, dd.y - lt)
            cands.append((bic(c2, dd.nfix(A) + 5, n), r))
        return min(cands, key=lambda x: x[0])[1]

    tau, sw = estimate_noise(d, resid_best)
    # 2. all models under the fixed noise model
    out = dict(n=n, n_sectors=len(d.sectors), span_d=d.span, tau_s=tau * 86400, white_scale=sw,
               n_spot_cols=int(d.spot_cols.shape[1]), has_sec=d.has_sec)
    bics = {}
    for b in ("lin", "quad", "aps"):
        A = d.base(b)
        bics[b] = bic(d.solve(A, d.y)[1], d.nfix(A), n)
    lt_fits = {}
    for b in ("lin", "quad"):
        p, c2 = fit_ltte(d, b)
        lt_fits[b] = (p, c2)
        bics["ltte_" + b] = bic(c2, d.nfix(d.base(b)) + 5, n)
    hi = MAX_P3_SPAN * d.span
    at_bound = {b: lt_fits[b][0][0] > 0.98 * hi for b in lt_fits}
    lb = "lin" if bics["ltte_lin"] <= bics["ltte_quad"] else "quad"
    if at_bound[lb] and not at_bound["quad" if lb == "lin" else "lin"]:
        lb = "quad" if lb == "lin" else "lin"
    p, c2 = lt_fits[lb]
    dbic = min(bics["lin"], bics["quad"], bics["aps"]) - bics["ltte_" + lb]
    # 3. wild cluster bootstrap
    A0 = d.base(lb)
    lt = el.ltte(d.t, *p)
    _, _, r_full, ranef = d.solve(A0, d.y - lt)
    model = d.y - r_full                      # fixed effects (base+spot) + LTTE
    resid = r_full                            # includes random sector offsets + white
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(boot):
        v = rng.choice([-1.0, 1.0], size=len(d.sectors))
        ystar = model + resid * (d.S @ v)
        f = fit_ltte(d, lb, y=ystar, p0=p, search=False)
        if f is not None:
            P3b, T3b, eb, wb, ab = f[0]
            draws.append([P3b, eb, el.amp_ltte(ab, eb, wb) * 86400, ab * 86400])
    draws = np.array(draws)

    def ci(col):
        if len(draws) < 10:
            return np.nan
        lo, hi_ = np.percentile(draws[:, col], [16, 84])
        return 0.5 * (hi_ - lo)

    P3, T3, e3, w3, asini = p
    A = el.amp_ltte(asini, e3, w3) * 86400
    sP, se, sA = ci(0), ci(1), ci(2)
    # pri/sec consistency at P3 (spot terms included)
    rho = ap_ = as_ = np.nan
    if d.has_sec:
        ph = 2 * np.pi * d.t / P3
        pri = (d.k > 0).astype(float)
        sec = 1 - pri
        X = np.column_stack([A0, np.sin(ph) * pri, np.cos(ph) * pri, np.sin(ph) * sec, np.cos(ph) * sec])
        c = d.solve(X, d.y)[0]
        vp, vs = c[-4:-2], c[-2:]
        rho = float(vp @ vs / (np.linalg.norm(vp) * np.linalg.norm(vs) + 1e-30))
        ap_, as_ = np.linalg.norm(vp) * 86400, np.linalg.norm(vs) * 86400
    # classification uses the raw bootstrap sigma; reported sigmas include the calibration
    snr = A / sA if sA and np.isfinite(sA) and sA > 0 else np.nan
    sP, se, sA = sP * CAL, se * CAL, sA * CAL
    flag = (dbic > DBIC_MIN) and (np.isfinite(snr) and snr > SNR_MIN) and not at_bound[lb]
    if d.has_sec and flag:
        ratio = ap_ / as_ if as_ > 0 else np.inf
        flag = flag and rho >= RHO_MIN and (1 / AMP_RATIO <= ratio <= AMP_RATIO)
    lab = "LTTE" if P3 / d.span <= LONG_SPAN else "LTTE-long"
    cls_v1 = lab if flag else "no-flag"           # thresholds fixed before the NCVZ test
    # v2 rules, added after vetting the NCVZ false-alarm set (disclosed in the README):
    #  (a) hierarchical: P3 >= MIN_P3_P1 * P1 (TIC 362104293: 'P3' = 29.5 d for P1 = 29.0 d)
    #  (b) A >= A_TAU_MIN * tau: the orbit must exceed the sector-to-sector jitter (spotted
    #      TIC 237116051 A/tau = 1.3, 199716496 A/tau = 2.3; all recovered B26 L1 >= 6.8)
    #  (c) season jackknife: common non-linear O-C still > DBIC_MIN with each season left out
    #      (TIC 377192658: signal carried by the first two sectors only)
    reasons = []
    if flag and P1 is not None and P3 < MIN_P3_P1 * P1:
        reasons.append("P3<5P1")
    if flag and A < A_TAU_MIN * tau * 86400:
        reasons.append("A<3tau")
    jk_min = np.nan
    if flag and not reasons:
        jk_min = season_jackknife(d, lb, p)
        if not (jk_min > DBIC_MIN):
            reasons.append("jackknife")
    cls = lab if (flag and not reasons) else "no-flag"
    out.update(cls_v1=cls_v1, v2_veto=";".join(reasons), jk_min_dBIC=float(jk_min))
    out.update(base=lb, dBIC=float(dbic), bics={k: float(v) for k, v in bics.items()},
               P3=float(P3), P3_sig=float(sP), e3=float(e3), e3_sig=float(se), A_s=float(A), A_sig=float(sA),
               asini_s=float(asini * 86400), f_m3=float(el.mass_function(asini, P3)),
               T3=float(T3 + d.tref), w3_deg=float(np.degrees(w3) % 360), P3_over_span=float(P3 / d.span),
               at_bound=bool(at_bound[lb]), snr=float(snr), rho=float(rho), A_pri_s=float(ap_), A_sec_s=float(as_),
               n_boot=int(len(draws)), boot_P3_med=float(np.median(draws[:, 0])) if len(draws) else np.nan,
               cls=cls, cal=CAL)
    return out
