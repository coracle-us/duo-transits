"""Alias analysis.

1. Event mid-times (from the joint 4-event fit) -> every period P = (T_last - T_first)/N, N=1..Nmax
   (P >= PMIN d), tested with a weighted linear ephemeris; keep P with chi2 < 13.8 (99.9%, 2 dof).
2. For each surviving P: every predicted epoch inside every dataset is tested with a fixed-shape
   template (the joint-fit transit model) whose amplitude A is fitted together with a local quadratic
   baseline (+-0.6 d). Noise on A is empirical: the same fit at ~150 random epochs in the dataset
   (away from the known events) -> robust std (includes red noise).
   Coverage = fraction of the T14 window that has cadences.
   Status per epoch/dataset:  'detected'  A>0.5 and A/sigma>5
                              'excluded'  coverage>=0.6, A<0.5 and (1-A)/sigma>5
                              'gap'       coverage<0.2 ;  'partial'/'ambiguous' otherwise.
   A period is EXCLUDED if any epoch is 'excluded' in >=2 reductions (or in the only reduction available)
   and no reduction calls it 'detected'.
"""
import os, json, sys
import numpy as np
from astropy.io import fits
import lcio, tfit

HERE = lcio.HERE
PMIN = 1.0
fr = json.load(open(os.path.join(HERE, "fit_results.json")))
J = fr["joint4_primary"]
T = np.array([x[0] for x in J["tc"]])
S = np.array([0.5 * (x[1] + x[2]) for x in J["tc"]])
ch = np.load(os.path.join(HERE, "cache", "chain_joint4_primary.npy"))
k_, b_, aR_ = np.median(ch[:, 0]), np.median(ch[:, 1]), np.exp(np.median(ch[:, 2]))
T14 = tfit.t14(tfit.PREF, aR_, k_, b_)
print("events", T, "+-", S, " template k=%.4f b=%.3f aR(Pref)=%.2f T14=%.3f h" % (k_, b_, aR_, T14 * 24))

# ---------------- 1. ephemeris-consistent periods
span = T[-1] - T[0]
surv = []
for N in range(1, int(span / PMIN) + 1):
    P = span / N
    E = np.round((T - T[0]) / P)
    A = np.vstack([np.ones_like(E), E]).T
    Wt = 1 / S ** 2
    cov = np.linalg.inv(A.T @ (A * Wt[:, None]))
    sol = cov @ (A.T @ (Wt * T))
    chi2 = np.sum(((T - A @ sol) / S) ** 2)
    if chi2 < 13.8:
        surv.append(dict(N=N, P=float(sol[1]), sP=float(np.sqrt(cov[1, 1])), T0=float(sol[0]), sT0=float(np.sqrt(cov[0, 0])), chi2=float(chi2), epochs=E.astype(int).tolist()))
print(f"{len(surv)} ephemeris-consistent periods with P >= {PMIN} d; N values:", [s["N"] for s in surv][:40], "...")
P1 = surv[0]["P"]
print("longest: P=%.5f +- %.5f d (N=%d, chi2=%.2f, epochs %s)" % (surv[0]["P"], surv[0]["sP"], surv[0]["N"], surv[0]["chi2"], surv[0]["epochs"]))
for s in surv:
    s["m"] = P1 / s["P"]

# ---------------- 2. datasets
HARD = 1 | 2 | 4 | 8 | 32 | 128 | (1 << 20)


def ev_noq(c):
    d = fits.getdata(lcio._find(f"hlsp_everest_k2_llc_211924561-c{c:02d}_kepler_v2.0_lc.fits"), 1)
    return lcio._norm(d["TIME"] - 2167.0, np.asarray(d["FLUX"], float), None, (d["QUALITY"] & (1 << 20)) == 0)


def k2sap(c):
    d = fits.getdata(lcio._find(f"ktwo211924561-c{c:02d}_llc.fits"), 1)
    return lcio._norm(d["TIME"] - 2167.0, np.asarray(d["SAP_FLUX"], float), None, (d["SAP_QUALITY"] & HARD) == 0)


KEXP = 0.020434
DS = {  # name: (loader, family, exptime, LD)
    "C5_K2SFF": (lambda: lcio.k2_k2sff(5), "K2C5", KEXP, (0.31, 0.31)),
    "C5_EVEREST": (lambda: ev_noq(5), "K2C5", KEXP, (0.31, 0.31)),
    "C16_EVEREST": (lambda: ev_noq(16), "K2C16", KEXP, (0.31, 0.31)),
    "C16_SAP": (lambda: k2sap(16), "K2C16", KEXP, (0.31, 0.31)),
    "C16_K2SFF": (lambda: lcio.k2_k2sff(16), "K2C16", KEXP, (0.31, 0.31)),
    "C18_EVEREST": (lambda: ev_noq(18), "K2C18", KEXP, (0.31, 0.31)),
    "C18_K2SFF": (lambda: lcio.k2_k2sff(18), "K2C18", KEXP, (0.31, 0.31)),
    "C18_SAP": (lambda: k2sap(18), "K2C18", KEXP, (0.31, 0.31)),
}
for s_ in (44, 45, 46, 72):
    DS[f"S{s_}_SPOC2min"] = (lambda s_=s_: lcio.tess_spoc2min(s_), f"S{s_}", None, (0.23, 0.31))
    DS[f"S{s_}_QLP"] = (lambda s_=s_: lcio.tess_qlp(s_), f"S{s_}", 600 / 86400. if s_ < 56 else 200 / 86400., (0.23, 0.31))
    DS[f"S{s_}_TESSSPOC"] = (lambda s_=s_: lcio.tess_tessspoc(s_), f"S{s_}", 600 / 86400. if s_ < 56 else 200 / 86400., (0.23, 0.31))
DS["S1751_TESSSPOC"] = (lambda: lcio.tess_tessspoc(1751), "S1751", 200 / 86400., (0.23, 0.31))
DS["S1751_SAP"] = (lambda: lcio.tess_tessspoc(1751, col="SAP_FLUX"), "S1751", 200 / 86400., (0.23, 0.31))

HWIN = 0.6
rng = np.random.default_rng(7)


def fit_amp(t, f, tc, exp, u):
    w = np.abs(t - tc) < HWIN
    if w.sum() < 10:
        return np.nan, 0.0
    tt, ff = t[w], f[w]
    m = tfit.model(tt, tc, k_, b_, aR_, u, exp) - 1  # negative in transit
    x = tt - tc
    A = np.vstack([np.ones_like(x), x, x ** 2, m]).T
    c, *_ = np.linalg.lstsq(A, ff, rcond=None)
    amp = c[3] / c[0]
    # coverage of the T14 window
    cad = np.median(np.diff(tt)) if len(tt) > 1 else 0.02
    grid = np.linspace(tc - T14 / 2, tc + T14 / 2, 50)
    cov = np.mean([np.any(np.abs(tt - g) < 0.75 * max(cad, 2 / 1440.)) for g in grid])
    return amp, cov


cache = {}
for name, (ld, fam, exp, u) in DS.items():
    try:
        d = ld()
    except FileNotFoundError:
        print("missing", name); continue
    t, f = d["t"], d["f"]
    o = np.argsort(t); t, f = t[o], f[o]
    # remove 6-sigma single-point upward outliers only
    r = f / np.convolve(f, np.ones(9) / 9, mode="same") - 1
    s = 1.4826 * np.median(np.abs(r))
    keep = r < 6 * s
    t, f = t[keep], f[keep]
    # empirical sigma of A at random epochs away from known events
    amps = []
    tries = 0
    while len(amps) < 150 and tries < 2000:
        tries += 1
        tc = rng.uniform(t.min() + 0.3, t.max() - 0.3)
        if np.min(np.abs(tc - T)) < 0.5:
            continue
        a, cv = fit_amp(t, f, tc, exp, u)
        if cv > 0.8 and np.isfinite(a):
            amps.append(a)
    amps = np.array(amps)
    sig = 1.4826 * np.median(np.abs(amps - np.median(amps))) if len(amps) > 10 else np.nan
    cache[name] = dict(t=t, f=f, exp=exp, u=u, fam=fam, sigA=sig, tmin=t.min(), tmax=t.max())
    print(f"{name:16s} span {t.min():8.2f}-{t.max():8.2f}  sigma_A(empirical)={sig:.3f}  (A=1 <-> full {1e3*(1-tfit.model(np.array([0.]),0,k_,b_,aR_,u)[0]):.2f} ppt)")

# ---------------- 3. test each period
table = []
for s in surv:
    P, T0 = s["P"], s["T0"]
    recs = []
    for name, c in cache.items():
        n0 = int(np.ceil((c["tmin"] - 0.2 - T0) / P)); n1 = int(np.floor((c["tmax"] + 0.2 - T0) / P))
        for n in range(n0, n1 + 1):
            tc = T0 + n * P
            if np.min(np.abs(tc - T)) < 0.05:
                known = True
            else:
                known = False
            a, cv = fit_amp(c["t"], c["f"], tc, c["exp"], c["u"])
            if not np.isfinite(a) or cv < 0.2:
                st = "gap"
            elif a > 0.5 and a / c["sigA"] > 5:
                st = "detected"
            elif cv >= 0.6 and a < 0.5 and (1 - a) / c["sigA"] > 5:
                st = "excluded"
            else:
                st = "ambiguous"
            if st != "gap" or cv > 0:
                recs.append(dict(ds=name, fam=c["fam"], epoch=n, tc=round(tc, 4), A=None if not np.isfinite(a) else round(float(a), 3),
                                 sigA=round(float(c["sigA"]), 3), cov=round(float(cv), 2), status=st, known=known))
    # aggregate per (fam, epoch)
    groups = {}
    for r in recs:
        groups.setdefault((r["fam"], r["epoch"]), []).append(r)
    excl_epochs, det_new, amb = [], [], []
    for key, rs in groups.items():
        sts = [r["status"] for r in rs]
        nex = sts.count("excluded"); nd = sts.count("detected"); nvalid = sum(st != "gap" for st in sts)
        if any(r["known"] for r in rs):
            continue
        if nd == 0 and nex >= 1 and (nex >= 2 or nvalid == 1):
            excl_epochs.append(dict(fam=key[0], tc=rs[0]["tc"], n_excl=nex, n_valid=nvalid, A=[r["A"] for r in rs], ds=[r["ds"] for r in rs]))
        elif nd >= 1:
            det_new.append(dict(fam=key[0], tc=rs[0]["tc"], A=[r["A"] for r in rs], ds=[r["ds"] for r in rs]))
        elif nvalid > 0:
            amb.append(dict(fam=key[0], tc=rs[0]["tc"], A=[r["A"] for r in rs], cov=[r["cov"] for r in rs], st=sts, ds=[r["ds"] for r in rs]))
    s.update(status="excluded" if excl_epochs else "survives", excluding_epochs=excl_epochs, new_detections=det_new, ambiguous=amb,
             n_epochs_tested=len(groups))
    table.append(s)

surv_final = [s for s in table if s["status"] == "survives"]
print("\nSURVIVING PERIODS:")
for s in surv_final:
    print("  m=%.3f N=%d P=%.5f +- %.5f d  T0=%.5f  ambiguous epochs: %s  new detections: %s" % (s["m"], s["N"], s["P"], s["sP"], s["T0"],
          [(a["fam"], a["tc"], a["A"], a["cov"]) for a in s["ambiguous"]], s["new_detections"]))
print("\nExcluded (first 25):")
for s in [s for s in table if s["status"] == "excluded"][:25]:
    e = s["excluding_epochs"]
    print("  m=%.2f P=%.4f : %d excluding epochs, e.g. %s" % (s["m"], s["P"], len(e), [(x["fam"], x["tc"], x["A"]) for x in e[:3]]))
json.dump(dict(events=T.tolist(), event_err=S.tolist(), template=dict(k=k_, b=b_, aR_at_Pref=aR_, T14_h=T14 * 24),
               periods=table), open(os.path.join(HERE, "alias_table.json"), "w"), indent=1, default=float)
