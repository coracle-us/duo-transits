"""Fit individual transit mid-times.

For a star with a trial ephemeris (P0, T00 from the search), every predicted epoch that
falls in data is examined. Timing source per sector (brief's rules):
  TESS-SPOC if present; else QLP for sectors < 94; TESScut for sectors >= 94 (and for
  sectors without any LC).  Other sources are also fitted for comparison only.

Method
  1. Each sector is detrended with a wotan biweight (window 3x duration, >= 0.75 d) with all
     predicted transits masked (+-0.9 duration).
  2. A joint batman fit over the timing-source transits gives a common shape
     (rp, a/R*, b; quadratic LD fixed at (0.40, 0.25)); P fixed at P0 (sets the shape only).
  3. Each transit is then fitted for (Tc, flux offset, slope, depth scale) with the shape
     otherwise fixed, in a window of +-1.5 durations. Errors: least-squares covariance with
     the flux errors rescaled to the window's out-of-transit scatter, times a beta factor
     for time-correlated noise (binned-rms method), and checked with a delta-chi2 = 1
     profile.  Coverage of the transit is recorded (ingress/egress present?).
Output: out/times_<tic>.csv and a per-transit plot out/transits_<tic>.png
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import TARGETS, OUT, load_lc_files, load_tesscut_phot, detrend, batman_model

U = (0.40, 0.25)


def timing_source_ok(lc, have_spoc):
    s, src = lc["sector"], lc["source"]
    if src == "TESS-SPOC":
        return True
    if src == "QLP":
        return (s < 94) and not have_spoc
    if src == "TESScut":
        return True  # only exists for S>=94 or sectors without LC (see s02)
    return False


def beta_factor(res, cad, dur):
    bs = []
    for tb in [dur / 8, dur / 4, dur / 2]:
        nb = int(round(tb / cad))
        if nb < 2:
            continue
        k = (len(res) // nb) * nb
        if k < 4 * nb:
            continue
        s1 = np.std(res)
        sb = np.std(res[:k].reshape(-1, nb).mean(1))
        bs.append(sb / (s1 / np.sqrt(nb)) * np.sqrt((k // nb) / (k // nb - 1)))
    return max(1.0, np.median(bs)) if bs else 1.0


def model(theta_shape, tc, t, cad, P):
    rp, a, b = theta_shape
    inc = np.degrees(np.arccos(np.clip(b / a, 0, 1)))
    return batman_model(t, tc, P, rp, a, inc, U, exptime=cad)


def gather(tic, P0, T00):
    tgt = TARGETS[tic]
    dur = tgt["dur"]
    lcs = load_lc_files(tic) + load_tesscut_phot(tic)
    spoc_secs = {lc["sector"] for lc in lcs if lc["source"] == "TESS-SPOC"}
    segs = []
    for lc in lcs:
        t, f, e = lc["t"], lc["f"], lc["e"]
        if len(t) < 100:
            continue
        n = np.arange(np.ceil((t.min() - T00 - dur) / P0), np.floor((t.max() - T00 + dur) / P0) + 1)
        pred = T00 + n * P0
        mask = np.zeros(len(t), bool)
        for tp in pred:
            mask |= np.abs(t - tp) < 0.9 * dur
        flat, _ = detrend(t, f, max(0.75, 3 * dur), mask=mask)
        cad = np.nanmedian(np.diff(t))
        for ep, tp in zip(n.astype(int), pred):
            w = np.abs(t - tp) < 1.5 * dur
            inn = np.abs(t - tp) < 0.5 * dur
            cov_in = inn.sum() * cad / dur
            # ingress/egress coverage
            ing = (np.abs(t - (tp - dur / 2)) < 0.15 * dur).sum() * cad / (0.3 * dur)
            egr = (np.abs(t - (tp + dur / 2)) < 0.15 * dur).sum() * cad / (0.3 * dur)
            segs.append(dict(sector=lc["sector"], source=lc["source"], epoch=ep, tpred=tp,
                             t=t[w], f=flat[w], e=e[w], cad=cad, cov=cov_in, ing=ing, egr=egr,
                             timing=timing_source_ok(lc, lc["sector"] in spoc_secs)))
    return segs


def joint_shape(segs, P0, dur, depth):
    use = [s for s in segs if s["timing"] and s["cov"] > 0.7 and s["ing"] > 0.5 and s["egr"] > 0.5]
    k0 = np.sqrt(depth)
    a0 = P0 / (np.pi * dur)  # b=0 circular approx
    x0 = [k0, a0, 0.3] + [s["tpred"] for s in use] + [0.0] * len(use)

    def resid(x):
        rp, a, b = x[:3]
        out = []
        for i, s in enumerate(use):
            m = model((rp, a, b), x[3 + i], s["t"], s["cad"], P0) * (1 + x[3 + len(use) + i])
            out.append((s["f"] - m) / (s["e"]))
        return np.concatenate(out)

    lo = [0.01, 1.5, 0.0] + [s["tpred"] - 0.3 * dur for s in use] + [-0.01] * len(use)
    hi = [0.5, 500, 1.2] + [s["tpred"] + 0.3 * dur for s in use] + [0.01] * len(use)
    best = None
    for b0 in [0.1, 0.5, 0.8]:
        x0[2] = b0
        r = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac")
        if best is None or r.cost < best.cost:
            best = r
    return best.x[:3], len(use)


def fit_one(s, shape, P0, dur):
    t, f, e, cad = s["t"], s["f"], s["e"], s["cad"]
    if len(t) < 20 or s["cov"] < 0.3:
        return None

    def resid(x, ee):
        tc, off, slope, scl = x
        m0 = model(shape, tc, t, cad, P0)
        m = (1 - scl * (1 - m0)) * (1 + off + slope * (t - s["tpred"]))
        return (f - m) / ee

    x0 = [s["tpred"], 0, 0, 1.0]
    lo = [s["tpred"] - 0.4 * dur, -0.02, -0.1, 0.1]
    hi = [s["tpred"] + 0.4 * dur, 0.02, 0.1, 3.0]
    # coarse grid in tc to avoid local minima
    grid = np.linspace(lo[0], hi[0], 41)
    c = [np.sum(resid([g, 0, 0, 1.0], e) ** 2) for g in grid]
    x0[0] = grid[int(np.argmin(c))]
    r = least_squares(resid, x0, bounds=(lo, hi), args=(e,), x_scale="jac")
    # rescale errors to out-of-transit scatter
    res = r.fun * e
    oot = np.abs(t - r.x[0]) > 0.6 * dur
    sig = np.std(res[oot]) if oot.sum() > 10 else np.std(res)
    ee = np.full(len(t), sig)
    r = least_squares(resid, r.x, bounds=(lo, hi), args=(ee,), x_scale="jac")
    J = r.jac
    cov = np.linalg.pinv(J.T @ J)
    err = np.sqrt(np.diag(cov))
    beta = beta_factor(r.fun * sig, cad, dur)
    # chi2 profile on tc (other params re-optimised)
    tcs = r.x[0] + np.linspace(-4, 4, 41) * err[0] * 1.0
    prof = []
    for tc in tcs:
        def rr(y):
            return resid([tc, *y], ee)
        q = least_squares(rr, r.x[1:], bounds=(lo[1:], hi[1:]))
        prof.append(2 * q.cost)
    prof = np.array(prof) - 2 * r.cost
    within = tcs[prof <= 1.0]
    prof_err = (within.max() - within.min()) / 2 if len(within) > 1 else err[0]
    tc_err = max(err[0], prof_err) * beta
    rchi2 = 2 * r.cost / (len(t) - 4)
    return dict(tc=r.x[0], tc_err=tc_err, tc_err_white=err[0], prof_err=prof_err, beta=beta,
                depth_scale=r.x[3], depth_scale_err=err[3] * beta, sig=sig, rchi2=rchi2, x=r.x)


def run(tic, P0, T00):
    tgt = TARGETS[tic]
    dur, dep = tgt["dur"], tgt["depth"]
    segs = gather(tic, P0, T00)
    shape, nused = joint_shape(segs, P0, dur, dep)
    rp, a, b = shape
    dep_model = 1 - model(shape, 0.0, np.array([0.0]), 0.0, P0)[0]
    print(f"TIC {tic}: joint shape from {nused} full transits: rp={rp:.4f} a/R*={a:.2f} b={b:.2f} "
          f"model depth={dep_model*1e6:.0f} ppm")
    rows = []
    for s in segs:
        if s["cov"] < 0.05:
            continue
        r = fit_one(s, shape, P0, dur)
        rows.append((s, r))
    # Per-epoch timing source: the first usable of TESS-SPOC > QLP (sector < 94) > TESScut,
    # where usable = fitted, in-transit coverage >= 0.4, ingress or egress covered, a sane fit
    # (Tc error < 0.1 duration, depth scale > 0.3).
    # (So QLP supplies a transit that SPOC's cadence mask drops, but never from S94 on.)
    prio = {"TESS-SPOC": 0, "QLP": 1, "TESScut": 2}
    for s, r in rows:
        s["timing"] = False
    for ep in sorted({s["epoch"] for s, _ in rows}):
        cands = [(prio[s["source"]], i) for i, (s, r) in enumerate(rows)
                 if s["epoch"] == ep and r is not None and s["cov"] >= 0.4 and max(s["ing"], s["egr"]) >= 0.5
                 and r["tc_err"] < 0.1 * dur and r["depth_scale"] > 0.3
                 and not (s["source"] == "QLP" and s["sector"] >= 94)]
        if cands:
            rows[min(cands)[1]][0]["timing"] = True
    # write
    fn = os.path.join(OUT, f"times_{tic}.csv")
    with open(fn, "w") as fh:
        fh.write("sector,source,epoch,timing_source,cov_in,ingress_cov,egress_cov,tc,tc_err,tc_err_white,"
                 "prof_err,beta,depth_scale,depth_scale_err,depth_ppm,oot_rms_ppm\n")
        for s, r in rows:
            if r is None:
                fh.write(f"{s['sector']},{s['source']},{s['epoch']},{int(s['timing'])},{s['cov']:.2f},"
                         f"{s['ing']:.2f},{s['egr']:.2f},,,,,,,,,\n")
                continue
            fh.write(f"{s['sector']},{s['source']},{s['epoch']},{int(s['timing'])},{s['cov']:.2f},"
                     f"{s['ing']:.2f},{s['egr']:.2f},{r['tc']:.6f},{r['tc_err']:.6f},{r['tc_err_white']:.6f},"
                     f"{r['prof_err']:.6f},{r['beta']:.2f},{r['depth_scale']:.3f},{r['depth_scale_err']:.3f},"
                     f"{r['depth_scale']*dep_model*1e6:.0f},{r['sig']*1e6:.0f}\n")
    json.dump(dict(rp=rp, a=a, b=b, P0=P0, T00=T00, depth=dep_model, u=U),
              open(os.path.join(OUT, f"shape_{tic}.json"), "w"), indent=1)
    # plot
    ok = [(s, r) for s, r in rows if r is not None]
    nc = 4
    nr = int(np.ceil(len(ok) / nc)) or 1
    fig, axs = plt.subplots(nr, nc, figsize=(16, 3 * nr), squeeze=False)
    for ax, (s, r) in zip(axs.flat, ok):
        t = s["t"]
        ax.plot((t - r["tc"]) * 24, s["f"], ".", ms=1.5, color="0.6")
        nb = max(1, int(round((20 / 1440) / s["cad"])))
        k = (len(t) // nb) * nb
        if k:
            ax.plot((t[:k].reshape(-1, nb).mean(1) - r["tc"]) * 24, s["f"][:k].reshape(-1, nb).mean(1), "k.", ms=3)
        tt = np.linspace(t.min(), t.max(), 600)
        m0 = model(shape, r["tc"], tt, s["cad"], P0)
        tc, off, slope, scl = r["x"]
        ax.plot((tt - tc) * 24, (1 - scl * (1 - m0)) * (1 + off + slope * (tt - s["tpred"])), "r-", lw=1)
        ax.set_title(f"S{s['sector']} {s['source']} E{s['epoch']}{' [T]' if s['timing'] else ''} "
                     f"Tc={r['tc']:.4f}±{r['tc_err']*1440:.1f}m d×{r['depth_scale']:.2f}", fontsize=8)
        ax.set_ylim(1 - 2 * dep, 1 + 0.8 * dep)
    for ax in list(axs.flat)[len(ok):]:
        ax.axis("off")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"transits_{tic}.png"), dpi=70); plt.close(fig)
    print(open(fn).read())


if __name__ == "__main__":
    tic, P0, T00 = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    run(tic, P0, T00)
