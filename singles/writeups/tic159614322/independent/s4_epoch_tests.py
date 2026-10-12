"""For every predicted (non-discovery) epoch: plot all products and fit the depth of a fixed-shape transit
at the predicted time with a local quadratic baseline (linear least squares).  Also report results using
quality-flagged cadences.  Output: epoch_tests.json + fig_epoch_<tc>.png
"""
import json, os
import numpy as np, batman
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lcio import load, OUT

E = json.load(open(os.path.join(OUT, "epochs.json")))
FK = json.load(open(os.path.join(OUT, "fit_known.json")))["primary"]["best"]
P_REF = 50.0
REF_DEPTH = 0.0030
WIN = 1.0  # days either side for baseline fit


def shape(t, tc, cad):
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


def fit_depth(t, f, tc, cad, deg=2):
    """Linear LSQ: f = poly(t) - A*S(t). Returns A, sigma_A (scaled by rms), rms, chi2 improvement."""
    x = t - tc
    S = shape(t, tc, cad)
    if S.max() <= 0 or np.sum(S > 0.5) < 2:
        return None
    X = np.vstack([x ** i for i in range(deg + 1)] + [-S]).T
    beta, *_ = np.linalg.lstsq(X, f, rcond=None)
    r = f - X @ beta
    rms = np.std(r, ddof=X.shape[1])
    cov = np.linalg.inv(X.T @ X) * rms ** 2
    A, sA = beta[-1], np.sqrt(cov[-1, -1])
    X0 = X[:, :-1]
    b0, *_ = np.linalg.lstsq(X0, f, rcond=None)
    r0 = f - X0 @ b0
    dchi2 = (np.sum(r0 ** 2) - np.sum(r ** 2)) / rms ** 2
    # chi2 of the hypothesis "transit with reference depth present"
    fr = f + REF_DEPTH * S
    br, *_ = np.linalg.lstsq(X0, fr, rcond=None)
    chi2_ref = np.sum((fr - X0 @ br) ** 2) / rms ** 2
    return dict(A=float(A), sA=float(sA), rms=float(rms), dchi2=float(dchi2), npts=int(t.size),
                n_in=int(np.sum(S > 0.05)), chi2_ref_minus_best=float(chi2_ref - np.sum(r ** 2) / rms ** 2))


def get(sector, product, tc, which="flux", strict=True):
    d = load(sector, product, which, strict_quality=strict)
    if d is None:
        return None
    m = np.abs(d.time - tc) < WIN
    if m.sum() < 10:
        return None
    t, f = d.time[m].values, d.f[m].values
    f = f / np.median(f)
    return t, f, d.cad.iloc[0], d[m]


uniq = {}
for r in E["epochs"]:
    if r["known"]:
        continue
    key = round(r["tc"], 2)
    uniq.setdefault(key, dict(tc=r["tc"], sig_tc=r["sig_tc"], sector=r["sector"], ns=[]))
    uniq[key]["ns"].append(r["n"])

results = []
for key, u in sorted(uniq.items()):
    tc, s = u["tc"], u["sector"]
    res = dict(tc=tc, sig_tc=u["sig_tc"], sector=s, ns=sorted(set(u["ns"])), fits={})
    fig, axs = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    off = 0
    for p, col in [("SPOC", "C0"), ("TESS-SPOC", "C1"), ("QLP", "C2")]:
        for which in ["flux", "sap"]:
            for strict in [True, False]:
                g = get(s, p, tc, which, strict)
                if g is None:
                    continue
                t, f, cad, d = g
                fd = fit_depth(t, f, tc, cad)
                res["fits"][f"{p}_{which}_{'q0' if strict else 'allq'}"] = fd
        # plot: all-quality flux in grey, q=0 in colour
        g_all = get(s, p, tc, "flux", False)
        g = get(s, p, tc, "flux", True)
        ax = axs[0] if p != "QLP" else axs[1]
        if p == "QLP":
            off = 0
        if g_all is not None:
            t, f, cad, d = g_all
            if p == "SPOC":
                tb = np.arange(tc - WIN, tc + WIN, 1 / 48)
                idx = np.digitize(t, tb)
                t = np.array([t[idx == k].mean() for k in np.unique(idx)])
                f = np.array([f[idx == k].mean() for k in np.unique(idx)])
            ax.plot(t - tc, f + off, "x", color="0.6", ms=3)
        if g is not None:
            t, f, cad, d = g
            tt = t
            if p == "SPOC":
                tb = np.arange(tc - WIN, tc + WIN, 1 / 48)
                idx = np.digitize(t, tb)
                tt = np.array([t[idx == k].mean() for k in np.unique(idx)])
                f = np.array([f[idx == k].mean() for k in np.unique(idx)])
            ax.plot(tt - tc, f + off, ".", color=col, ms=4, label=f"{p} (q=0)")
            # background and centroid
            bk = d["sap_bkg"].values
            axs[2].plot(t - tc if p != "SPOC" else d.time.values - tc, bk / np.nanmedian(np.abs(bk)), ".", ms=2, color=col, label=f"{p} bkg (norm)")
            if "mom_centr1" in d:
                axs[3].plot(d.time.values - tc, d["mom_centr1"].values - np.nanmedian(d["mom_centr1"]), ".", ms=2, color=col, label=f"{p} col")
                axs[3].plot(d.time.values - tc, d["mom_centr2"].values - np.nanmedian(d["mom_centr2"]) + 0.05, "+", ms=2, color=col, label=f"{p} row+0.05")
        mod = 1 - REF_DEPTH * shape(np.linspace(tc - WIN, tc + WIN, 800), tc, 120)
        ax.plot(np.linspace(-WIN, WIN, 800), mod + off, "k-", lw=0.6, alpha=0.6)
        if p != "QLP":
            off -= 0.006
    for ax in axs[:2]:
        ax.axvspan(-E["T14_d"] / 2, E["T14_d"] / 2, color="r", alpha=0.08)
        ax.legend(fontsize=7)
    axs[0].set_ylim(0.983, 1.004); axs[1].set_ylim(0.993, 1.004); axs[0].set_ylim(0.983, 1.004); axs[2].set_ylim(-3, 4)
    axs[2].legend(fontsize=7); axs[3].legend(fontsize=7)
    axs[3].set_xlabel(f"t - {tc:.3f} [d]")
    axs[0].set_title(f"Predicted epoch {tc:.3f} (S{s}), n = {res['ns']}; black = 3.0 ppt model; grey x = flagged cadences")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, f"fig_epoch_{tc:.2f}.png"), dpi=100)
    plt.close()
    results.append(res)
    print(f"tc={tc:.3f} S{s} n={res['ns']}")
    for k, v in res["fits"].items():
        if v:
            print(f"   {k:24s} A={v['A']*1e3:6.2f}±{v['sA']*1e3:5.2f} ppt  n_in={v['n_in']:4d}  dchi2={v['dchi2']:7.1f} chi2(3ppt)-chi2(best)={v['chi2_ref_minus_best']:7.1f}")
        else:
            print(f"   {k:24s} no in-transit data")

json.dump(results, open(os.path.join(OUT, "epoch_tests.json"), "w"), indent=1)
