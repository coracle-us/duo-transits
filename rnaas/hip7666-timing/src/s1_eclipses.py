"""Binary template, pulsation subtraction and eclipse times (RULES.md sections 2-3).

Usage: python src/s1_eclipses.py [main|qlp|sap] [nboot]   (sap = S58 SAP_FLUX check)
Writes results/template_S58.npz, results/pulsation_fit_<kind>.csv, results/eclipse_times_<kind>.csv
and results/dropped_eclipses_<kind>.csv.
"""
import sys
import csv
import pickle

import numpy as np

import hip

KIND = sys.argv[1] if len(sys.argv) > 1 else "main"
NBOOT = int(sys.argv[2]) if len(sys.argv) > 2 else 500
FN = list(hip.FREQS.values())
SECTORS = [58] if KIND == "sap" else [18, 58, 85]


def build_template(t, y):
    tpl = hip.Template(hip.phase(t), y)
    tpl.recentre()
    return tpl


def prepare(sector, tpl=None, kind=KIND):
    """Return dict with t, flux, err, trend, pulsation model, corrected flux y, template."""
    t, f, e, exptime = hip.load(kind, sector)
    oot = ~hip.in_eclipse(t)
    trend = hip.trend_spline(t, f, oot, spacing=5.0)               # pass 0: slow trends only
    own = tpl is None
    if own:
        tpl = build_template(t, f / trend)                         # pass 1
    puls = np.zeros_like(t)
    for _ in range(2):                                             # passes 2-3
        shape = hip.smear(tpl, hip.phase(t), 0.0) if exptime <= 130 else \
            hip.smear(lambda tt: tpl(hip.phase(tt)), t, exptime)
        trend = hip.trend_spline(t, (f - puls * trend) / shape, oot, spacing=1.5)
        r = f / trend - shape               # pulsation is additive at full amplitude (RULES 2)
        coef, cov, rms = hip.fit_sines(t[oot], r[oot], FN)
        puls = hip.sines_model(t, FN, coef)
        if own:
            tpl = build_template(t, f / trend - puls)
    y = f / trend - puls
    return dict(t=t, f=f, e=e, exptime=exptime, trend=trend, puls=puls, y=y, tpl=tpl,
                coef=coef, cov=cov, rms=rms, oot=oot)


def eclipse_list(t, tpl):
    e0 = int(np.floor((t.min() - hip.T0_LEE) / hip.P_LEE)) - 1
    e1 = int(np.ceil((t.max() - hip.T0_LEE) / hip.P_LEE)) + 1
    out = []
    for E in range(e0, e1 + 1):
        for which, off in (("p", 0.0), ("s", 0.5)):
            tc = hip.T0_LEE + (E + off + tpl.shift[which]) * hip.P_LEE
            out.append((E, which, tc))
    return out


def time_eclipses(sector, d, writer, dropped, nboot=NBOOT):
    t, y, e, exptime, tpl = d["t"], d["y"], d["e"], d["exptime"], d["tpl"]
    nexp = 2 * hip.ECL_HALF * hip.P_LEE * 86400 / exptime
    for E, which, tc in eclipse_list(t, tpl):
        dph = (t - tc) / hip.P_LEE
        win = np.abs(dph) < hip.FIT_HALF
        core = np.abs(dph) < hip.ECL_HALF
        if core.sum() == 0:
            continue
        ing = ((dph > -hip.ECL_HALF) & (dph < -0.015)).sum()
        egr = ((dph < hip.ECL_HALF) & (dph > 0.015)).sum()
        nside = (hip.ECL_HALF - 0.015) * hip.P_LEE * 86400 / exptime
        cov = core.sum() / nexp
        if cov < 0.8 or ing < 0.5 * nside or egr < 0.5 * nside:
            dropped.writerow([sector, which, E, f"{tc:.5f}", f"{cov:.2f}", ing, egr])
            continue
        prof = tpl.eclipse_profile(which)
        tc_ref = hip.T0_LEE + (E + (0.0 if which == "p" else 0.5) + tpl.shift[which]) * hip.P_LEE
        best, r = hip.fit_eclipse(t[win], y[win], e[win], tc_ref, prof, exptime)
        sig = hip.block_bootstrap(t[win], y[win], e[win], tc_ref, prof, exptime, best,
                                  nboot=nboot, seed=1000 * sector + 2 * E + (which == "s"))
        tobs = tc_ref + best[0]
        res = y[win] - hip.eclipse_model(best, t[win], tc_ref, prof, exptime)
        # O-C against Lee Eq. 1 (secondary expected at phase 0.5)
        tlee = hip.T0_LEE + (E + (0.0 if which == "p" else 0.5)) * hip.P_LEE
        writer.writerow([sector, which, E, f"{tobs:.7f}", f"{sig:.7f}", f"{best[1]:.4f}",
                         int(win.sum()), f"{np.std(res):.6f}", f"{(tobs - tlee) * 86400:.2f}"])


def main():
    # S58 main product defines the template for every sector and product
    with open(hip.OUTDIR / "template_S58.pkl", "rb") if KIND != "main" else open(
            hip.OUTDIR / "template_S58.pkl", "wb") as fh:
        if KIND == "main":
            d58 = prepare(58, kind="main")
            pickle.dump(d58["tpl"], fh)
            tpl = d58["tpl"]
        else:
            tpl = pickle.load(fh)
    print("template recentring (phase):", tpl.shift,
          "-> secondary offset from 0.5 P: %.1f s" % ((tpl.shift["s"] - tpl.shift["p"]) * hip.P_LEE * 86400))
    fe = open(hip.OUTDIR / f"eclipse_times_{KIND}.csv", "w", newline="")
    fd = open(hip.OUTDIR / f"dropped_eclipses_{KIND}.csv", "w", newline="")
    fp = open(hip.OUTDIR / f"pulsation_fit_{KIND}.csv", "w", newline="")
    we, wd, wp = csv.writer(fe), csv.writer(fd), csv.writer(fp)
    we.writerow(["sector", "type", "E", "t_btjd", "err_d", "depth_scale", "n", "rms", "oc_lee_s"])
    wd.writerow(["sector", "type", "E", "tc_pred", "coverage", "n_ingress", "n_egress"])
    wp.writerow(["sector", "mode", "freq", "amp_ppt", "phase_rad", "phase_err_rad", "rms_ppt"])
    prepared = {}
    for s in SECTORS:
        d = prepare(s, tpl=tpl, kind=KIND)
        prepared[s] = d
        for i, (name, fq) in enumerate(hip.FREQS.items()):
            a, b = d["coef"][2 * i], d["coef"][2 * i + 1]
            A, ph = hip.amp_phase(a, b)
            va, vb, cab = d["cov"][2 * i, 2 * i], d["cov"][2 * i + 1, 2 * i + 1], d["cov"][2 * i, 2 * i + 1]
            dphi = np.sqrt((b ** 2 * va + a ** 2 * vb - 2 * a * b * cab)) / A ** 2
            wp.writerow([s, name, fq, f"{A * 1e3:.4f}", f"{ph:.4f}", f"{dphi:.4f}", f"{d['rms'] * 1e3:.4f}"])
        time_eclipses(s, d, we, wd)
        fe.flush()
        print(f"S{s} {KIND}: {len(d['t'])} cadences, exptime {d['exptime']:.0f} s, oot rms {d['rms']*1e3:.3f} ppt")
    for fh in (fe, fd, fp):
        fh.close()
    with open(hip.OUTDIR / f"prepared_{KIND}.pkl", "wb") as fh:
        pickle.dump({s: {k: v for k, v in d.items() if k != "tpl"} for s, d in prepared.items()}, fh)


if __name__ == "__main__":
    main()
