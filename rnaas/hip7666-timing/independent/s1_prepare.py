"""Step 1: detrend, build the S58 binary template (two passes), fit pulsations per sector.

Outputs (independent/out/):
  template_s58.npz            spline coefficients of the final template (via pickle)
  prepared_sXX.npz            t, flux (detrended, pulsations subtracted), err, puls model, trend
  prep_log.txt                diagnostics
"""
import os, pickle, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import hiplib as H

OUT = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUT, exist_ok=True)
log = open(os.path.join(OUT, "prep_log.txt"), "a")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()


def build_s58():
    t, f, e = H.load(58)
    ph = H.phase(t)
    ooe = ~H.in_eclipse(t)
    # pass 0: 5 d knot spline, no template
    tr = H.spline_trend(t, f, ooe, 5.0)
    fd = f / tr
    spl, _, _ = H.fold_template(ph, fd)
    T1 = H.Template(spl)
    # detrend with template, 1.5 d knots
    tr = H.spline_trend(t, f / T1.full(ph), ooe, 1.5)
    fd = f / tr
    spl, _, _ = H.fold_template(ph, fd)
    T1 = H.Template(spl)
    say("S58 pass1 template centres:", T1.centers)
    # pulsations from out-of-eclipse residuals
    r = fd - T1.full(ph)
    coef, P1 = H.fit_pulsations(t, r, ooe)
    amps = np.hypot(coef[1::2], coef[2::2])
    say("S58 pass1 pulsation amplitudes (ppt):", np.round(amps * 1e3, 3))
    say("S58 rms OOE before/after pulsations (ppt): %.3f %.3f" %
        (np.std(r[ooe]) * 1e3, np.std((r - P1)[ooe]) * 1e3))
    # pass 2: rebuild template from data minus pulsations
    spl, _, _ = H.fold_template(ph, fd - P1)
    T2 = H.Template(spl)
    # re-detrend with full model, refit pulsations, final template
    tr = H.spline_trend(t, f / (T2.full(ph) + P1), ooe, 1.5)
    fd = f / tr
    coef, P2 = H.fit_pulsations(t, fd - T2.full(ph), ooe)
    spl, xc, yc = H.fold_template(ph, fd - P2)
    T = H.Template(spl)
    r = fd - P2 - T.full(ph)
    say("S58 pass2 template centres (phase):", T.centers,
        " -> offsets from 0 / 0.5 in s: %.2f %.2f" % (T.centers["P"] * H.P_LEE * 86400,
                                                     (T.centers["S"] - 0.5) * H.P_LEE * 86400))
    say("S58 final rms OOE (ppt) %.3f, in-eclipse %.3f" % (np.std(r[ooe]) * 1e3, np.std(r[~ooe]) * 1e3))
    say("S58 template depth: primary min %.4f secondary min %.4f" %
        (T.eclipse("P", 0.0), T.eclipse("S", 0.0)))
    # binned residual structure in eclipse (diagnostic of template adequacy)
    for kind, c in (("P", 0.0), ("S", 0.5)):
        d = H.wrap(ph - c)
        m = np.abs(d) < H.FIT_HALF
        bins = np.linspace(-H.FIT_HALF, H.FIT_HALF, 25)
        idx = np.digitize(d[m], bins)
        br = [np.mean(r[m][idx == i]) * 1e3 for i in range(1, len(bins))]
        say(" binned residual %s (ppt):" % kind, np.round(br, 3))
    np.savez(os.path.join(OUT, "prepared_s58.npz"), t=t, flux=fd - P2, err=e, puls=P2, trend=tr,
             fd=fd)
    np.savez(os.path.join(OUT, "fold_s58.npz"), xc=xc, yc=yc)
    with open(os.path.join(OUT, "template_s58.pkl"), "wb") as fh:
        pickle.dump(spl, fh)
    return T


def prep_other(sector, T, data=None, texp=None, tag=""):
    texp = H.MAIN[sector]["texp"] if texp is None else texp
    t, f, e = H.load(sector) if data is None else data
    ph = H.phase(t)
    ooe = ~H.in_eclipse(t)
    B = T.full_integrated(ph, texp)
    b = 1.0
    P = np.zeros_like(t)
    for it in range(3):
        model = 1.0 + b * (B - 1.0) + P
        tr = H.spline_trend(t, f / model, ooe, 1.5)
        fd = f / tr
        # sector-wide dilution factor b: (fd - P - 1) = b (B - 1)
        y = fd - P - 1.0
        b = np.sum(y * (B - 1.0)) / np.sum((B - 1.0) ** 2)
        r = fd - (1.0 + b * (B - 1.0))
        coef, P = H.fit_pulsations(t, r, ooe, texp=texp)
    r = fd - (1.0 + b * (B - 1.0))
    amps = np.hypot(coef[1::2], coef[2::2])
    X = H.puls_design(t[ooe], texp=texp)
    say("%sS%d: N=%d, texp=%g s, dilution scale b=%.4f, design cond=%.1f" %
        (tag, sector, len(t), texp, b, np.linalg.cond(X)))
    say("S%d pulsation amplitudes (ppt, as observed):" % sector, np.round(amps * 1e3, 3))
    say("S%d rms OOE before/after pulsations (ppt): %.3f %.3f" %
        (sector, np.std(r[ooe]) * 1e3, np.std((r - P)[ooe]) * 1e3))
    np.savez(os.path.join(OUT, "prepared_%ss%02d.npz" % (tag, sector)), t=t, flux=fd - P, err=e, puls=P,
             trend=tr, fd=fd, b=b)


if __name__ == "__main__":
    os.chdir(os.path.join(os.path.dirname(__file__), ".."))
    T = build_s58()
    for s in (18, 85):
        prep_other(s, T)
