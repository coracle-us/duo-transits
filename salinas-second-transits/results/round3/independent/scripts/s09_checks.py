"""Vetting checks per star, written to out/checks_<tic>.txt:
 - odd/even depth (depth scale from the per-transit fits, timing-source rows)
 - depth consistency (chi2 of depth scales about their weighted mean, per source)
 - secondary eclipse at phase 0.5 (primary shape, depth scale fitted at T0+(k+0.5)P on raw flux)
 - stellar density implied by the fitted a/R* (circular orbit) vs TIC
 - QLP minus TESScut mid-times in S>=94 (the clock-offset caveat)
"""
import sys, os, csv, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import TARGETS, OUT, load_lc_files, load_tesscut_phot
from s07_aliases import depth_at

TIC_RHO = {275647095: 0.6338, 241550976: 1.28278, 107113048: 0.913414, 230133616: 0.0725086}  # solar units


def wmean(x, e):
    w = 1 / e**2
    m = np.sum(w * x) / np.sum(w)
    return m, 1 / np.sqrt(np.sum(w)), np.sum(w * (x - m) ** 2)


def run(tic):
    L = []
    pr = L.append
    eph = json.load(open(os.path.join(OUT, f"ephem_{tic}.json")))
    sh = json.load(open(os.path.join(OUT, f"shape_{tic}.json")))
    P, T0 = eph["P"], eph["T0_E0"]
    rows = list(csv.DictReader(open(os.path.join(OUT, f"times_{tic}.csv"))))
    tr = [r for r in rows if r["tc"] and r["timing_source"] == "1"]
    E = np.array([int(r["epoch"]) for r in tr])
    s = np.array([float(r["depth_scale"]) for r in tr])
    se = np.array([float(r["depth_scale_err"]) for r in tr])
    src = np.array([r["source"] for r in tr])
    pr(f"TIC {tic}  P={P:.6f}  T0={T0:.5f}  model depth {sh['depth']*1e6:.0f} ppm")
    if (E % 2 == 0).any() and (E % 2 == 1).any():
        mo, eo, _ = wmean(s[E % 2 == 1], se[E % 2 == 1])
        me, ee, _ = wmean(s[E % 2 == 0], se[E % 2 == 0])
        pr(f"odd/even depth scale: odd {mo:.3f}±{eo:.3f} (n={np.sum(E%2==1)})  even {me:.3f}±{ee:.3f} "
           f"(n={np.sum(E%2==0)})  diff = {(mo-me)/np.hypot(eo,ee):+.2f} sigma")
    else:
        pr("odd/even: only one parity observed")
    for so in sorted(set(src)):
        m = src == so
        if m.sum() >= 2:
            mm, ee, c2 = wmean(s[m], se[m])
            pr(f"depth scale [{so}] n={m.sum()}: mean {mm:.3f}±{ee:.3f}, chi2={c2:.1f}/{m.sum()-1}, "
               f"min {s[m].min():.2f} max {s[m].max():.2f}")
        else:
            pr(f"depth scale [{so}] n=1: {s[m][0]:.3f}±{se[m][0]:.3f}")
    # secondary eclipse at phase 0.5
    lcs = [lc for lc in load_lc_files(tic) + load_tesscut_phot(tic) if len(lc["t"]) > 100]
    tmin = min(lc["t"].min() for lc in lcs); tmax = max(lc["t"].max() for lc in lcs)
    sec = []
    for k in np.arange(np.ceil((tmin - T0) / P - 0.5), np.floor((tmax - T0) / P - 0.5) + 1):
        m = depth_at(lcs, sh, T0 + (k + 0.5) * P, TARGETS[tic]["dur"])
        if m is not None and m[2] >= 0.5:
            sec.append(m[:2])
    if sec:
        a = np.array(sec)
        mm, ee, c2 = wmean(a[:, 0], a[:, 1])
        pr(f"phase-0.5 eclipse (circular): {len(sec)} epochs, depth = {mm*sh['depth']*1e6:.0f} ± "
           f"{ee*sh['depth']*1e6:.0f} ppm ({mm:.3f}±{ee:.3f} of primary)")
    # density
    rho = 3 * np.pi * (sh["a"]) ** 3 / (6.674e-8 * (P * 86400) ** 2) / 1.41
    pr(f"shape fit: Rp/R*={sh['rp']:.4f} a/R*={sh['a']:.2f} b={sh['b']:.2f} -> rho*(circ) = {rho:.2f} rho_sun; "
       f"TIC rho = {TIC_RHO[tic]:.2f} rho_sun")
    # QLP - TESScut
    by = {}
    for r in rows:
        if r["tc"] and int(r["sector"]) >= 94:
            by.setdefault((r["sector"], r["epoch"]), {})[r["source"]] = (float(r["tc"]), float(r["tc_err"]), float(r["cov_in"]))
    for (sct, ep), d in sorted(by.items()):
        if "QLP" in d and "TESScut" in d:
            q, c = d["QLP"], d["TESScut"]
            pr(f"S{sct} E{ep}: QLP - TESScut = {(q[0]-c[0])*1440:+.2f} ± {np.hypot(q[1], c[1])*1440:.2f} min "
               f"(cov QLP {q[2]:.2f}, TESScut {c[2]:.2f})")
    txt = "\n".join(L)
    open(os.path.join(OUT, f"checks_{tic}.txt"), "w").write(txt + "\n")
    print(txt + "\n")


if __name__ == "__main__":
    for tic in ([int(a) for a in sys.argv[1:]] or TARGETS):
        run(tic)
