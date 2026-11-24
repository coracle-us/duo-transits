"""Step 5 (cross-check): QLP SAP_FLUX for S18/S58/S85, same method as the main products
(S58 SPOC template, per-sector dilution scale, pulsations fitted per sector, same timer,
same error scaling). Compares sector-mean O-C with the main products.
Writes independent/out/eclipse_times_qlp.csv and independent/out/qlp_check.txt
"""
import os, sys, pickle, io, glob
import numpy as np
from astropy.io import fits
sys.path.insert(0, os.path.dirname(__file__))
import hiplib as H
import s1_prepare as S1
from s2_times import time_sector, scale_errors

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
DAY_S = 86400.0
QLP = {18: ("lc/hlsp_qlp_tess_ffi_s0018-0000000308953703_tess_v01_llc.fits", 1800.0),
       58: ("lc/hlsp_qlp_tess_ffi_s0058-0000000308953703_tess_v01_llc.fits", 200.0),
       85: ("lc/hlsp_qlp_tess_ffi_s0085-0000000308953703_tess_v01_llc.fits", 200.0)}


def load_qlp(path):
    d = fits.getdata(path, 1)
    cols = d.columns.names
    ecol = "KSPSAP_FLUX_ERR" if "KSPSAP_FLUX_ERR" in cols else "DET_FLUX_ERR"
    t = np.asarray(d["TIME"], float); f = np.asarray(d["SAP_FLUX"], float)
    e = np.asarray(d[ecol], float); q = np.asarray(d["QUALITY"])
    m = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e > 0) & (q == 0)
    t, f, e = t[m], f[m], e[m]
    med = np.median(f)
    return t, f / med, e / med, ecol


def wmean(x, e):
    w = 1 / e ** 2
    return np.sum(w * x) / w.sum(), 1 / np.sqrt(w.sum())


if __name__ == "__main__":
    os.chdir(os.path.join(HERE, ".."))
    tmpl = H.Template(pickle.load(open(os.path.join(OUT, "template_s58.pkl"), "rb")))
    log = io.StringIO()
    rows = []
    for s, (path, texp) in QLP.items():
        t, f, e, ecol = load_qlp(path)
        print("QLP S%d: N=%d, error column %s" % (s, len(t), ecol))
        S1.prep_other(s, tmpl, data=(t, f, e), texp=texp, tag="qlp_")
        d = np.load(os.path.join(OUT, "prepared_qlp_s%02d.npz" % s))
        r, dropped = time_sector(s, tmpl, log, data=(d["t"], d["flux"], d["err"]), texp=texp)
        print("  kept P=%d S=%d, dropped %d" % (sum(x["type"] == "P" for x in r), sum(x["type"] == "S" for x in r), len(dropped)))
        rows += r
    rows = scale_errors(rows, log)
    with open(os.path.join(OUT, "eclipse_times_qlp.csv"), "w") as fh:
        fh.write("sector,type,cycle,t_btjd,err_d,err_d_unscaled\n")
        for r in sorted(rows, key=lambda x: x["t_btjd"]):
            fh.write("%d,%s,%.1f,%.7f,%.7f,%.7f\n" % (r["sector"], r["type"], r["cycle"], r["t_btjd"], r["err_d"], r["err_raw"]))
    # compare with main products, O-C against the M0 ephemeris from the main products
    import csv
    m0 = np.load(os.path.join(OUT, "m0_fit.npz"))
    c = m0["c"]
    main = list(csv.DictReader(open(os.path.join(HERE, "eclipse_times.csv"))))
    out = open(os.path.join(OUT, "qlp_check.txt"), "w")

    def oc(cyc, t, sec):
        return (t - (c[0] + c[1] * (cyc - 3420.0) + c[2] * sec)) * DAY_S

    for s in (18, 58, 85):
        for k, lab in (("P", "primary"), ("S", "secondary")):
            q = [r for r in rows if r["sector"] == s and r["type"] == k]
            mm = [r for r in main if int(r["sector"]) == s and r["type"] == lab]
            qo = oc(np.array([r["cycle"] for r in q]), np.array([r["t_btjd"] for r in q]), k == "S")
            qe = np.array([r["err_d"] for r in q]) * DAY_S
            mo = oc(np.array([float(r["cycle"]) for r in mm]), np.array([float(r["t_btjd"]) for r in mm]), lab == "secondary")
            me = np.array([float(r["err_d"]) for r in mm]) * DAY_S
            a, ae = wmean(qo, qe); b, be = wmean(mo, me)
            dif = a - b; de = np.hypot(ae, be)
            line = ("S%d %-9s QLP N=%2d O-C=%+7.2f+-%5.2f | main N=%2d O-C=%+7.2f+-%5.2f | diff %+6.2f+-%5.2f (%.1f sigma)%s"
                    % (s, lab, len(q), a, ae, len(mm), b, be, dif, de, abs(dif) / de, "  >3 sigma: systematic" if abs(dif) / de > 3 else ""))
            print(line); out.write(line + "\n")
