"""
Mid-times of the three TIC 159614322 transits and the two linear ephemerides that fit them.

Each transit: transit shape fixed (rp, b, a/R* from the S14+S41 joint fit at P = 89.97 d, b = 0;
see analysis.json), free mid-time, depth scale and a linear baseline; error from the
Jacobian, scaled by sqrt(reduced chi2) if > 1.
  S14  TESS-SPOC PDC (30-min FFI)
  S21  TESS-SPOC SAP (30-min FFI); egress falls in the orbit-1 downlink gap, so the window
       runs from 1883.00 (after the 1882.97 momentum dump) to the end of data at 1883.97
  S41  SPOC 2-min PDC
S21 cross-check: QLP and our own TESScut photometry (s21_check.py, own_s21 below).

Ephemerides (revised after second-reviewer item 5, 2026-10-11):
  ADOPTED  P and T0 from the two FULL transits (S14, S41) only: P_n = dT/n, n = 4 and 8.
           S21 (partial, no egress) is used only to select n = 4/8 and as a consistency check:
           its O-C against the two-transit ephemeris is reported per S21 reduction.
  three-transit fit (cycles 0,1,4 / 0,2,8) kept for comparison; its chi2 (1 dof) sets an error
  inflation factor sqrt(chi2/dof) used for the PLANNING sigma of predicted transits
  (sigma_plan = sigma_formal * max(1, sqrt(chi2))).
Weighted least squares; chi2 and O-C per transit.
Predictions: every transit of each ephemeris 2026-10-11 .. 2029-01-01, with sigma, plus
night-time visibility (Sun < -12 deg, target altitude > 20 deg for the whole transit
including 1 h of baseline on each side, or at least > 50% of it) at a few observatory sites.

Output: ephemeris.json, predictions.csv
Run:    python singles/writeups/tic159614322/ephemeris.py
"""
import json
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import vet_candidate as vc  # noqa: E402

TIC, RA, DEC = 159614322, 220.00977, 65.484041
D = 8.773 / 24
SITES = {
    "Boston, MA": (42.36, -71.06, 20),
    "Chicago, IL": (41.88, -87.63, 180),
    "Kitt Peak, AZ": (31.96, -111.60, 2096),
    "La Palma": (28.76, -17.88, 2326),
    "Mt Wilson, CA": (34.22, -118.06, 1742),
}


def shape():
    b = json.load(open(HERE / "analysis.json"))["best_b_P90"]
    return b["rp"], b["b"], b["aR"]


def fit_time(p, Tg, col, lo, hi, rp, b, aR):
    t = p["time"]
    m = (t > lo) & (t < hi) & np.isfinite(p[col])
    t, f = t[m], p[col][m]
    f = f / np.median(f)
    P = 89.9651

    def model(x):
        T0, s, c0, c1 = x
        return (c0 + c1 * (t - Tg)) * (1 - s * (1 - vc.model(t, T0, P, rp, b, aR, p["exptime"])))

    oot = np.abs(t - Tg) > D / 2 + 0.05
    sig = np.std(f[oot] - np.polyval(np.polyfit(t[oot], f[oot], 1), t[oot]))
    best = None
    for dT in np.arange(-0.1, 0.101, 0.01):
        s_ = least_squares(lambda x: (f - model(x)) / sig, [Tg + dT, 1, 1, 0],
                           bounds=([Tg - 0.3, 0, 0.9, -1], [Tg + 0.3, 3, 1.1, 1]))
        if best is None or s_.cost < best.cost:
            best = s_
    J = best.jac
    red = 2 * best.cost / (len(t) - 4)
    cov = np.linalg.inv(J.T @ J) * max(red, 1)
    return dict(product=p["product"], column=col, sector=int(p["sector"]), T=float(best.x[0]),
                T_err=float(np.sqrt(cov[0, 0])), depth_scale=float(best.x[1]),
                depth_scale_err=float(np.sqrt(cov[1, 1])), red_chi2=float(red), n=int(len(t)))


def ephem(times, errs, cycles):
    A = np.column_stack([np.ones(len(cycles)), cycles])
    W = np.diag(1 / np.asarray(errs) ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    T0, P = cov @ A.T @ W @ np.asarray(times)
    oc = np.asarray(times) - (T0 + P * np.asarray(cycles))
    return dict(T0=float(T0), T0_err=float(np.sqrt(cov[0, 0])), P=float(P), P_err=float(np.sqrt(cov[1, 1])),
                cov=cov.tolist(), chi2=float(np.sum((oc / np.asarray(errs)) ** 2)),
                O_minus_C_min=[float(v * 1440) for v in oc])


def visibility(Tmid_btjd):
    tgt = SkyCoord(RA * u.deg, DEC * u.deg)
    out = {}
    tt = Time(Tmid_btjd + 2457000 + np.linspace(-(D / 2 + 1 / 24), D / 2 + 1 / 24, 25), format="jd", scale="tdb")
    sun_all = get_sun(tt)
    for name, (lat, lon, h) in SITES.items():
        loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=h * u.m)
        fr = AltAz(obstime=tt, location=loc)
        alt = tgt.transform_to(fr).alt.deg
        sun = sun_all.transform_to(fr).alt.deg
        ok = (alt > 20) & (sun < -12)
        out[name] = round(float(ok.mean()), 2)
    return out


def main():
    rp, b, aR = shape()
    prods = vc.all_products(TIC)
    get = lambda s, pr: next(q for q in prods if q["sector"] == s and q["product"] == pr)  # noqa: E731
    tr = [fit_time(get(14, "TESS-SPOC"), 1703.828, "flux", 1703.828 - 1.2, 1703.828 + 1.2, rp, b, aR),
          fit_time(get(21, "TESS-SPOC"), 1883.77, "sap_flux", 1883.00, 1883.975, rp, b, aR),
          fit_time(get(41, "SPOC_120"), 2423.549, "flux", 2423.549 - 1.2, 2423.549 + 1.2, rp, b, aR)]
    cross = [fit_time(get(21, "QLP"), 1883.77, "flux", 1883.00, 1883.975, rp, b, aR),
             fit_time(get(21, "TESS-SPOC"), 1883.77, "flux", 1883.00, 1883.975, rp, b, aR),
             fit_time(get(14, "QLP"), 1703.828, "flux", 1703.828 - 1.2, 1703.828 + 1.2, rp, b, aR),
             fit_time(get(41, "TESS-SPOC"), 2423.549, "flux", 2423.549 - 1.2, 2423.549 + 1.2, rp, b, aR)]
    T = [r["T"] for r in tr]
    E = [r["T_err"] for r in tr]
    res = dict(transits=tr, cross_checks=cross, dT=T[2] - T[0])
    res["ephemerides_3transit"] = {"n4": ephem(T, E, [0, 1, 4]), "n8": ephem(T, E, [0, 2, 8])}
    # adopted: two full transits only
    res["ephemerides"] = {"n4": ephem([T[0], T[2]], [E[0], E[2]], [0, 4]),
                          "n8": ephem([T[0], T[2]], [E[0], E[2]], [0, 8])}
    infl = float(max(1.0, np.sqrt(res["ephemerides_3transit"]["n8"]["chi2"])))
    res["sigma_inflation"] = infl
    e8 = res["ephemerides"]["n8"]
    c8 = np.array(e8["cov"])
    Tp = e8["T0"] + 2 * e8["P"]
    sp = float(np.sqrt(c8[0, 0] + 4 * c8[1, 1] + 4 * c8[0, 1]))
    res["S21_check"] = dict(predicted=Tp, predicted_err=sp, reductions=[
        dict(product=r["product"], column=r["column"], T=r["T"], T_err=r["T_err"],
             O_minus_C_min=(r["T"] - Tp) * 1440, sigma=(r["T"] - Tp) / np.hypot(r["T_err"], sp))
        for r in [tr[1]] + [c for c in cross if c["sector"] == 21]])
    # the other aliases, for completeness: S21 time vs their prediction (none has a transit there)
    res["S21_offset_from_dT_over_n_days"] = {n: float(((T[1] - T[0]) / ((T[2] - T[0]) / n)) % 1) for n in (1, 2, 4, 5, 8)}
    rows = []
    t0, t1 = Time("2026-10-11").jd - 2457000, Time("2029-01-01").jd - 2457000
    for key, n in (("n4", 4), ("n8", 8)):
        e = res["ephemerides"][key]
        c = np.array(e["cov"])
        for k in range(int(np.ceil((t0 - e["T0"]) / e["P"])), int(np.floor((t1 - e["T0"]) / e["P"])) + 1):
            Tk = e["T0"] + k * e["P"]
            sig = np.sqrt(c[0, 0] + k ** 2 * c[1, 1] + 2 * k * c[0, 1]) * 1440
            row = dict(alias=f"P = {e['P']:.3f} d (n={n})", cycle=k, T_btjd=round(Tk, 4),
                       mid_UTC=Time(Tk + 2457000, format="jd", scale="tdb").utc.iso[:16], sigma_min=round(float(sig)),
                       sigma_plan_min=round(float(sig * infl)))
            row.update({f"vis_{s}": v for s, v in visibility(Tk).items()})
            rows.append(row)
    pr = pd.DataFrame(rows)
    pr.to_csv(HERE / "predictions.csv", index=False)
    (HERE / "ephemeris.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items()}, indent=1)[:4000])
    pd.set_option("display.width", 250)
    print(pr.to_string(index=False))


if __name__ == "__main__":
    main()
