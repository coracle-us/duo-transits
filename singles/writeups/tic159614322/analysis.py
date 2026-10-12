"""
Numbers for the TIC 159614322 (CTOI 159614322.01) write-up.

  1. Gaia DR3 astrometry + astrophysical parameters (by position; the stored source_id is
     float-corrupted, see gaia_flags.py)
  2. Alias map: for n = 1..12, every predicted epoch inside TESS data, with the verdict of the
     window test (excluded = flat, unconstrained, no-data) from the same grid the pipeline used
  3. Transit shape: joint fit of both transits for a grid of impact parameters b (rp, a/R*, T0s
     free; P fixed at each alias); chi2 profile in b and the implied circular-orbit stellar
     density rho_circ, compared with the TIC/Gaia density
  4. Predicted transits of each surviving alias, 2026-10 .. 2029-01

Output: analysis.json next to this script.
Run:    python singles/writeups/tic159614322/analysis.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
import alias_pipeline as ap  # noqa: E402
import pilot  # noqa: E402
import vet_candidate as vc  # noqa: E402

TIC, RA, DEC = 159614322, 220.00977, 65.484041
DEPTH, DUR_H = 3130e-6, 8.773
G_CGS = 6.674e-8


def gaia():
    from astroquery.gaia import Gaia
    q = f"""SELECT g.source_id, g.phot_g_mean_mag, g.bp_rp, g.parallax, g.parallax_error, g.ruwe,
            g.ipd_frac_multi_peak, g.ipd_gof_harmonic_amplitude, g.non_single_star, g.radial_velocity,
            g.radial_velocity_error, g.rv_chisq_pvalue, g.rv_renormalised_gof, g.rv_nb_transits,
            a.teff_gspphot, a.logg_gspphot, a.mh_gspphot, a.radius_gspphot, a.radius_flame,
            a.mass_flame, a.age_flame, a.lum_flame, a.distance_gspphot,
            DISTANCE(POINT({RA},{DEC}),POINT(g.ra,g.dec))*3600 AS sep
            FROM gaiadr3.gaia_source AS g
            LEFT JOIN gaiadr3.astrophysical_parameters AS a ON a.source_id = g.source_id
            WHERE 1=CONTAINS(POINT(g.ra,g.dec),CIRCLE({RA},{DEC},0.0015)) ORDER BY g.phot_g_mean_mag"""
    r = Gaia.launch_job(q).get_results().to_pandas().iloc[0].to_dict()
    r["source_id"] = str(int(Gaia.launch_job(q).get_results()["source_id"][0]))
    return {k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in r.items()}


def alias_map(lc, epoch, D):
    g = pilot.grid_windows(lc, D, DEPTH, epoch)
    Tg = g.T_pred.to_numpy()
    rows = []
    dT = 2423.5486 - 1703.8276
    for n in range(1, 13):
        P = dT / n
        ks = np.arange(np.ceil((Tg.min() - epoch) / P), np.floor((Tg.max() - epoch) / P) + 1)
        for k in ks:
            T = epoch + k * P
            i = np.argmin(np.abs(Tg - T))
            if abs(Tg[i] - T) > D / 4:
                continue  # no data at all near this epoch
            rows.append(dict(n=n, P=P, k=int(k), T=round(float(T), 3), sector=int(g.sector.iloc[i]),
                             product=g["product"].iloc[i], verdict=g.verdict.iloc[i],
                             depth_ppm=round(float(g.depth.iloc[i]) * 1e6) if "depth" in g and np.isfinite(g.depth.iloc[i]) else None,
                             err_ppm=round(float(g.depth_err.iloc[i]) * 1e6) if "depth_err" in g and np.isfinite(g.depth_err.iloc[i]) else None))
    return rows


def b_profile(P):
    """chi2 of the two-transit joint fit at fixed b, P; free rp, aR, T0 per epoch."""
    prods = vc.all_products(TIC)
    D = DUR_H / 24
    data = []
    for Tm, sec, pref in ((1703.8276, 14, "TESS-SPOC"), (2423.5486, 41, "SPOC_120")):
        p = next(q for q in prods if q["sector"] == sec and q["product"] == pref)
        w = vc.window(p, Tm, D)
        data.append((w[0], w[1], w[2], p["exptime"], Tm))
    out = []
    for b in np.round(np.arange(0.0, 1.01, 0.05), 2):
        def resid(x):
            rp, aR, t1, t2 = x
            r = []
            for (t, f, s, ex, Tm), T0 in zip(data, (t1, t2)):
                r.append((f - vc.model(t, T0, P, rp, b, aR, ex)) / s)
            return np.concatenate(r)
        aR0 = P * 24 / (np.pi * DUR_H) * np.sqrt(max(1.05 ** 2 - b ** 2, 0.02))
        best = None
        for scale in (0.7, 1.0, 1.4):
            try:
                s = least_squares(resid, [np.sqrt(DEPTH), aR0 * scale, 1703.8276, 2423.5486],
                                  bounds=([0.01, 1.5, 1703.6, 2423.3], [0.3, 5000, 1704.1, 2423.8]))
            except ValueError:
                continue
            if best is None or s.cost < best.cost:
                best = s
        rp, aR, t1, t2 = best.x
        rho = 3 * np.pi / (G_CGS * (P * 86400) ** 2) * aR ** 3
        out.append(dict(b=float(b), chi2=float(2 * best.cost), rp=float(rp), aR=float(aR), rho_circ=float(rho),
                        T1=float(t1), T2=float(t2), n=int(sum(len(d[0]) for d in data))))
    return out


def predictions(periods, T2, sigT2_d, dT_err):
    rows = []
    t0, t1 = Time("2026-10-11").jd - 2457000, Time("2029-01-01").jd - 2457000
    for n, P in periods:
        sP = dT_err / n
        for k in range(int(np.ceil((t0 - T2) / P)), int(np.floor((t1 - T2) / P)) + 1):
            T = T2 + k * P
            sig = np.hypot(sigT2_d, k * sP) * 24 * 60
            rows.append(dict(n=n, P=round(P, 4), cycle_after_S41=k, T_btjd=round(T, 4),
                             mid_UTC=Time(T + 2457000, format="jd", scale="tdb").utc.iso[:16],
                             sigma_min=round(float(sig))))
    return rows


def main():
    lc = ap.load_lightcurves(TIC)
    D = DUR_H / 24
    res = {}
    res["sectors"] = sorted({int(s) for s in lc.sector.unique()})
    res["gaia"] = gaia()
    res["alias_map"] = alias_map(lc, 2423.5486, D)
    # shape: profile in b at the longest and shortest surviving aliases
    res["b_profile"] = {str(n): b_profile(719.6572 / n) for n in (1, 8)}
    # T0 errors from the vet fit are NaN; take the timing error from the chi2 curvature
    # in the joint fit at the best b (approx. 1-sigma from delta chi2 = 1 along T0 alone)
    bp = res["b_profile"]["8"]
    bb = min(bp, key=lambda r: r["chi2"])
    res["best_b_P90"] = bb
    T1, T2 = bb["T1"], bb["T2"]
    res["dT"] = T2 - T1
    periods = [(n, (T2 - T1) / n) for n in (1, 2, 4, 5, 8)]
    res["periods"] = periods
    sigT = 5.0 / 60 / 24   # conservative 5-min per-epoch timing error (2 transits, SNR ~20-28)
    res["timing_sigma_assumed_min"] = 5.0
    res["predictions"] = predictions(periods, T2, sigT, np.sqrt(2) * sigT)
    (HERE / "analysis.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k not in ("alias_map", "predictions", "b_profile")},
                     indent=1, default=float))
    am = pd.DataFrame(res["alias_map"])
    print(am.to_string(index=False))
    for n, prof in res["b_profile"].items():
        print("n", n, [(r["b"], round(r["chi2"] - min(x["chi2"] for x in prof), 1), round(r["rho_circ"], 2)) for r in prof])
    print(pd.DataFrame(res["predictions"]).to_string(index=False))


if __name__ == "__main__":
    main()
