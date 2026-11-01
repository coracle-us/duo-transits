#!/usr/bin/env python3
"""Round 2, section 6: QLP mid-times from Sector 94 on.

None of the shortlist's QLP S94+ sectors has a TESS-SPOC or SPOC light curve,
so the comparison is against our own TESScut photometry. Two measurements:

1. Clock check. A TESScut cutout gives spacecraft time per FFI cadence
   (TIME - TIMECORR). The target's own barycentric (Roemer) correction comes
   from the JPL Horizons TESS ephemeris (this matches TESS-SPOC's per-target
   TIMECORR to 0.001 s). Matching QLP cadences by CADENCENO gives QLP's time
   offset per cadence. TESScut's own TIMECORR is for another direction (the
   CCD centre, presumably), so its offset is reported and removed too.
2. Transit refit. A trapezoid plus linear baseline fitted to the QLP light
   curve (QLP times) and to TESScut aperture photometry (corrected times).

  python round2/s6_qlp_times.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import st2_pipeline as pl  # noqa: E402
import st2core as core  # noqa: E402
import tesscut_phot as tp  # noqa: E402

WORK, OUT = HERE / "work", HERE / "results" / "round2"
C_AU_S = 499.004784
# (tic, sector): shortlist QLP S94+ events from the brief, plus one pre-S94 control
TARGETS = [(152464302, 99), (107113048, 101), (438588300, 101), (241550976, 100), (241550976, 101),
           (241550976, 102), (275647095, 100), (275647095, 101), (302377820, 95), (241550976, 65)]


def roemer_s(jd_tdb, ra, dec):
    """Light travel time (s) from TESS to the solar-system barycentre along the target direction."""
    from astroquery.jplhorizons import Horizons
    grid = np.arange(np.nanmin(jd_tdb) - 0.5, np.nanmax(jd_tdb) + 1.0, 0.5)
    vec = []
    for chunk in np.array_split(grid, int(np.ceil(grid.size / 50))):
        o = Horizons(id="-95", location="@0", epochs=[round(float(x), 6) for x in chunk]).vectors(refplane="earth")
        vec.append(np.column_stack([np.asarray(o[k], float) for k in "xyz"]))
    r = np.vstack(vec)
    n = np.array([np.cos(np.radians(dec)) * np.cos(np.radians(ra)),
                  np.cos(np.radians(dec)) * np.sin(np.radians(ra)), np.sin(np.radians(dec))])
    return np.interp(jd_tdb, grid, r @ n * C_AU_S)


def trapezoid(p, t):
    tc, depth, t14, t12, a, b = p
    t12 = np.clip(t12, 1e-4, t14 / 2)
    x = np.abs(t - tc)
    m = np.ones_like(t)
    full = x <= t14 / 2 - t12
    ing = (x > t14 / 2 - t12) & (x < t14 / 2)
    m[full] -= depth
    m[ing] -= depth * (t14 / 2 - x[ing]) / t12
    return m * (a + b * (t - tc))


def fit_transit(t, f, tc0, dur0, depth0):
    w = np.abs(t - tc0) < max(2.0 * dur0, 0.5)
    t, f = t[w], f[w]
    if t.size < 20:
        return dict(tc=np.nan, tc_err=np.nan, depth=np.nan, t14=np.nan, beta=np.nan, n=int(t.size))
    p0 = [tc0, depth0, dur0, dur0 / 8, 1.0, 0.0]
    lo = [tc0 - dur0 / 2, 0, dur0 / 3, 1e-3, 0.9, -1]
    hi = [tc0 + dur0 / 2, 5 * depth0, 3 * dur0, dur0, 1.1, 1]
    r = least_squares(lambda p: trapezoid(p, t) - f, p0, bounds=(lo, hi))
    res = r.fun
    J = r.jac
    s2 = np.sum(res ** 2) / max(res.size - len(p0), 1)
    cov = np.linalg.pinv(J.T @ J) * s2
    # inflate for correlated noise: rms of residuals binned at ~ingress length vs white expectation
    nb = max(int(round(max(r.x[3], 15 / 1440) / np.median(np.diff(t)))), 1)
    beta = 1.0
    if nb > 1 and res.size > 4 * nb:
        k = res.size // nb
        binned = res[: k * nb].reshape(k, nb).mean(1)
        beta = max(1.0, np.std(binned) / (np.std(res) / np.sqrt(nb)))
    return dict(tc=r.x[0], tc_err=np.sqrt(cov[0, 0]) * beta, depth=r.x[1], t14=r.x[2], beta=beta, n=t.size)


def qlp_lightcurve(tic, sector):
    """Every row of the QLP file in order (time, normalised flux, quality == 0 and finite)."""
    import lightkurve as lk
    sr = lk.search_lightcurve(f"TIC {tic}", sector=sector, author="QLP")
    lc = sr[0].download(download_dir=str(WORK / "lc_cache"), quality_bitmask="none")
    t = np.asarray(lc.time.value, float)
    f = np.asarray(lc.flux.value, float)
    q = np.asarray(getattr(lc.quality, "value", lc.quality))
    ok = (q == 0) & np.isfinite(t) & np.isfinite(f)
    return t, f / np.nanmedian(f[ok]), ok


def qlp_minus_spacecraft(tq_all, okq, t_sc, good):
    """QLP time minus spacecraft time (s) for each good QLP row, and how its branch was fixed.

    CADENCENO is no help (QLP numbers its own cadences, TESScut's are 0). Each
    QLP time is compared with the nearest FFI exposure and the difference is
    unwrapped across the sector, which gives it exactly up to a whole FFI
    cadence. When both files have the same number of rows (one per exposure),
    row i is the same exposure and fixes that whole number. Otherwise the
    branch with |median| under half a cadence is taken and flagged: the
    absolute value is then known only modulo one cadence.
    """
    ts = np.sort(t_sc[good])
    cad = np.median(np.diff(ts)) * 86400
    tq = tq_all[okq]
    j = np.clip(np.searchsorted(ts, tq), 1, ts.size - 1)
    j = np.where(np.abs(ts[j - 1] - tq) < np.abs(ts[j] - tq), j - 1, j)
    r = np.unwrap((tq - ts[j]) * 86400, period=cad)
    if tq_all.size == t_sc.size:
        exact = (tq_all - t_sc)[okq] * 86400
        ok = np.isfinite(exact)
        k = np.round(np.median(exact[ok] - r[ok]) / cad)
        how = "rows aligned"
    else:
        k = -np.round(np.median(r) / cad)
        how = "row counts differ: known modulo one cadence"
    return r + k * cad, ts[j], how, cad


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cov = pd.read_csv(WORK / "coverage.csv").set_index("tic")
    ev = pd.read_csv(WORK / "events_all.csv")
    ev = ev[ev.passes == True]
    rows = []
    for tic, sector in TARGETS:
        try:
            rows += one_target(tic, sector, cov, ev)
        except Exception as e:
            print(f"TIC {tic} S{sector}: FAILED {e}")
            rows.append(dict(tic=tic, sector=sector, error=str(e)[:200]))
        pd.DataFrame(rows).to_csv(OUT / "qlp_s94_timing.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


def one_target(tic, sector, cov, ev):
    rows = []
    if True:
        c = cov.loc[tic]
        evs = ev[(ev.tic == tic) & (ev.sector == sector)].sort_values("snr", ascending=False)
        evs = evs.loc[~evs.t2.round(1).duplicated()]
        if evs.empty:
            print(f"TIC {tic} S{sector}: no passing event")
            return rows
        tpf = tp.cutout(c.ra, c.dec, sector, WORK / "tesscut")
        d = tpf.hdu[1].data
        good = np.isfinite(d["TIME"]) & np.isfinite(d["TIMECORR"]) & (d["TIME"] > 0)
        t_sc = d["TIME"] - d["TIMECORR"]  # spacecraft time (TJD), days
        rm = np.full(t_sc.size, np.nan)
        rm[good] = roemer_s(t_sc[good] + 2457000, c.ra, c.dec)
        t_true = t_sc + rm / 86400
        cut_off_s = (d["TIME"] - t_true) * 86400
        tq_all, fq_all, okq = qlp_lightcurve(tic, sector)
        qlp_minus_sc, tsc_near, how, cad = qlp_minus_spacecraft(tq_all, okq, t_sc, good)
        tq, fq = tq_all[okq], fq_all[okq]
        o = np.argsort(t_sc[good])
        rm_q = np.interp(tsc_near, t_sc[good][o], rm[good][o])
        q_off = qlp_minus_sc - rm_q  # QLP minus target-corrected time, s
        m = np.isfinite(q_off)
        tracks = np.ptp(q_off[m]) < 0.3 * np.ptp(rm[good])
        if how != "rows aligned":
            if tracks:  # QLP applied the target's correction: offset is constant, take the branch nearest 0
                kk = -np.round(np.median(q_off[m]) / cad)
                q_off, qlp_minus_sc = q_off + kk * cad, qlp_minus_sc + kk * cad
                how += "; QLP follows the target's correction, branch nearest 0 taken"
            else:
                how += "; QLP correction constant over the sector, absolute offset unknown modulo one cadence"
        q_off_s = q_off[m]
        # TESScut photometry on the corrected clock
        t_c, f_c, e_c, info = tp.sector_lightcurve(tic, c.ra, c.dec, sector, WORK / "tesscut")
        t_c = np.interp(t_c, d["TIME"][good], t_true[good])
        for h in evs.itertuples():
            near_q = m & (np.abs(tq - h.t2) < 1.0)
            near_c = good & (np.abs(d["TIME"] - h.t2) < 1.0)
            fq_ = fit_transit(tq, fq, h.t2, h.dur_ref_h / 24, h.depth_ref_ppm * 1e-6)
            fc_ = fit_transit(t_c, f_c, h.t2, h.dur_ref_h / 24, h.depth_ref_ppm * 1e-6)
            qo = q_off[near_q]
            rows.append(dict(
                tic=tic, sector=sector, cand_id=h.cand_id, t2_box=h.t2, snr=h.snr,
                n_qlp_rows=int(tq_all.size), n_ffi_rows=int(t_sc.size), branch=how,
                qlp_cadence_s=float(np.nanmedian(np.diff(tq_all)) * 86400),
                ffi_cadence_spacecraft_s=float(cad),
                qlp_minus_spacecraft_median_s=float(np.nanmedian(qlp_minus_sc)),
                qlp_minus_spacecraft_spread_s=float(np.nanpercentile(qlp_minus_sc, 99) - np.nanpercentile(qlp_minus_sc, 1)),
                target_roemer_range_s=float(np.nanmax(rm) - np.nanmin(rm)),
                qlp_clock_offset_at_event_s=float(np.median(qo)) if qo.size else np.nan,
                qlp_clock_offset_sector_min_s=float(np.min(q_off_s)) if q_off_s.size else np.nan,
                qlp_clock_offset_sector_max_s=float(np.max(q_off_s)) if q_off_s.size else np.nan,
                qlp_follows_target_correction=bool(tracks),
                tesscut_clock_offset_at_event_s=float(np.median(cut_off_s[near_c])),
                tc_qlp=fq_["tc"], tc_qlp_err_min=fq_["tc_err"] * 1440,
                tc_tesscut=fc_["tc"], tc_tesscut_err_min=fc_["tc_err"] * 1440,
                tc_qlp_clock_corrected=(fq_["tc"] - np.median(qo) / 86400
                                        if qo.size and ("aligned" in how or tracks) else np.nan),
                qlp_minus_tesscut_min=(fq_["tc"] - fc_["tc"]) * 1440,
                qlp_minus_tesscut_err_min=np.hypot(fq_["tc_err"], fc_["tc_err"]) * 1440,
                depth_qlp_ppm=fq_["depth"] * 1e6, depth_tesscut_ppm=fc_["depth"] * 1e6,
                t14_qlp_h=fq_["t14"] * 24, t14_tesscut_h=fc_["t14"] * 24, tesscut_n_ap=info["n_ap"]))
            r = rows[-1]
            print(f"TIC {tic} S{sector} T~{h.t2:.3f}: QLP clock {r['qlp_clock_offset_at_event_s']:+.1f} s "
                  f"(sector {r['qlp_clock_offset_sector_min_s']:+.1f} to {r['qlp_clock_offset_sector_max_s']:+.1f}), "
                  f"TESScut clock {r['tesscut_clock_offset_at_event_s']:+.1f} s; fit QLP-TESScut "
                  f"{r['qlp_minus_tesscut_min']:+.2f} +/- {r['qlp_minus_tesscut_err_min']:.2f} min "
                  f"(QLP - spacecraft {r['qlp_minus_spacecraft_median_s']:+.1f} s, spread {r['qlp_minus_spacecraft_spread_s']:.1f} s)")
    return rows


if __name__ == "__main__":
    main()
