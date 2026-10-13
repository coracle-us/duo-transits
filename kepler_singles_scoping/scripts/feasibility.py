"""
Feasibility of recovering periods of Kepler single/duo transits with TESS.

For every candidate in data/kepler_singles_targets.csv:
  * status: KOI disposition (cumulative table), confirmed planets on the host (pscomppars), TOI / CTOI on the TIC
    (repo tables in data/raw/, read only), KM19 'flagged' table, depth > 2.5 %, U16 'EB-like' -> class
  * expected per-transit TESS SNR = depth / sigma_1h(Tmag) * sqrt(dur_h) with the empirical noise model
    from data/noise_model.txt (falls back to 220 ppm * 10^(0.2 (T-9.1)) if absent); duration capped at 30 h
  * Kepler coverage = quarters the star was observed (MAST) x approximate quarter windows (BKJD)
  * TESS coverage   = tess-point sectors x orbit windows (singles/data/tess_orbit_times_latest.csv)
  * singles: Monte Carlo over P with prior dN/dlnP ∝ P^(-2/3) (flat occurrence x transit probability)
    from 50 d to 7300 d (20 yr); reject P that puts another transit inside Kepler data; p_repeat = fraction of the
    surviving prior mass with >= 1 transit centre inside TESS data observed by 2026-10-12 (and through the last
    sector tess-point knows about). Same with the catalogue's own period posterior (lognormal approx) where given.
    For simulated repeats, also count the aliases P = dT/n left after the Kepler + TESS flat-data cuts.
  * duos: aliases P_n = dT/n (P_n >= 10 d) with no extra transit in Kepler data; weight ∝ P^(-2/3);
    p_repeat = weighted fraction of aliases with a transit in TESS data. 1 surviving alias = period already known.

Outputs: data/feasibility.csv, data/feasibility_summary.txt
"""
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent
RNG = np.random.default_rng(42)
TODAY_BTJD = (pd.Timestamp(date.today()).to_julian_date()) - 2457000.0
TESS_EFF = 0.85   # fraction of in-sector time with usable data (quality flags, scattered light)
DUR_CAP_H = 30.0  # longest duration that a 13-d TESS orbit + detrending can use at full weight
SNR_MIN = 7.0

# Kepler quarter windows in BKJD (approximate, from the Kepler data release notes; monthly downlinks ignored)
KQ = {0: (120.53, 130.30), 1: (131.51, 164.98), 2: (169.76, 258.47), 3: (260.23, 349.49), 4: (352.38, 442.20),
      5: (443.49, 538.17), 6: (539.46, 629.30), 7: (630.18, 719.55), 8: (735.38, 802.34), 9: (808.52, 905.92),
      10: (906.85, 1001.20), 11: (1002.25, 1098.32), 12: (1099.42, 1181.97), 13: (1182.76, 1273.06),
      14: (1274.15, 1371.32), 15: (1373.48, 1471.13), 16: (1472.12, 1558.40), 17: (1559.21, 1591.00)}


def noise_fn():
    f = HERE / "data" / "noise_model.txt"
    if f.exists():
        line = f.read_text().splitlines()[0]
        rhs = line.split("] =")[1]
        a = float(rhs.split("+")[0])
        b = float(rhs.split("+")[1].split("(")[0])
        return (lambda T: 10 ** (a + b * (T - 13))), line
    return (lambda T: 220 * 10 ** (0.2 * (T - 9.1))), "old model 220 ppm*10^(0.2(T-9.1))"


def kepler_windows(quarters):
    qs = [int(q) for q in str(quarters).split()] if isinstance(quarters, str) and quarters.strip() else list(KQ)
    return np.array([KQ[q] for q in qs if q in KQ]) - 2167.0  # -> BTJD


def tess_windows(sectors, orb):
    if not isinstance(sectors, str) or not sectors.strip():
        return np.zeros((0, 2))
    s = [int(x) for x in sectors.split()]
    w = orb[orb.Sector.isin(s)][["t0", "t1"]].values
    return w


def any_in(tc, win, margin):
    """tc: (N, K) transit centres; True where any centre lies inside a window shrunk by margin."""
    if len(win) == 0:
        return np.zeros(tc.shape[0], bool)
    lo, hi = win[:, 0] + margin, win[:, 1] - margin
    x = tc[..., None]
    return ((x > lo) & (x < hi)).any(axis=(1, 2))


def epochs_in_range(t0, P, tmin, tmax):
    kmin = np.ceil((tmin - t0) / P).astype(int)
    kmax = np.floor((tmax - t0) / P).astype(int)
    return kmin, kmax


def transit_times(t0, P, tmin, tmax, exclude_t=None):
    """Matrix of transit centres in [tmin, tmax] for each P (padded with nan)."""
    kmin, kmax = epochs_in_range(t0, P, tmin, tmax)
    nmax = int(np.max(kmax - kmin + 1)) if len(P) else 0
    nmax = max(nmax, 1)
    k = kmin[:, None] + np.arange(nmax)[None, :]
    tc = t0 + k * P[:, None]
    tc = np.where(k <= kmax[:, None], tc, np.nan)
    if exclude_t is not None:
        for te in exclude_t:
            tc = np.where(np.abs(tc - te) < 0.5, np.nan, tc)
    return tc


def single_sim(t0, kwin, twin_obs, twin_all, dur_d, prior="geo", Pcat=None, Pup=None, Plo=None, n=40000):
    if prior == "geo":
        # dN/dlnP ∝ P^(-2/3) on [50, 7300] -> sample u = P^(-2/3)
        a, b = 50.0 ** (-2 / 3), 7300.0 ** (-2 / 3)
        P = (b + (a - b) * RNG.random(n)) ** (-3 / 2)
    else:
        s = np.where(RNG.random(n) < 0.5, -np.abs(RNG.normal(0, np.log1p(Plo / Pcat), n)),
                     np.abs(RNG.normal(0, np.log1p(Pup / Pcat), n)))
        P = Pcat * np.exp(s)
        P = P[(P > 20) & (P < 1e5)]
    kmin_t, kmax_t = kwin.min(), kwin.max()
    tc_k = transit_times(t0, P, kmin_t - 1, kmax_t + 1, exclude_t=[t0])
    seen_k = any_in(tc_k, kwin, 0.0)
    P = P[~seen_k]
    if len(P) == 0:
        return dict(p_obs=np.nan, p_all=np.nan, P_min_eff=np.nan, med_alias=np.nan)
    tcO = transit_times(t0, P, 1600, 5200.0)
    hit_obs = any_in(tcO, twin_obs, 0.25 * dur_d)
    hit_all = any_in(tcO, twin_all, 0.25 * dur_d)
    # alias count after a simulated repeat (first TESS hit): P = dT/n allowed by Kepler + TESS flat data
    nal = []
    idx = np.flatnonzero(hit_obs)[:300]
    for i in idx:
        row = tcO[i][np.isfinite(tcO[i])]
        inw = [t for t in row if ((t > twin_obs[:, 0] + 0.25 * dur_d) & (t < twin_obs[:, 1] - 0.25 * dur_d)).any()]
        tT = inw[0]
        dT = tT - t0
        nn = np.arange(1, int(dT / 50) + 1)
        Pn = dT / nn
        okk = ~any_in(transit_times(t0, Pn, kmin_t - 1, kmax_t + 1, exclude_t=[t0]), kwin, 0.0)
        okt = ~any_in(transit_times(t0, Pn, 1600, TODAY_BTJD, exclude_t=[tT]), twin_obs, 0.25 * dur_d)
        nal.append(int((okk & okt).sum()))
    return dict(p_obs=TESS_EFF * hit_obs.mean(), p_all=TESS_EFF * hit_all.mean(), P_min_eff=float(P.min()),
                P_med_post=float(np.median(P)), med_alias=float(np.median(nal)) if nal else np.nan)


def duo_sim(t0, dT, kwin, twin_obs, twin_all, dur_d):
    nmax = int(dT / 10)
    n = np.arange(1, nmax + 1)
    Pn = dT / n
    tck = transit_times(t0, Pn, kwin.min() - 1, kwin.max() + 1, exclude_t=[t0, t0 + dT])
    ok = ~any_in(tck, kwin, 0.0)
    note = ""
    if not ok.any():
        # the catalogue period itself is 'rejected' only because our quarter windows ignore intra-quarter gaps
        ok = n == 1
        note = "n=1 kept: approximate Kepler windows reject every alias"
    Pn = Pn[ok]
    w = Pn ** (-2 / 3)
    w /= w.sum()
    tcO = transit_times(t0, Pn, 1600, 5200.0)
    hit_obs = any_in(tcO, twin_obs, 0.25 * dur_d)
    hit_all = any_in(tcO, twin_all, 0.25 * dur_d)
    return dict(n_alias_kepler=len(Pn), p_obs=TESS_EFF * float((w * hit_obs).sum()),
                p_all=TESS_EFF * float((w * hit_all).sum()),
                aliases=" ".join(f"{p:.2f}" for p in Pn[:8]) + (" ..." if len(Pn) > 8 else ""),
                n_alias_tess_flat=int((~hit_obs).sum()), coverage_note=note)


def status_tables():
    raw = ROOT / "data" / "raw"
    koi = pd.read_csv(HERE / "data" / "raw" / "koi_cumulative.csv")
    ps = pd.read_csv(raw / "nasa_pscomppars.csv")
    ps["tic"] = pd.to_numeric(ps.tic_id.astype(str).str.replace("TIC ", ""), errors="coerce")
    toi = pd.read_csv(raw / "exofop_toi.csv")
    ctoi = pd.read_csv(raw / "exofop_ctoi.csv")
    return koi, ps, toi, ctoi


def main():
    tg = pd.read_csv(HERE / "data" / "kepler_singles_targets.csv")
    sig_fn, noise_line = noise_fn()
    orb = pd.read_csv(ROOT / "singles" / "data" / "tess_orbit_times_latest.csv", comment="#")
    orb["t0"] = pd.to_datetime(orb["Start of Orbit"]).map(lambda x: x.to_julian_date()) - 2457000
    orb["t1"] = pd.to_datetime(orb["End of Orbit"]).map(lambda x: x.to_julian_date()) - 2457000
    orb_obs = orb[orb.t1 < TODAY_BTJD]
    # future sectors that tess-point (max_sector 134) knows pointings for but the orbit-times file lacks:
    # assume contiguous 27.4-d sectors after the last listed one, 1 d lost at each end (approximation)
    last = orb.groupby("Sector").t1.max()
    s0, tend = int(last.index.max()), float(last.max())
    fut = [dict(Sector=s, t0=tend + (s - s0 - 1) * 27.4 + 1.0, t1=tend + (s - s0) * 27.4 - 1.0)
           for s in range(s0 + 1, 135)]
    orb = pd.concat([orb, pd.DataFrame(fut)], ignore_index=True)
    koi, ps, toi, ctoi = status_tables()
    out = []
    for r in tg.itertuples():
        dur = r.dur_h if np.isfinite(r.dur_h) else 15.0
        dur_d = dur / 24
        kwin = kepler_windows(r.quarters)
        tw_obs = tess_windows(r.tess_sectors_observed, orb_obs)
        tw_all = tess_windows(r.tess_sectors_all, orb)
        sig = sig_fn(r.Tmag)
        snr = r.depth_ppm / sig * np.sqrt(min(dur, DUR_CAP_H))
        rec = dict(KIC=r.KIC, TIC=r.TIC, t1_btjd=r.t1_btjd, sigma1h_ppm=sig, snr_tess=snr,
                   dur_used_h=dur, n_kepler_quarters=len(kwin))
        if r.n_kepler_transits == 2 and np.isfinite(r.dT_duo):
            rec.update(duo_sim(r.t1_btjd, r.dT_duo, kwin, tw_obs, tw_all, dur_d))
        else:
            g = single_sim(r.t1_btjd, kwin, tw_obs, tw_all, dur_d, "geo")
            rec.update({f"{k}": v for k, v in g.items()})
            if np.isfinite(r.P_cat_single):
                # use the KM19/FM16/W15/HZW19 central value with a factor-2 lognormal width (their posteriors
                # are wide and asymmetric; only the per-source medians are merged here)
                c = single_sim(r.t1_btjd, kwin, tw_obs, tw_all, dur_d, "cat", r.P_cat_single,
                               r.P_cat_single, r.P_cat_single * 0.5)
                rec.update(p_obs_catprior=c["p_obs"], p_all_catprior=c["p_all"])
        # status
        kk = koi[koi.kepid == r.KIC]
        rec["koi_on_star"] = " ".join(f"{a}:{b[:4]}:{p:.0f}d" for a, b, p in
                                      zip(kk.kepoi_name, kk.koi_disposition, kk.koi_period.fillna(0)))
        pp = ps[ps.tic == r.TIC]
        rec["confirmed_on_star"] = " ".join(f"{a}({p:.0f}d)" for a, p in zip(pp.pl_name, pp.pl_orbper.fillna(0)))
        rec["confirmed_long_P"] = " ".join(f"{a}({p:.0f}d)" for a, p in zip(pp.pl_name, pp.pl_orbper.fillna(0))
                                           if p > 250)
        tt = toi[toi["TIC ID"] == r.TIC]
        rec["TOI"] = " ".join(f"{a}:{d}" for a, d in zip(tt.TOI, tt["TFOPWG Disposition"].fillna("")))
        cc = ctoi[ctoi["TIC ID"] == r.TIC]
        rec["CTOI"] = " ".join(str(a) for a in cc.CTOI)
        out.append(rec)
    f = tg.merge(pd.DataFrame(out), on=["KIC", "TIC", "t1_btjd"])
    # classification
    fp_koi = f.koi_on_star.str.contains("FALS") & ~f.koi_on_star.str.contains("CAND|CONF")
    # U16 put the two KOI-1168 events (KIC 10460629, BKJD 608.26 / 1133.28) in its 'EB-like' group
    u16_eb = (f.KIC == 10460629) & f.sources.eq("U16")
    eb_like = (f.km19_table.fillna("").str.contains("flagged")) | (f.depth_ppm > 25000) | \
        f["flags"].fillna("").str.contains("EB|Binary", case=False) | u16_eb
    f["planet_like"] = ~eb_like
    # literature notes (prior follow-up, see REPORT.txt) that change the status of a candidate
    notes = {5951458: "GOT'EM IV (arXiv:2401.03021): RV P~3378 d, likely M-dwarf companion",
             3239945: "Kepler-167 e, P known (Kepler duo)", 8800954: "Kepler-421 b, P known",
             8410697: "Kepler-1654 b, P known", 7906827: "Kepler-1708 b, P known",
             8631160: "Hansen & Dittmann 2024 checked TESS S14-S74 by eye: no repeat",
             4918810: "Quinn+2021 (arXiv:2107.00027): substellar companion, P ~ 10 yr from RVs"}
    f["lit_note"] = f.KIC.map(notes).fillna("")
    f["period_known"] = (f.n_kepler_transits == 2) & (f.n_alias_kepler == 1)
    f["snr_ok"] = f.snr_tess >= SNR_MIN
    f.to_csv(HERE / "data" / "feasibility.csv", index=False)

    L = [f"Run {date.today()}  (TESS data observed before BTJD {TODAY_BTJD:.1f}); noise: {noise_line}",
         f"candidates: {len(f)} on {f.KIC.nunique()} stars; singles {(f.n_kepler_transits == 1).sum()}, "
         f"duos {(f.n_kepler_transits == 2).sum()}",
         f"TESS-observed (>=1 sector): {(f.n_tess_observed > 0).sum()}; sectors per star median "
         f"{f.drop_duplicates('TIC').n_tess_observed.median():.0f} (range {f.n_tess_observed.min()}-{f.n_tess_observed.max()})",
         f"planet-like (not KM19-flagged, depth<=2.5%, not EB-tagged): {f.planet_like.sum()}",
         f"per-transit TESS SNR >= {SNR_MIN:.0f}: all {f.snr_ok.sum()}, planet-like {(f.snr_ok & f.planet_like).sum()}",
         ""]
    for lab, m in [("singles, planet-like", (f.n_kepler_transits == 1) & f.planet_like),
                   ("singles, planet-like, SNR>=7", (f.n_kepler_transits == 1) & f.planet_like & f.snr_ok),
                   ("singles, all, SNR>=7", (f.n_kepler_transits == 1) & f.snr_ok),
                   ("duos, planet-like", (f.n_kepler_transits == 2) & f.planet_like),
                   ("duos, planet-like, SNR>=7", (f.n_kepler_transits == 2) & f.planet_like & f.snr_ok),
                   ("duos, planet-like, SNR>=7, period NOT already unique", (f.n_kepler_transits == 2) & f.planet_like & f.snr_ok & ~f.period_known)]:
        s = f[m]
        line = (f"{lab:55s} N={len(s):3d}  expected TESS repeats so far: geo-prior {s.p_obs.sum():5.2f}"
                f" | incl. future sectors to S134 (approx.) {s.p_all.sum():5.2f}")
        if "p_obs_catprior" in s:
            line += f" | catalogue-P prior {s.p_obs_catprior.fillna(s.p_obs).sum():5.2f}"
        L.append(line)
    # red-noise check: day-long transits do not gain the full sqrt(duration) in TESS FFI data; cap at 12 h
    f["snr_tess_cap12"] = f.depth_ppm / f.sigma1h_ppm * np.sqrt(np.minimum(f.dur_used_h, 12.0))
    cons = f[(f.n_kepler_transits == 1) & f.planet_like & (f.snr_tess_cap12 >= SNR_MIN)]
    L.append(f"conservative (duration capped at 12 h): planet-like singles with SNR>=7: {len(cons)}, "
             f"expected repeats so far {cons.p_obs.sum():.2f}")
    old = 220 * 10 ** (0.2 * (f.Tmag - 9.1))
    oldsnr = f.depth_ppm / old * np.sqrt(np.minimum(f.dur_used_h, DUR_CAP_H))
    L.append(f"old noise model (220 ppm*10^(0.2(T-9.1))): SNR>=7 all {int((oldsnr >= 7).sum())}, planet-like singles "
             f"{int(((oldsnr >= 7) & f.planet_like & (f.n_kepler_transits == 1)).sum())}")
    f.to_csv(HERE / "data" / "feasibility.csv", index=False)
    L.append("\nPer source catalogue (a candidate counts for every catalogue that lists it):")
    L.append(f"{'source':7s} {'N':>4s} {'singles':>7s} {'duos':>5s} {'planet-like':>11s} {'SNR>=7':>6s} "
             f"{'pl+SNR':>6s} {'E[repeat] all':>13s} {'E[repeat] pl+SNR':>16s} {'new-P-able pl+SNR':>17s}")
    for src in ["KM19", "W15", "FM16", "HZW19", "U16", "HC24"]:
        s = f[f.sources.str.contains(src)]
        good = s[s.planet_like & s.snr_ok]
        newp = good[~good.period_known & (good.lit_note == "")]
        L.append(f"{src:7s} {len(s):4d} {(s.n_kepler_transits == 1).sum():7d} {(s.n_kepler_transits == 2).sum():5d} "
                 f"{s.planet_like.sum():11d} {s.snr_ok.sum():6d} {len(good):6d} {s.p_obs.sum():13.2f} "
                 f"{good.p_obs.sum():16.2f} {newp.p_obs.sum():17.2f}")
    sing = f[(f.n_kepler_transits == 1) & f.planet_like & f.snr_ok]
    L.append(f"\nsingles (planet-like, SNR>=7): median aliases left after a simulated repeat = "
             f"{sing.med_alias.median():.0f}; median effective P_min = {sing.P_min_eff.median():.0f} d")
    duo = f[(f.n_kepler_transits == 2)]
    L.append(f"duos: Kepler coverage leaves 1 alias (period already known) for {(duo.n_alias_kepler == 1).sum()}"
             f"/{len(duo)}; >1 alias for {(duo.n_alias_kepler > 1).sum()}")
    cols = ["KIC", "TIC", "sources", "n_kepler_transits", "t1_bkjd", "dT_duo", "n_alias_kepler", "depth_ppm",
            "dur_used_h", "Tmag", "snr_tess", "n_tess_observed", "p_obs", "p_all", "planet_like", "koi_on_star",
            "confirmed_long_P", "TOI", "CTOI", "lit_note"]
    L.append("\nTop candidates by expected repeat probability (planet-like, SNR>=7):")
    top = f[f.planet_like & f.snr_ok].sort_values("p_obs", ascending=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 40):
        L.append(top[cols].round(2).head(40).to_string(index=False))
    txt = "\n".join(L)
    (HERE / "data" / "feasibility_summary.txt").write_text(txt + "\n", encoding="utf-8")
    print(txt)


if __name__ == "__main__":
    main()
