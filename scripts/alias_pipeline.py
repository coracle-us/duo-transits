"""
Plan step 3: recompute the period aliases of each duotransit from the actual TESS coverage.

For a candidate with transits at Tc1, Tc2 (dT = Tc2 - Tc1), alias n has period P_n = dT / n.
For every n (P_n >= P_FLOOR) we predict the epochs T_k = Tc1 + k * P_n (k != 0, n) that fall
inside the light curves, and test each predicted window:

  1. Local baseline: a low-order polynomial fitted to out-of-transit data within
     +/- (D/2 + 3 sigma_T + BASE_W) of the prediction (in-transit window masked).
  2. Depth: mean depression of the flat part of the transit (|t - T| < FLAT * D/2), scanned
     over timing offsets within +/- 3 sigma_T (sigma_T propagated from the Tc errors).
  3. Noise: white scatter of the baseline residuals / sqrt(N_in), inflated by a red-noise
     factor beta measured by binning the residuals on the in-transit timescale.

Per-epoch verdict (delta = expected depth from Tc1/Tc2, a conservative threshold
DEPTH_FRAC * delta is what we test against):
  excluded       at EVERY trial offset: enough in-transit coverage AND d + K_SIG*sigma < DEPTH_FRAC*delta
                 AND |d| < K_SIG*sigma (flat: no significant dip or brightening)
                 AND beta <= BETA_MAX_EXCL
  supported      at SOME offset: d > DEPTH_FRAC*delta and d/sigma > SNR_MIN
  unconstrained  not enough data in the window at some offset (could hide in a gap)
  ambiguous      covered, but neither excluded nor supported
Per-alias verdict (see summarize_aliases): excluded (a flat window; 'conflict' if it also has
informative support), excluded-by-transits (misses a corroborated detected dip), supported,
or allowed (not ruled out). Survivors = supported + allowed.

The known transits (k = 0, n) are tested the same way as a built-in sanity check: they
should come out 'supported' whenever the data covers them.

Light curves: one product per sector, preferring SPOC 2-min PDCSAP > TESS-SPOC PDCSAP >
QLP SAP (local detrending makes SAP adequate) > SPOC 20-s. quality == 0 only.

Usage:
  python scripts/alias_pipeline.py TIC [TIC ...]          # all data
  python scripts/alias_pipeline.py --max-sector 39 TIC     # paper-era data only
Outputs results/aliases/<label>/<TIC>_epochs.csv and <TIC>_aliases.csv
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
LC = ROOT / "data" / "raw" / "lightcurves"
CAT = ROOT / "data" / "catalog"
RES = ROOT / "results" / "aliases"

P_FLOOR = 5.0          # shortest period tested (d)
DEFAULT_TC_ERR = 0.005  # d, when the catalog has no Tc uncertainty
BASE_W = 0.75          # d of baseline either side of the masked window
FLAT = 0.7             # fraction of the duration treated as in-transit ("flat bottom")
DEPTH_FRAC = 0.5       # test against half the expected depth (conservative)
K_SIG = 3.0            # sigma margin for excluded/supported
SNR_MIN = 5.0          # detection significance required for 'supported'
MIN_COVER = 0.3        # min fraction of the flat window that must have data
MIN_IN = 3             # min number of in-transit points
BETA_MAX_EXCL = 3.0    # no exclusion from windows with red-noise factor above this
# A "detected transit" can exclude other aliases, so it must be robust:
DETECT_SNR = 8.0       # SNR for a supported epoch to count as a detected transit
DEPTH_MAX_FRAC = 2.0   # depth between DEPTH_FRAC and this times the expected depth
DETECT_COVER = 0.8     # in-transit coverage
DETECT_SIDE = 0.5      # baseline on BOTH sides, >= this fraction of a full BASE_W of cadences
PRODUCT_RANK = ["SPOC_120", "TESS-SPOC", "QLP", "SPOC_20"]
SEG_GAP = 0.5          # d; a gap longer than this starts a new data segment
QLP_TRIM = 0.5         # d trimmed from the start of each QLP segment (raw SAP ramps; newer
                       # QLP files carry no detrended KSPSAP flux)


# ------------------------------------------------------------------ light curves
def load_lightcurves(tic, max_sector=None, products=PRODUCT_RANK):
    """Best product per sector, quality==0, finite, normalised by the sector median."""
    files = {}
    for f in (LC / str(tic)).glob("s*.csv.gz"):
        m = re.fullmatch(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", f.name)
        if not m:
            continue
        sec, author, exp = int(m.group(1)), m.group(2), int(m.group(3))
        if sec > 200 or (max_sector is not None and sec > max_sector):
            continue
        key = f"SPOC_{exp}" if author == "SPOC" else author
        files.setdefault(sec, {})[key] = f
    parts = []
    for sec, avail in sorted(files.items()):
        key = next((k for k in products if k in avail), None)
        if key is None:
            continue
        df = pd.read_csv(avail[key], usecols=["time", "flux", "flux_err", "quality"])
        df = df[np.isfinite(df.time) & np.isfinite(df.flux) & (df.quality == 0)]
        if len(df) < 50:
            continue
        df = df.sort_values("time")
        seg = np.concatenate([[0], np.cumsum(np.diff(df.time.values) > SEG_GAP)])
        if key == "QLP":
            seg_start = pd.Series(df.time.values).groupby(seg).transform("min").values
            keep = df.time.values - seg_start >= QLP_TRIM
            df, seg = df[keep], seg[keep]
        med = np.nanmedian(df.flux)
        err = df.flux_err / med if np.isfinite(df.flux_err).any() else np.nan
        parts.append(pd.DataFrame({"time": df.time.values, "flux": df.flux.values / med,
                                   "err": np.asarray(err, dtype=float), "sector": sec, "product": key,
                                   "segment": sec * 1000 + seg}))
    if not parts:
        return pd.DataFrame(columns=["time", "flux", "err", "sector", "product", "segment"])
    return pd.concat(parts, ignore_index=True).sort_values("time").reset_index(drop=True)


# ------------------------------------------------------------------ one window
def test_window(lc_t, lc_f, T, sigT, D, delta, lc_seg=None):
    """Test one predicted transit time. lc_t/lc_f/lc_seg: arrays (sorted). Returns dict.
    The baseline is a polynomial in time plus a free offset for each data segment, so jumps
    across data gaps are not mistaken for dips; segments without baseline points are dropped."""
    half = D / 2
    mask_half = half + 3 * sigT + 0.02          # region masked from the baseline fit
    lo, hi = T - mask_half - BASE_W, T + mask_half + BASE_W
    i0, i1 = np.searchsorted(lc_t, [lo, hi])
    t, f = lc_t[i0:i1], lc_f[i0:i1]
    sg = lc_seg[i0:i1] if lc_seg is not None else np.zeros(len(t), dtype=int)
    oot = np.abs(t - T) > mask_half
    anchored = np.isin(sg, np.unique(sg[oot])) if oot.any() else np.zeros(len(t), bool)
    t, f, sg, oot = t[anchored], f[anchored], sg[anchored], oot[anchored]
    res = dict(T_pred=T, sigma_T=sigT, n_window=len(t))
    if len(t) == 0:
        return {**res, "verdict": "no-data"}
    n_left, n_right = np.sum(oot & (t < T)), np.sum(oot & (t > T))
    res.update(n_oot_left=int(n_left), n_oot_right=int(n_right))
    if n_left + n_right < 10:
        # in-window data but no usable baseline
        flat_any = np.sum(np.abs(t - T) < FLAT * half)
        return {**res, "verdict": "unconstrained" if flat_any else "no-data"}
    deg = 2 if (n_left >= 10 and n_right >= 10) else (1 if min(n_left, n_right) >= 5 else 0)
    segs = np.unique(sg)
    X = np.column_stack([(t - T) ** p for p in range(1, deg + 1)]
                        + [(sg == u).astype(float) for u in segs]) if (deg or len(segs)) else None
    coef, *_ = np.linalg.lstsq(X[oot], f[oot], rcond=None)
    model = X @ coef
    r = f / model - 1.0
    r_oot = r[oot]
    sig_w = 1.4826 * np.median(np.abs(r_oot - np.median(r_oot)))  # robust white scatter
    cad = np.median(np.diff(t)) if len(t) > 1 else 0.0
    n_flat_full = max(1, int(round(FLAT * D / cad))) if cad > 0 else 1
    beta = red_noise_beta(r_oot, n_flat_full)

    offsets = np.unique(np.concatenate([[0.0], np.arange(-3 * sigT, 3 * sigT + 1e-9,
                                                         max(D / 10, 1e-3))]))
    best = None
    all_excl, any_uncon = True, False
    for tau in offsets:
        inn = np.abs(t - (T + tau)) < FLAT * half
        n_in = int(inn.sum())
        cover = n_in * cad / (FLAT * D) if cad > 0 else 0.0
        if n_in < MIN_IN or cover < MIN_COVER:
            any_uncon = True
            all_excl = False
            continue
        d = -np.mean(r[inn])
        s = beta * sig_w / np.sqrt(n_in)
        # Absence of a transit is only credible in a clean, FLAT window: low red noise, no
        # significant brightening (systematics, e.g. an edge ramp) and no significant dip
        # (a real but shallower/shifted event is evidence of something, not of absence).
        excl = ((d + K_SIG * s < DEPTH_FRAC * delta) and (abs(d) < K_SIG * s)
                and (beta <= BETA_MAX_EXCL))
        supp = (d > DEPTH_FRAC * delta) and (d / s > SNR_MIN)
        all_excl &= excl
        cand = dict(offset=tau, depth=d, depth_err=s, snr=d / s, n_in=n_in, cover=min(cover, 1.0),
                    supp=supp)
        # keep the most transit-like offset (highest SNR) as the reported one
        if best is None or cand["snr"] > best["snr"]:
            best = cand
    res.update(beta=beta, sigma_white=sig_w, baseline_deg=deg, cadence=cad)
    if best is None:
        return {**res, "verdict": "unconstrained"}
    res.update({k: v for k, v in best.items() if k != "supp"})
    if best["supp"]:
        v = "supported"
    elif all_excl and not any_uncon:
        v = "excluded"
    elif any_uncon:
        v = "unconstrained"
    else:
        v = "ambiguous"
    return {**res, "verdict": v}


def red_noise_beta(r, n_bin):
    """Ratio of binned scatter to white-noise expectation (>= 1)."""
    if n_bin <= 1 or len(r) < 3 * n_bin:
        return 1.0
    nb = len(r) // n_bin
    rb = r[: nb * n_bin].reshape(nb, n_bin).mean(axis=1)
    sw = np.std(r, ddof=1)
    if sw == 0 or nb < 3:
        return 1.0
    return float(max(1.0, np.std(rb, ddof=1) / (sw / np.sqrt(n_bin)) * np.sqrt(nb / (nb - 1))))


# ------------------------------------------------------------------ one candidate
def run_candidate(c, lc, p_floor=P_FLOOR):
    """c: row of candidates.csv. lc: output of load_lightcurves. Returns (epochs, aliases)."""
    tc1, tc2, dT = c.Tc1, c.Tc2, c.Tc2 - c.Tc1
    s1 = c.Tc1_err if np.isfinite(c.Tc1_err) else DEFAULT_TC_ERR
    s2 = c.Tc2_err if np.isfinite(c.Tc2_err) else DEFAULT_TC_ERR
    D = c.Dur_h / 24.0
    delta = 1e-3 * np.nanmean([c.Depth1_ppt, c.Depth2_ppt])
    t, f = lc.time.to_numpy(), lc.flux.to_numpy()
    sec = lc.sector.to_numpy()
    sg_all = lc.segment.to_numpy()
    if len(t) == 0:
        return pd.DataFrame(), pd.DataFrame()
    tmin, tmax = t.min() - 1, t.max() + 1
    n_max = int(np.floor(dT / p_floor))
    ep_rows = []
    for n in range(1, n_max + 1):
        P = dT / n
        k_lo = int(np.ceil((tmin - tc1) / P))
        k_hi = int(np.floor((tmax - tc1) / P))
        for k in range(k_lo, k_hi + 1):
            T = tc1 + k * P
            # does any data lie near this epoch? (cheap pre-check)
            j = np.searchsorted(t, T)
            near = min(abs(t[min(j, len(t) - 1)] - T), abs(t[max(j - 1, 0)] - T))
            if near > D / 2 + 0.1:
                continue
            u = k / n
            sigT = np.hypot((1 - u) * s1, u * s2)
            r = test_window(t, f, T, sigT, D, delta, sg_all)
            jj = min(j, len(t) - 1)
            r.update(n=n, P=P, k=k, known=(k in (0, n)), sector=int(sec[jj]))
            ep_rows.append(r)
    ep = pd.DataFrame(ep_rows)
    return ep, summarize_aliases(ep, c, n_max, D, s1, s2)


def summarize_aliases(ep, c, n_max, D, s1, s2):
    """Combine epoch verdicts into one verdict per alias.

    - Epochs at Tc1 + m*dT (k a multiple of n) are predicted by EVERY alias, so support there
      is uninformative ('common'); only other epochs count as informative support.
    - Flat windows (epoch verdict 'excluded') rule an alias out: 'excluded'.
    - Detected dips: supported epochs passing the stricter DETECT_* cuts (SNR, depth,
      coverage, baseline on both sides), merged by time, excluding the Tc1 + m*dT grid.
    - Corroboration: a dip is corroborated if it lies on the ephemeris of an alias that
      survives the flat-window test AND that ephemeris also passes through at least one other
      detected dip (two dips aligning by chance with one alias is ~1% vs ~40% for one).
      Aliases that miss a corroborated dip are 'excluded-by-transits' (they cannot explain an
      observed transit, e.g. sub-harmonics at 2x, 3x the true period).
    - Uncorroborated dips never exclude; they are listed in `uncorroborated_dips` for vetting
      (systematics, a second planet, or an eclipsing binary).
    """
    dT = c.Tc2 - c.Tc1
    delta = 1e-3 * np.nanmean([c.Depth1_ppt, c.Depth2_ppt])
    detected = np.array([])
    if len(ep):
        ep["common"] = (ep.k % ep.n) == 0
        side_min = DETECT_SIDE * BASE_W / ep.cadence
        det = ep[(ep.verdict == "supported") & (ep.snr >= DETECT_SNR) & ~ep.known
                 & (ep.depth <= DEPTH_MAX_FRAC * delta) & (ep.cover >= DETECT_COVER)
                 & (ep.n_oot_left >= side_min) & (ep.n_oot_right >= side_min)]
        merged = []
        for tt in np.sort((det.T_pred + det.offset).to_numpy()):
            if not merged or tt - merged[-1] > D:
                merged.append(tt)
        # drop dips on the common Tc1 + m*dT grid (they discriminate nothing)
        detected = np.array([tt for tt in merged
                             if abs(tt - (c.Tc1 + np.round((tt - c.Tc1) / dT) * dT)) > D])

    def on_grid(tt, n):
        P = dT / n
        k = np.round((tt - c.Tc1) / P)
        u = k / n
        return abs(tt - (c.Tc1 + k * P)) <= D / 2 + 3 * np.hypot((1 - u) * s1, u * s2)

    per = {}
    for n in range(1, n_max + 1):
        e = ep[(ep.n == n) & ~ep.known] if len(ep) else pd.DataFrame(
            columns=["verdict", "common", "T_pred", "sector"])
        per[n] = e
    flat_ok = [n for n in per if not (per[n].verdict == "excluded").any()]
    corroborated = set()
    for n in flat_ok:
        hits = [tt for tt in detected if on_grid(tt, n)]
        if len(hits) >= 2:
            corroborated.update(hits)
    uncorroborated = [tt for tt in detected if tt not in corroborated]

    rows = []
    for n, e in per.items():
        inf = e[~e.common] if len(e) else e
        n_ex = int((e.verdict == "excluded").sum())
        n_su_inf = int((inf.verdict == "supported").sum())
        missed = [tt for tt in sorted(corroborated) if not on_grid(tt, n)]
        if n_ex and n_su_inf:
            verdict = "conflict"
        elif n_ex:
            verdict = "excluded"
        elif missed:
            verdict = "excluded-by-transits"
        elif n_su_inf:
            verdict = "supported"
        else:
            verdict = "allowed"
        fmt = lambda d: " ".join(f"{T:.2f}(s{sec})" for T, sec in zip(d.T_pred, d.sector))
        rows.append(dict(
            TIC=c.TIC, n=n, P=dT / n, verdict=verdict, n_excluded=n_ex, n_supported_informative=n_su_inf,
            n_supported_common=int(((e.verdict == "supported") & e.common).sum()) if len(e) else 0,
            n_ambiguous=int((e.verdict == "ambiguous").sum()),
            n_unconstrained=int((e.verdict == "unconstrained").sum()),
            n_informative_tested=int(inf.verdict.isin(["excluded", "supported", "ambiguous"]).sum()),
            excluded_by=fmt(e[e.verdict == "excluded"]) if n_ex else "",
            misses_corroborated=" ".join(f"{tt:.2f}" for tt in missed),
            uncorroborated_dips=" ".join(f"{tt:.2f}" for tt in uncorroborated),
            supported_at=fmt(inf[inf.verdict == "supported"]) if n_su_inf else ""))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tics", nargs="*", type=int, help="TIC IDs (default: all 85)")
    ap.add_argument("--max-sector", type=int, default=None,
                    help="use only sectors <= this (plus each candidate's own Tc sectors)")
    ap.add_argument("--label", default=None, help="output subdirectory (default: all / maxsNN)")
    args = ap.parse_args()
    cand = pd.read_csv(CAT / "candidates.csv").set_index("TIC", drop=False)
    tics = args.tics or list(cand.index)
    label = args.label or ("all" if args.max_sector is None else f"maxs{args.max_sector}")
    out = RES / label
    out.mkdir(parents=True, exist_ok=True)
    for tic in tics:
        c = cand.loc[tic]
        ms = args.max_sector
        lc = load_lightcurves(tic, max_sector=max(ms, int(c.S2)) if ms is not None else None)
        ep, al = run_candidate(c, lc)
        ep.to_csv(out / f"{tic}_epochs.csv", index=False)
        al.to_csv(out / f"{tic}_aliases.csv", index=False)
        vc = al.verdict.value_counts().to_dict() if len(al) else {}
        # verdict on the two known transits (identical for every alias; take n=1)
        known = ep[ep.known & (ep.n == 1)].verdict.tolist() if len(ep) else []
        print(f"TIC {tic}: {len(lc)} pts, {lc.sector.nunique()} sectors | aliases {vc} | "
              f"known transits (Tc1, Tc2): {known}")


if __name__ == "__main__":
    main()
