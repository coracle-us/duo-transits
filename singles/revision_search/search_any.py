"""
Unconstrained dip search (RNAAS AAS82207 revision, editor's question Q5).

The note's repeat search looked only for events of the CATALOGUED depth and duration. This
script searches the same light curves for ANY transit- or eclipse-like dip, with a free depth
and a grid of durations, and runs the identical search on inverted light curves (2 - flux) as a
false-alarm control.

Data (same as singles/scripts/pilot.py): data/raw/lightcurves/<TIC>/s<NNN>_<product>_<cad>s.csv.gz,
products SPOC 2-min, TESS-SPOC, QLP and our own TESScut photometry (OWN, sectors >= 100 without
pipeline products). quality == 0 only; QLP segments trimmed by 0.5 d at the start (as
alias_pipeline.load_lightcurves). SPOC 20-s is used only if a sector has nothing else.
The search runs on the best product per sector (pipeline ranking); every detection is then
re-measured in each other product of that sector.

RULES (fixed before any result was looked at):
  durations   1, 2, 4, 8, 16, 24 h, plus the catalogued duration D_cat if not within 10% of these
  fast pass   per data segment (gap > 0.5 d splits), box of width 0.7 D on a D/4 grid against a
              linear baseline from the two flanks (each max(0.75 d, 2 D) long); white-noise SNR
              >= 5 local maxima go to the full test (at most 1500 per duration per star)
  full test   within the event's own segment: quadratic baseline (linear if a side has < 10
              points) fitted out of |t - T| > D/2 + D/4 + 0.02 d over a window of max(0.75, 2D) d
              each side; centre scanned over +/- D/4; depth = mean residual over |t - T| < 0.35 D;
              sigma = beta * sigma_white / sqrt(n_in), beta = red-noise factor at the in-event
              timescale (must be measurable: >= 3 bins in the baseline)
  DETECTION   SNR >= 8 (red-noise inflated); in-event coverage >= 0.8; n_in >= 3;
              baseline on both sides >= max(0.25 d, 0.5 D) of data in the same segment;
              event edges >= 0.25 d from the segment (sector/orbit/gap) edges;
              momentum dumps (quality bit 32, also bits 1 and 4; flags from every product of the
              star and its TESScut cutouts): v1 rejected any event with a dump within D/2 + 0.05 d
              of the centre. That vetoed 37 of the 272 catalogued events (up to SNR 76), because
              long events often straddle a dump, so v2 (used for all results) rejects an event
              only if a dump lies in the in-event core (|t - T| < 0.35 D + 0.02 d) and the core
              data on either side of the dump do not BOTH show >= 0.5x the depth (>= 3 points
              each). Every event records dump_in_window, so the v1 rule can be re-applied;
              median in-event depth >= 0.5 x mean depth (not driven by a few points);
              depth >= 3 x the local 1-h binned scatter
  merging     detections whose [T - D/2, T + D/2] windows overlap are one event (highest SNR kept)
  2nd product for each event and each other product covering it (SPOC-120, TESS-SPOC, QLP, OWN):
              same test at the same D (centre scanned +/- D/8); 'confirmed' if SNR >= 4 and depth
              within 0.5-2x of the primary depth
  control     the identical search on 2 - flux (primary and alternative products)

Usage:
  python singles/revision_search/search_any.py                  # 306 classified stars
  python singles/revision_search/search_any.py --tics 159614322 # subset
  python singles/revision_search/search_any.py --ebs            # Prsa+2022 EB positive controls
Outputs: singles/revision_search/data/per_star/<TIC>.json (or data/eb_per_star/)
"""
import argparse
import glob
import json
import re
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LC_DIR = ROOT / "data" / "raw" / "lightcurves"
TESSCUT_DIR = ROOT / "data" / "raw" / "tesscut"
RES = ROOT / "singles" / "results"

DURS_H = [1, 2, 4, 8, 16, 24]
RANK = ["SPOC_120", "TESS-SPOC", "QLP", "SPOC_20", "OWN"]  # alias_pipeline.PRODUCT_RANK + OWN
ALT_PRODUCTS = ["SPOC_120", "TESS-SPOC", "QLP", "OWN"]
SPOC_FAMILY = {"SPOC_120", "TESS-SPOC", "SPOC_20"}  # same pixels/pipeline family
SEG_GAP, QLP_TRIM = 0.5, 0.5
DUMP_BITS = 32 | 1 | 4
FAST_SNR, MAX_FAST = 5.0, 1500
FLAT = 0.7
DET_SNR, DET_COVER, DET_NIN = 8.0, 0.8, 3
EDGE = 0.25
DUMP_PAD = 0.05
MIN_DEPTH_1H = 3.0
CONF_SNR, CONF_LO, CONF_HI = 4.0, 0.5, 2.0


# ------------------------------------------------------------------ data
def load_products(tic, lc_dir=LC_DIR):
    """All products per sector. Returns (prods, dumps): prods[sector][product] = dict(t, f, seg),
    dumps = sorted array of flagged (dump/attitude) times from every file and TESScut cutout."""
    files = {}
    for fn in (Path(lc_dir) / str(tic)).glob("s*.csv.gz"):
        m = re.fullmatch(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", fn.name)
        if not m:
            continue
        sec, author, exp = int(m.group(1)), m.group(2), int(m.group(3))
        if sec > 200:
            continue
        key = f"SPOC_{exp}" if author == "SPOC" else author
        if key not in RANK:
            continue
        files.setdefault(sec, {})[key] = fn
    prods, dumps = {}, []
    for sec, avail in sorted(files.items()):
        for key, fn in avail.items():
            if key == "SPOC_20" and len(set(avail) - {"SPOC_20"}) > 0:
                continue
            df = pd.read_csv(fn, usecols=["time", "flux", "quality"])
            df = df[np.isfinite(df.time)]
            q = df.quality.fillna(0).astype(int).to_numpy()
            dumps.append(df.time.to_numpy()[(q & DUMP_BITS) > 0])
            df = df[np.isfinite(df.flux) & (df.quality == 0)].sort_values("time")
            if len(df) < 50:
                continue
            t, f = df.time.to_numpy(float), df.flux.to_numpy(float)
            seg = np.concatenate([[0], np.cumsum(np.diff(t) > SEG_GAP)])
            if key == "QLP":
                s0 = pd.Series(t).groupby(seg).transform("min").to_numpy()
                keep = t - s0 >= QLP_TRIM
                t, f, seg = t[keep], f[keep], seg[keep]
            if len(t) < 50:
                continue
            f = f / np.nanmedian(f)
            prods.setdefault(sec, {})[key] = dict(t=t, f=f, seg=sec * 1000 + seg)
    for fn in glob.glob(str(TESSCUT_DIR / str(tic) / "tesscut" / "*.fits")):
        try:
            from astropy.io import fits
            with fits.open(fn) as h:
                q, tt = h[1].data["QUALITY"], h[1].data["TIME"]
            dumps.append(np.asarray(tt, float)[(np.asarray(q, int) & DUMP_BITS) > 0])
        except Exception:  # noqa: BLE001
            pass
    d = np.concatenate(dumps) if dumps else np.array([])
    return prods, np.unique(np.round(d[np.isfinite(d)], 4))


def primary(prods):
    """Best product per sector -> concatenated arrays + per-sector product name."""
    parts, used = [], {}
    for sec, avail in sorted(prods.items()):
        key = next((k for k in RANK if k in avail), None)
        if key is None:
            continue
        used[sec] = key
        parts.append((avail[key]["t"], avail[key]["f"], avail[key]["seg"], np.full(len(avail[key]["t"]), sec)))
    if not parts:
        return None, used
    t = np.concatenate([p[0] for p in parts])
    o = np.argsort(t)
    return dict(t=t[o], f=np.concatenate([p[1] for p in parts])[o], seg=np.concatenate([p[2] for p in parts])[o],
                sec=np.concatenate([p[3] for p in parts])[o]), used


# ------------------------------------------------------------------ tests
def red_noise_beta(r, n_bin):
    """Binned/white scatter ratio (>= 1); NaN if fewer than 3 bins (k2lib.red_noise_beta_k2)."""
    if n_bin <= 1:
        return 1.0
    nb = len(r) // n_bin
    if nb < 3:
        return np.nan
    rb = r[: nb * n_bin].reshape(nb, n_bin).mean(axis=1)
    sw = np.std(r, ddof=1)
    if sw == 0:
        return np.nan
    return float(max(1.0, np.std(rb, ddof=1) / (sw / np.sqrt(n_bin)) * np.sqrt(nb / (nb - 1))))


def fast_pass(t, f, D):
    """White-noise box SNR on a D/4 grid inside one segment; returns candidate centres."""
    if len(t) < 20 or t[-1] - t[0] < D + 2 * EDGE:
        return np.array([])
    T = np.arange(t[0] + D / 2 + EDGE, t[-1] - D / 2 - EDGE + 1e-9, D / 4)
    if len(T) == 0:
        return T
    B = max(0.75, 2 * D)
    cf, ct = np.concatenate([[0], np.cumsum(f)]), np.concatenate([[0], np.cumsum(t)])

    def box(a, b):
        i, j = np.searchsorted(t, a), np.searchsorted(t, b)
        n = j - i
        with np.errstate(invalid="ignore", divide="ignore"):
            return n, (cf[j] - cf[i]) / n, (ct[j] - ct[i]) / n

    n_in, m_in, _ = box(T - FLAT * D / 2, T + FLAT * D / 2)
    mh = D / 2 + 0.02
    nL, mL, tL = box(T - mh - B, T - mh)
    nR, mR, tR = box(T + mh, T + mh + B)
    ok = (n_in >= 2) & (nL >= 3) & (nR >= 3)
    with np.errstate(invalid="ignore", divide="ignore"):
        base = mL + (mR - mL) * (T - tL) / (tR - tL)
        d = base - m_in
        sp = 1.4826 * np.median(np.abs(np.diff(f) - np.median(np.diff(f)))) / np.sqrt(2)
        snr = d / (sp * np.sqrt(1 / n_in + 0.25 / nL + 0.25 / nR))
    snr = np.where(ok & np.isfinite(snr), snr, 0.0)
    pk, _ = find_peaks(np.concatenate([[0], snr, [0]]), height=FAST_SNR, distance=2)
    pk = pk - 1
    pk = pk[np.argsort(-snr[pk])][:MAX_FAST]
    return T[pk]


def test_free(t, f, T, D, scan):
    """Free-depth test of a dip of duration D near T (scan = +/- range of centre, d).
    t, f: data of ONE segment (sorted). Returns dict or None."""
    B = max(0.75, 2 * D)
    mh = D / 2 + scan + 0.02
    i0, i1 = np.searchsorted(t, [T - mh - B, T + mh + B])
    tt, ff = t[i0:i1], f[i0:i1]
    if len(tt) < 20:
        return None
    oot = np.abs(tt - T) > mh
    nL, nR = int(np.sum(oot & (tt < T))), int(np.sum(oot & (tt > T)))
    if nL < 5 or nR < 5:
        return None
    deg = 2 if (nL >= 10 and nR >= 10) else 1
    x = tt - T
    c = np.polyfit(x[oot], ff[oot], deg)
    r = ff / np.polyval(c, x) - 1.0
    ro = r[oot]
    sw = 1.4826 * np.median(np.abs(ro - np.median(ro)))
    cad = float(np.median(np.diff(tt)))
    if sw <= 0 or cad <= 0:
        return None
    beta = red_noise_beta(ro, max(1, int(round(FLAT * D / cad))))
    b = 1.0 if np.isnan(beta) else beta
    nb1 = max(1, int(round((1 / 24) / cad)))
    k = len(ro) // nb1
    if k >= 5:
        rb = ro[: k * nb1].reshape(k, nb1).mean(axis=1)
        s1h = 1.4826 * np.median(np.abs(rb - np.median(rb)))
    else:
        s1h = np.nan
    taus = np.arange(-scan, scan + 1e-9, max(D / 20, cad / 2)) if scan > 0 else np.array([0.0])
    best = None
    for tau in taus:
        inn = np.abs(tt - (T + tau)) < FLAT * D / 2
        n_in = int(inn.sum())
        if n_in < 2:
            continue
        d = -float(np.mean(r[inn]))
        snr = d / (b * sw / np.sqrt(n_in))
        if best is None or snr > best["snr"]:
            best = dict(T=float(T + tau), depth=d, depth_med=-float(np.median(r[inn])), snr=float(snr), n_in=n_in,
                        cover=float(min(1.0, n_in * cad / (FLAT * D))), _tin=tt[inn], _rin=r[inn])
    if best is None:
        return None
    Tc = best["T"]
    side_L = float(np.sum(oot & (tt < Tc - D / 2)) * cad)
    side_R = float(np.sum(oot & (tt > Tc + D / 2)) * cad)
    best.update(D=D, beta=float(b), beta_measured=bool(not np.isnan(beta)), sigma_white=float(sw),
                sigma_1h=float(s1h), cadence=cad, side_L=side_L, side_R=side_R, seg_start=float(t[0]),
                seg_end=float(t[-1]))
    return best


def passes(r, dumps):
    """Detection rules (see module docstring). Returns (bool, reason)."""
    D = r["D"]
    if r["snr"] < DET_SNR:
        return False, "snr"
    if not r["beta_measured"]:
        return False, "beta-unmeasured"
    if r["cover"] < DET_COVER or r["n_in"] < DET_NIN:
        return False, "coverage"
    side_min = max(0.25, 0.5 * D)
    if r["side_L"] < side_min or r["side_R"] < side_min:
        return False, "baseline-side"
    if r["T"] - D / 2 - r["seg_start"] < EDGE or r["seg_end"] - (r["T"] + D / 2) < EDGE:
        return False, "edge"
    if r["depth_med"] < 0.5 * r["depth"]:
        return False, "outlier-driven"
    if not np.isfinite(r["sigma_1h"]) or r["depth"] < MIN_DEPTH_1H * r["sigma_1h"]:
        return False, "depth-vs-scatter"
    dd = dumps[np.abs(dumps - r["T"]) < D / 2 + DUMP_PAD] if len(dumps) else []
    r["dump_in_window"] = bool(len(dd))
    for td in dd:
        # v2 (after the recovery check, see REPORT): a dump is only fatal if it lies in the
        # in-event core and the core data on either side of it do not both show the dip
        if abs(td - r["T"]) >= FLAT * D / 2 + 0.02:
            continue
        tin, rin = r["_tin"], r["_rin"]
        for side in (rin[tin < td], rin[tin > td]):
            if len(side) < 3 or -np.mean(side) < 0.5 * r["depth"]:
                return False, "momentum-dump"
    return True, "ok"


def segments(lc):
    out = {}
    for u in np.unique(lc["seg"]):
        m = lc["seg"] == u
        out[int(u)] = (lc["t"][m], lc["f"][m])
    return out


def search(lc, used, dumps, durs, sign=1.0):
    """Full search on the primary light curve; sign=-1 searches 2 - flux (control).
    Returns (merged events, counts of near-misses by reason)."""
    segs = segments(lc)
    dets, why = [], {}
    for D in durs:
        for u, (t, f) in segs.items():
            ff = f if sign > 0 else 2.0 - f
            for T in fast_pass(t, ff, D):
                r = test_free(t, ff, T, D, scan=D / 4)
                if r is None:
                    continue
                ok, reason = passes(r, dumps)
                r.pop("_tin", None)
                r.pop("_rin", None)
                if r["snr"] >= DET_SNR:
                    why[reason] = why.get(reason, 0) + 1
                if ok:
                    sec = int(u // 1000)
                    r.update(sector=sec, product=used[sec], segment=int(u))
                    dets.append(r)
    return merge(dets), why


def merge(dets):
    dets = sorted(dets, key=lambda r: -r["snr"])
    out = []
    for r in dets:
        hit = None
        for m in out:
            if abs(r["T"] - m["T"]) < (r["D"] + m["D"]) / 2:
                hit = m
                break
        if hit is None:
            r = dict(r)
            r["durs_detected"] = [round(r["D"] * 24, 2)]
            out.append(r)
        elif round(r["D"] * 24, 2) not in hit["durs_detected"]:
            hit["durs_detected"].append(round(r["D"] * 24, 2))
    return sorted(out, key=lambda r: r["T"])


def cross_check(ev, prods, sign=1.0):
    """Re-measure an event in every other product of its sector."""
    res = {}
    sec = ev["sector"]
    for key in ALT_PRODUCTS:
        if key == ev["product"] or key not in prods.get(sec, {}):
            continue
        p = prods[sec][key]
        j = np.searchsorted(p["t"], ev["T"])
        if j <= 0 or j >= len(p["t"]):
            res[key] = dict(status="no-coverage")
            continue
        u = p["seg"][min(j, len(p["t"]) - 1)]
        m = p["seg"] == u
        t, f = p["t"][m], p["f"][m] if sign > 0 else 2.0 - p["f"][m]
        r = test_free(t, f, ev["T"], ev["D"], scan=ev["D"] / 8)
        if r is None or r["cover"] < 0.5:
            res[key] = dict(status="no-coverage")
            continue
        conf = r["snr"] >= CONF_SNR and CONF_LO * ev["depth"] <= r["depth"] <= CONF_HI * ev["depth"]
        res[key] = dict(status="confirmed" if conf else "not-confirmed", depth_ppm=round(r["depth"] * 1e6, 1),
                        snr=round(r["snr"], 1))
    return res


D_FINE = np.exp(np.linspace(np.log(0.5 / 24), np.log(30 / 24), 41))


def refine(ev, segs, sign=1.0):
    """Best-SNR duration on a fine log grid (0.5-30 h) at the event, centre scanned +/- D/8.
    Used only to describe the event (classification), not for detection."""
    t, f = segs[ev["segment"]]
    f = f if sign > 0 else 2.0 - f
    best = None
    for D in D_FINE:
        r = test_free(t, f, ev["T"], D, scan=D / 8)
        if r is None or r["cover"] < DET_COVER or not r["beta_measured"]:
            continue
        if best is None or r["snr"] > best["snr"]:
            best = r
    if best is not None:
        ev.update(T_fit=best["T"], D_fit=best["D"], depth_fit=best["depth"], snr_fit=best["snr"])
    else:
        ev.update(T_fit=ev["T"], D_fit=ev["D"], depth_fit=ev["depth"], snr_fit=ev["snr"])
    return ev


def annotate(evs, prods, sign=1.0, segs=None):
    for ev in evs:
        if segs is not None:
            refine(ev, segs, sign)
        cc = cross_check(ev, prods, sign)
        ev["alt"] = cc
        conf = [k for k, v in cc.items() if v.get("status") == "confirmed"]
        tested = [k for k, v in cc.items() if v.get("status") in ("confirmed", "not-confirmed")]
        ev["n_reductions"] = 1 + len(conf)
        ev["n_alt_tested"] = len(tested)
        fam = lambda k: "SPOC" if k in SPOC_FAMILY else k  # noqa: E731
        ev["independent_confirm"] = any(fam(k) != fam(ev["product"]) for k in conf)
    return evs


# ------------------------------------------------------------------ one star
def run_star(spec, lc_dir, out_dir):
    tic = int(spec["TIC"])
    t0 = time.time()
    prods, dumps = load_products(tic, lc_dir)
    lc, used = primary(prods)
    rep = dict(spec, n_star_sectors=len(used), sectors=sorted(used), products_used=used,
               n_dump_flags=int(len(dumps)))
    if lc is None:
        rep["status"] = "no light curves"
        (out_dir / f"{tic}.json").write_text(json.dumps(rep, indent=1, default=float))
        return rep
    Dcat = spec.get("dur_h")
    durs = list(DURS_H)
    if Dcat and np.isfinite(Dcat) and all(abs(Dcat - d) / d > 0.1 for d in DURS_H):
        durs.append(float(Dcat))
    durs = sorted(d / 24 for d in durs)
    evs, why = search(lc, used, dumps, durs, +1.0)
    ctl, why_c = search(lc, used, dumps, durs, -1.0)
    segs_all = segments(lc)
    rep["events"] = annotate(evs, prods, +1.0, segs_all)
    rep["control_events"] = annotate(ctl, prods, -1.0, segs_all)
    rep["near_miss_reasons"] = why
    rep["near_miss_reasons_control"] = why_c
    # free measurement of the catalogued event at the catalogued duration
    Tev = spec.get("event_T")
    if Tev is not None and Dcat and np.isfinite(Dcat):
        segs = segments(lc)
        j = np.searchsorted(lc["t"], Tev)
        j = min(max(j, 0), len(lc["t"]) - 1)
        u = int(lc["seg"][j])
        t, f = segs[u]
        r = test_free(t, f, Tev, Dcat / 24, scan=Dcat / 24 / 4)
        if r is not None:
            ok, reason = passes(r, dumps)
            rep["catalogued_event_free"] = dict(T=r["T"], depth_ppm=r["depth"] * 1e6, snr=r["snr"], cover=r["cover"],
                                                beta=r["beta"], passes=ok, reason=reason)
    rep["runtime_s"] = round(time.time() - t0, 1)
    rep["status"] = "ok"
    (out_dir / f"{tic}.json").write_text(json.dumps(rep, indent=1, default=float))
    return rep


def run_star_safe(spec, lc_dir, out_dir):
    try:
        return run_star(spec, lc_dir, out_dir)
    except Exception:  # noqa: BLE001
        rep = dict(spec, status="error", error=traceback.format_exc()[-800:])
        (out_dir / f"{int(spec['TIC'])}.json").write_text(json.dumps(rep, indent=1, default=float))
        return rep


# ------------------------------------------------------------------ target lists
def sample_specs():
    cl = pd.read_csv(RES / "final_classification_2026-10-11_rnaas.csv")
    ct = pd.read_csv(ROOT / "singles" / "data" / "ctoi_singles.csv")
    specs = []
    for _, r in cl.iterrows():
        j = json.loads((RES / r.run / f"{int(r.TIC)}.json").read_text())
        crow = ct[ct.CTOI == r.CTOI]
        epoch = float(crow.epoch_btjd.iloc[0]) if len(crow) else float(j.get("epoch_btjd"))
        specs.append(dict(TIC=int(r.TIC), CTOI=float(r.CTOI), cls=r.cls, run=r.run, depth_ppm=float(r.depth_ppm),
                          dur_h=float(r.dur_h), epoch_btjd=epoch, event_T=float(j.get("event_T", epoch)),
                          known_repeats=[float(x["T"]) for x in j.get("repeats", [])],
                          periods=str(r.periods) if isinstance(r.periods, str) else ""))
    return specs


def eb_specs():
    """Prša+2022 EB controls: CTOI-dropped EBs fetched by fetch_ebs.py + 53847727 and 396720998
    (already in data/raw)."""
    lines = [ln for ln in open(ROOT / "data" / "raw" / "tess_eb_prsa2022.tsv") if not ln.startswith("#") and ln.strip()]
    e = pd.DataFrame([ln.rstrip("\n").split("\t") for ln in lines[3:]], columns=["TIC", "Per", "Morph", "Sectors"])
    e["TIC"] = pd.to_numeric(e.TIC, errors="coerce")
    e = e.dropna(subset=["TIC"]).astype({"TIC": int})
    e["Per"] = pd.to_numeric(e.Per, errors="coerce")
    dr = pd.read_csv(ROOT / "singles" / "data" / "ctoi_singles_dropped.csv")
    specs = []
    for tic, lc_dir in [(92349924, "eb"), (77320223, "eb"), (219771660, "eb"), (362104293, "eb"),
                        (53847727, "raw"), (396720998, "raw")]:
        P = float(e[e.TIC == tic].Per.iloc[0])
        d = dr[dr["TIC ID"] == tic]
        if len(d):
            d = d.iloc[0]
            ep = float(d["Transit Epoch (BJD)"]) - 2457000
            spec = dict(TIC=tic, prsa_P=P, depth_ppm=float(d["Depth ppm"]), dur_h=float(d["Duration (hrs)"]),
                        epoch_btjd=ep, event_T=ep, known_repeats=[], lc="eb" if lc_dir == "eb" else "raw")
        else:
            spec = dict(TIC=tic, prsa_P=P, depth_ppm=np.nan, dur_h=np.nan, epoch_btjd=None, event_T=None,
                        known_repeats=[], lc="raw")
        specs.append(spec)
    return specs


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--tics", nargs="*", type=int)
    a.add_argument("--ebs", action="store_true")
    a.add_argument("--workers", type=int, default=14)
    args = a.parse_args()
    if args.ebs:
        specs = eb_specs()
        out_dir = HERE / "data" / "eb_per_star"
    else:
        specs = sample_specs()
        out_dir = HERE / "data" / "per_star"
    if args.tics:
        specs = [s for s in specs if s["TIC"] in args.tics]
    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for s in specs:
            lc_dir = HERE / "data" / "eb_lightcurves" if s.get("lc") == "eb" else LC_DIR
            jobs.append(ex.submit(run_star_safe, s, lc_dir, out_dir))
        for k, fu in enumerate(as_completed(jobs), 1):
            r = fu.result()
            print(f"[{k}/{len(jobs)}] TIC {r['TIC']} {r.get('status')} sectors {r.get('n_star_sectors')} "
                  f"events {len(r.get('events', []))} control {len(r.get('control_events', []))} "
                  f"{r.get('runtime_s')} s {r.get('error', '')[-200:] if r.get('error') else ''}", flush=True)


if __name__ == "__main__":
    main()
