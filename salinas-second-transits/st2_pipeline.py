#!/usr/bin/env python3
"""Search the Salinas et al. (arXiv 2610.24882) TESS candidates for a second
transit in Sector 61 onward.

Steps (run in order, or `all`):
  coverage  TIC lookup, light curves on MAST and tess-point pointings for S>=61
  search    box search for a transit matching each candidate's depth and duration
  toi       match against the ExoFOP TOI and CTOI lists
  aliases   period alias set for each hit, checked against all TESS data
  report    shortlist.csv, shortlist.md and plots
  cnn       optional: rerun the Salinas CNN on each new event (needs their repo)

Example:
  python st2_pipeline.py all --candidates candidates.csv --work work
"""
import argparse
import json
import re
import sys
import threading
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

import st2core as core

warnings.filterwarnings("ignore")

FIRST_NEW_SECTOR = 61
# Pointings MAST and tess-point number outside the sector sequence. Their data
# are real, so the alias step still uses them, but they are not searched as
# new sectors.
SPECIAL_SECTORS = {1751: "3I/ATLAS special pointing, 2026 Jan 15-22 (7.2 d)"}


def is_new_sector(s):
    return s >= FIRST_NEW_SECTOR and s not in SPECIAL_SECTORS
AUTHOR_RANK = {"SPOC": 0, "TESS-SPOC": 1, "QLP": 2}
# lightkurve's download() swaps sys.stdout for devnull and back, which is not
# thread-safe (a thread can restore another's closed devnull and kill print).
# Holding this lock also stops two candidates on one star fetching the same file.
DOWNLOAD_LOCK = threading.Lock()
# wotan's numba code aborts the process when two threads enter it under numba's
# workqueue threading layer, the only one a pip install gets on macOS.
DETREND_LOCK = threading.Lock()

COLUMN_ALIASES = {
    "tic": ["tic", "tic_id", "ticid", "tic id", "id"],
    "t0": ["t0", "t0_btjd", "tmid", "t_mid", "t_mid_pred_btjd", "tc", "epoch", "t_c"],
    "depth_ppm": ["depth_ppm", "depth", "delta", "depth (ppm)"],
    "duration_h": ["duration_h", "duration", "dur_h", "t14", "t14_h", "duration (h)"],
    "sector": ["sector", "sec", "sectors"],
    "type": ["type", "class", "category", "flag", "kind"],
    "period_d": ["period_d", "period", "p", "per"],
    "depth_err": ["depth_err", "sigma_depth", "e_depth", "depth_unc"],
    "duration_err": ["duration_err", "sigma_duration", "e_duration", "duration_unc"],
}


# ---------------------------------------------------------------- inputs

def load_candidates(path):
    df = pd.read_csv(path)
    lower = {c.strip().lower(): c for c in df.columns}
    out = pd.DataFrame()
    for key, names in COLUMN_ALIASES.items():
        for n in names:
            if n in lower:
                out[key] = df[lower[n]]
                break
    for req in ("tic", "t0", "depth_ppm", "duration_h"):
        if req not in out:
            sys.exit(f"{path}: no column for '{req}'. Columns seen: {list(df.columns)}")
    out["tic"] = out["tic"].astype(str).str.extract(r"(\d+)")[0].astype(np.int64)
    out["t0"] = out["t0"].astype(float)
    out.loc[out["t0"] > 2e6, "t0"] -= 2457000.0          # BJD -> BTJD
    out.loc[out["t0"] > 2e5, "t0"] -= 2450000.0 - 7000.0  # HJD-2450000 style -> BTJD
    if out["duration_h"].median() < 1.0:
        print("WARNING: median duration < 1; are durations in days? Expecting hours.")
    if "sector" in out:
        out["sector"] = out["sector"].astype(str).str.extract(r"(\d+)")[0].astype(float)
    for c in ("type", "period_d", "depth_err", "duration_err"):
        if c not in out:
            out[c] = np.nan
    out["planet"] = df[lower["planet"]].fillna("").astype(str) if "planet" in lower else ""
    out["cand_id"] = [f"{t}{p}_{x:.3f}" for t, p, x in zip(out["tic"], out["planet"], out["t0"])]
    return out


def sector_of(mission):
    m = re.search(r"Sector\s*(\d+)", str(mission))
    return int(m.group(1)) if m else None


# ---------------------------------------------------------------- MAST

def tic_info(tic):
    from astroquery.mast import Catalogs
    r = Catalogs.query_criteria(catalog="Tic", ID=int(tic))
    if len(r) == 0:
        return {}
    r = r[0]
    def g(k):
        try:
            v = float(r[k])
            return v if np.isfinite(v) else np.nan
        except Exception:
            return np.nan
    return dict(ra=g("ra"), dec=g("dec"), tmag=g("Tmag"), teff=g("Teff"),
                rstar=g("rad"), mstar=g("mass"), contratio=g("contratio"))


def lc_products(tic):
    import lightkurve as lk
    sr = lk.search_lightcurve(f"TIC {tic}", mission="TESS")
    rows = []
    for i in range(len(sr)):
        s = sector_of(sr.mission[i])
        a = str(sr.author[i])
        if s is None or a not in AUTHOR_RANK:
            continue
        tr = (np.nan, np.nan)
        try:  # observation time span (MJD) -> BTJD, if MAST gave it
            tr = (float(sr.table["t_min"][i]) - 56999.5, float(sr.table["t_max"][i]) - 56999.5)
        except Exception:
            pass
        rows.append(dict(sector=s, author=a, exptime=float(sr.exptime[i].value), idx=i,
                         tmin=tr[0], tmax=tr[1]))
    return sr, rows


def best_per_sector(rows):
    best = {}
    for r in rows:
        key = (AUTHOR_RANK[r["author"]], r["exptime"] if r["exptime"] >= 60 else 1e9)
        if r["sector"] not in best or key < best[r["sector"]][0]:
            best[r["sector"]] = (key, r)
    return {s: v[1] for s, v in best.items()}


def pointing_sectors(tic, ra, dec):
    try:
        import tess_stars2px as tp
        res = tp.tess_stars2px_function_entry(int(tic), float(ra), float(dec))
        return sorted(set(int(s) for s in res[3] if s > 0))
    except Exception:
        return []


def load_sector_lc(sr, row, cache_dir):
    """Download one sector, return (t, f) normalised, quality-masked, NaN-free."""
    with DOWNLOAD_LOCK:
        lc = sr[row["idx"]].download(quality_bitmask="default", download_dir=str(cache_dir))
    if row["author"] == "QLP":
        # QLP sets its own flags (bits 29-30) that lightkurve's TESS mask ignores;
        # those cadences can sit at ~10% flux. QLP's advice is QUALITY == 0.
        lc = lc[lc.quality == 0]
    lc = lc.remove_nans().normalize()
    t = np.asarray(lc.time.value, float)
    f = np.asarray(lc.flux.value, float)
    m = np.isfinite(t) & np.isfinite(f)
    t, f = t[m], f[m]
    cad = core.cadence_days(t)
    if cad < 9 / 1440:  # bin anything faster than ~10 min to 10 min
        edges = np.arange(t.min(), t.max() + 10 / 1440, 10 / 1440)
        idx = np.digitize(t, edges)
        tb = np.bincount(idx, t) / np.bincount(idx)
        fb = np.bincount(idx, f) / np.bincount(idx)
        ok = np.bincount(idx) > 0
        t, f = tb[ok], fb[ok]
    return t, f


def detrend(t, f, dur_d):
    from wotan import flatten
    w = max(4.0 * dur_d, 1.0)  # wide enough to keep most of the transit depth
    with DETREND_LOCK:
        flat = flatten(t, f, method="biweight", window_length=w, edge_cutoff=0.0,
                       break_tolerance=0.5)
    m = np.isfinite(flat)
    # clip upward outliers only (flares), keep dips
    r = flat[m] - 1
    s = core.robust_std(r)
    keep = r < 5 * s
    return t[m][keep], flat[m][keep]


# ---------------------------------------------------------------- steps

def step_coverage(cands, work, workers):
    out = work / "coverage.csv"
    tics = sorted(cands["tic"].unique())
    rows = []

    def one(tic):
        info = tic_info(tic)
        _, prods = lc_products(tic)
        best = best_per_sector(prods)
        pts = pointing_sectors(tic, info.get("ra", np.nan), info.get("dec", np.nan)) if info else []
        new_lc = sorted(s for s in best if is_new_sector(s))
        new_pt = sorted(s for s in pts if is_new_sector(s))
        return dict(tic=tic, **info,
                    lc_sectors=" ".join(f"{s}:{best[s]['author']}" for s in sorted(best)),
                    lc_sectors_ge61=" ".join(map(str, new_lc)),
                    pointed_sectors_ge61=" ".join(map(str, new_pt)),
                    pointed_no_lc_ge61=" ".join(str(s) for s in new_pt if s not in best),
                    special_sectors=" ".join(str(s) for s in sorted(best) if s in SPECIAL_SECTORS),
                    n_new_lc=len(new_lc))

    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(one, t): t for t in tics}
        for k, fu in enumerate(as_completed(futs), 1):
            tic = futs[fu]
            try:
                rows.append(fu.result())
            except Exception as e:
                rows.append(dict(tic=tic, error=str(e)))
            print(f"[coverage {k}/{len(tics)}] TIC {tic}")
    cov = pd.DataFrame(rows).sort_values("tic")
    cov.to_csv(out, index=False)
    print(f"wrote {out}: {int((cov.get('n_new_lc', 0) > 0).sum())} of {len(cov)} TICs have S>={FIRST_NEW_SECTOR} light curves")
    return cov


def search_one(c, cov_row, work, args):
    """Search all S>=61 sectors of one candidate; also measure the reference event."""
    tic = int(c.tic)
    cache = work / "lc_cache"
    sr, prods = lc_products(tic)
    best = best_per_sector(prods)
    dur_d = c.duration_h / 24.0
    depth0 = c.depth_ppm * 1e-6
    res = dict(cand_id=c.cand_id, tic=tic, t0=c.t0, depth0_ppm=c.depth_ppm,
               dur0_h=c.duration_h, type=c.type, events=[])

    # reference event: sector containing t0
    ref = None
    old = [(s, r) for s, r in best.items() if s < FIRST_NEW_SECTOR]
    # sector whose time span holds t0 first, so we usually download just one
    old.sort(key=lambda sr_: 0 if (sr_[1].get("tmin", np.nan) - 1 <= c.t0 <= sr_[1].get("tmax", np.nan) + 1) else 1)
    for s, row in old:
        if np.isfinite(row.get("tmin", np.nan)) and not (row["tmin"] - 1 <= c.t0 <= row["tmax"] + 1):
            continue
        try:
            t, f = load_sector_lc(sr, row, cache)
        except Exception as e:  # fall back to the paper's depth and duration
            res["ref_error"] = f"S{s}: {e}"
            continue
        if t.min() - 0.5 <= c.t0 <= t.max() + 0.5:
            t, f = detrend(t, f, dur_d)
            ref = core.refine_event(t, f, c.t0, dur_d)
            if ref:
                ref.update(sector=s, author=row["author"])
            np.savez_compressed(work / "lc" / f"{tic}_s{s:03d}.npz", t=t, f=f)
            break
    res["ref"] = ref
    # use our own measurement of the reference depth if it is solid
    if ref and ref["snr"] >= 5 and 0.33 < ref["depth"] / depth0 < 3:
        depth_ref, dur_ref = ref["depth"], ref["dur_d"]
        derr_ref = ref["depth_err"]
    else:
        depth_ref, dur_ref = depth0, dur_d
        derr_ref = (c.depth_err * 1e-6) if np.isfinite(c.depth_err) else 0.2 * depth0
    res.update(depth_ref=depth_ref, dur_ref_d=dur_ref, depth_ref_err=derr_ref)

    for s in sorted(x for x in best if is_new_sector(x)):
        row = best[s]
        try:
            t, f = load_sector_lc(sr, row, cache)
            t, f = detrend(t, f, dur_ref)
        except Exception as e:
            res["events"].append(dict(sector=s, error=str(e)))
            continue
        np.savez_compressed(work / "lc" / f"{tic}_s{s:03d}.npz", t=t, f=f)
        cen, dep, snr, frac, sig = core.box_search(t, f, dur_ref)
        if not np.any(np.isfinite(snr)):
            continue
        anti = -np.nanmin(snr)
        for i in core.pick_events(cen, snr, dur_ref, n_max=3):
            ev = core.refine_event(t, f, cen[i], dur_ref)
            if not ev:
                continue
            ev.update(sector=s, author=row["author"], noise_ppm=sig * 1e6,
                      anti_snr=float(anti),
                      depth_ratio=ev["depth"] / depth_ref,
                      dur_ratio=ev["dur_d"] / dur_ref,
                      edge_dist_d=core.gap_edge_distance(t, ev["tc"]),
                      n_sector_points=int(t.size))
            res["events"].append(ev)
    return res


def passes(ev, args):
    if "error" in ev or ev.get("snr") is None:
        return False
    return (ev["snr"] >= args.snr_min
            and args.depth_lo <= ev["depth_ratio"] <= args.depth_hi
            and args.dur_lo <= ev["dur_ratio"] <= args.dur_hi
            and ev["snr"] >= args.anti_ratio * max(ev["anti_snr"], 1e-9))


def alias_events(passing):
    """One event per candidate for the alias step and the report.

    The highest-SNR clean passing event (core.CLEAN_EDGE_D from a gap or edge,
    core.CLEAN_FRAC of the window covered); the top-SNR event only if none is
    clean. Adds n_clean_events and alias_event_reason.
    """
    p = passing.copy()
    p["clean"] = (p["edge_dist_d"] >= core.CLEAN_EDGE_D) & (p["frac"] >= core.CLEAN_FRAC)
    top = p.sort_values("snr", ascending=False).drop_duplicates("cand_id").set_index("cand_id")
    hits = p.sort_values(["clean", "snr"], ascending=False).drop_duplicates("cand_id").copy()
    n_clean = p.groupby("cand_id")["clean"].sum()
    reasons = []
    for h in hits.itertuples():
        tp = top.loc[h.cand_id]
        q = f"edge {h.edge_dist_d:.2f} d, {h.frac:.0%} covered"
        if not h.clean:
            reasons.append(f"no clean event; top-SNR event used ({q})")
        elif abs(tp.t2 - h.t2) > 1e-6:
            reasons.append(f"highest-SNR clean event ({q}); top-SNR event S{int(tp.sector)} "
                           f"T={tp.t2:.4f} SNR {tp.snr:.1f} is not clean (edge "
                           f"{tp.edge_dist_d:.2f} d, {tp.frac:.0%} covered)")
        else:
            reasons.append(f"top-SNR event, clean ({q})")
    hits["n_clean_events"] = [int(n_clean[c]) for c in hits["cand_id"]]
    hits["alias_event_reason"] = reasons
    return hits.drop(columns="clean")


def step_search(cands, work, args):
    (work / "lc").mkdir(parents=True, exist_ok=True)
    (work / "per_cand").mkdir(parents=True, exist_ok=True)
    cov = pd.read_csv(work / "coverage.csv").set_index("tic")
    todo = []
    for c in cands.itertuples():
        if c.tic in cov.index and cov.loc[c.tic].get("n_new_lc", 0) > 0:
            p = work / "per_cand" / f"{c.cand_id}.json"
            if p.exists() and not args.force:
                continue
            todo.append(c)
    print(f"searching {len(todo)} candidates with new-sector light curves")

    def run(c):
        r = search_one(c, cov.loc[c.tic], work, args)
        (work / "per_cand" / f"{c.cand_id}.json").write_text(json.dumps(r, default=float, indent=1))
        return c.cand_id

    with ThreadPoolExecutor(args.workers) as ex:
        futs = {ex.submit(run, c): c for c in todo}
        for k, fu in enumerate(as_completed(futs), 1):
            c = futs[fu]
            try:
                fu.result()
                print(f"[search {k}/{len(todo)}] {c.cand_id}")
            except Exception as e:
                print(f"[search {k}/{len(todo)}] {c.cand_id} ERROR {e}")

    rows = []
    for p in sorted((work / "per_cand").glob("*.json")):
        r = json.loads(p.read_text())
        for ev in r["events"]:
            if ev["sector"] in SPECIAL_SECTORS:  # searched by an older version
                continue
            if "error" in ev:
                rows.append(dict(cand_id=r["cand_id"], tic=r["tic"], sector=ev["sector"], error=ev["error"]))
                continue
            rows.append(dict(cand_id=r["cand_id"], tic=r["tic"], type=r["type"], t1=r["t0"],
                             depth_ref_ppm=r["depth_ref"] * 1e6, dur_ref_h=r["dur_ref_d"] * 24,
                             sector=ev["sector"], author=ev["author"], t2=ev["tc"],
                             depth_ppm=ev["depth"] * 1e6, depth_err_ppm=ev["depth_err"] * 1e6,
                             dur_h=ev["dur_d"] * 24, snr=ev["snr"], anti_snr=ev["anti_snr"],
                             depth_ratio=ev["depth_ratio"], dur_ratio=ev["dur_ratio"],
                             edge_dist_d=ev["edge_dist_d"], frac=ev["frac"],
                             passes=passes(ev, args)))
    ev = pd.DataFrame(rows)
    ev.to_csv(work / "events_all.csv", index=False)
    n = int(ev["passes"].sum()) if "passes" in ev else 0
    print(f"wrote {work/'events_all.csv'}: {n} events pass the cuts")
    return ev


# Later-sector observations the paper itself mentions (Section 7)
PAPER_NOTES = {
    302377820: "paper already reports a transit in S68 matching candidate b",
    404736082: "paper mentions further transits in S70, association unclear",
    318706605: "paper notes two extra dips in S42 needing more sectors",
}

TOI_URLS = {
    "toi": "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv",
    "ctoi": "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv",
}
TOI_FALLBACK = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
                "select+tid,toi,tfopwg_disp,pl_orbper,pl_tranmid,pl_trandurh+from+toi&format=csv")


def _find(cols, *keys):
    for c in cols:
        lc = c.lower()
        if all(k in lc for k in keys):
            return c
    return None


def step_toi(cands, work, args):
    import requests
    lists = {}
    for name, url in TOI_URLS.items():
        p = work / f"exofop_{name}.csv"
        if not p.exists() or args.force:
            try:
                r = requests.get(url, timeout=120)
                r.raise_for_status()
                p.write_bytes(r.content)
            except Exception as e:
                print(f"could not download {name} list: {e}")
                if name == "toi":
                    r = requests.get(TOI_FALLBACK, timeout=120)
                    p.write_bytes(r.content)
        if p.exists():
            lists[name] = pd.read_csv(p, comment="#", low_memory=False)
    out = []
    for c in cands.itertuples():
        row = dict(cand_id=c.cand_id, tic=c.tic, toi="", toi_disp="", toi_period="",
                   toi_same_signal="", ctoi="", ctoi_period="")
        for name, df in lists.items():
            tcol = _find(df.columns, "tic") or _find(df.columns, "tid")
            if tcol is None:
                continue
            m = df[pd.to_numeric(df[tcol], errors="coerce") == c.tic]
            if m.empty:
                continue
            idc = _find(df.columns, name) if name == "ctoi" else (_find(df.columns, "toi") or "toi")
            pcol = _find(df.columns, "period") or _find(df.columns, "orbper")
            ecol = _find(df.columns, "epoch") or _find(df.columns, "midpoint") or _find(df.columns, "tranmid")
            ids = ";".join(str(x) for x in m[idc]) if idc in m else "yes"
            pers = ";".join(f"{float(x):.4f}" if pd.notna(x) else "nan"
                            for x in pd.to_numeric(m[pcol], errors="coerce")) if pcol else ""
            if name == "toi":
                dcol = _find(df.columns, "tfopwg") or _find(df.columns, "disp")
                row.update(toi=ids, toi_period=pers,
                           toi_disp=";".join(str(x) for x in m[dcol]) if dcol else "")
                # does a TOI ephemeris pass through the paper's event?
                same = []
                if pcol and ecol:
                    for P, E in zip(pd.to_numeric(m[pcol], errors="coerce"),
                                    pd.to_numeric(m[ecol], errors="coerce")):
                        if not (np.isfinite(P) and np.isfinite(E)) or P <= 0:
                            same.append("?")
                            continue
                        E = E - 2457000.0 if E > 2e6 else E
                        ph = ((c.t0 - E) / P + 0.5) % 1 - 0.5
                        same.append("yes" if abs(ph * P) < c.duration_h / 24 else "no")
                row["toi_same_signal"] = ";".join(same)
            else:
                row.update(ctoi=ids, ctoi_period=pers)
        out.append(row)
    df = pd.DataFrame(out)
    df.to_csv(work / "toi_match.csv", index=False)
    print(f"wrote {work/'toi_match.csv'}: {int((df['toi'] != '').sum())} on the TOI list, "
          f"{int((df['ctoi'] != '').sum())} on the CTOI list")
    return df


def all_lc_for_tic(tic, work, dur_d):
    """Every sector (old and new) for alias checks, detrended, cached as npz."""
    sr, prods = lc_products(tic)
    best = best_per_sector(prods)
    T, F, spans = [], [], []
    for s, row in sorted(best.items()):
        p = work / "lc" / f"{tic}_s{s:03d}.npz"
        if p.exists():
            d = np.load(p)
            t, f = d["t"], d["f"]
        else:
            try:
                t, f = load_sector_lc(sr, row, work / "lc_cache")
                t, f = detrend(t, f, dur_d)
                np.savez_compressed(p, t=t, f=f)
            except Exception as e:
                print(f"  TIC {tic} S{s}: {e}")
                continue
        T.append(t)
        F.append(f)
        spans.append((s, row["author"], float(t.min()), float(t.max())))
    return np.concatenate(T), np.concatenate(F), sorted(best), spans


def load_hypotheses(path):
    """{tic: [quoted periods]} from period_hypotheses.csv, if it exists."""
    if not path or not Path(path).exists():
        return {}
    out = {}
    for r in pd.read_csv(path, comment="#").itertuples():
        out.setdefault(int(r.tic), []).append(float(r.period_d))
    return out


def epoch_rows(h, P, source, t, f, spans, dur_d, depth, depth_err, sig, edges):
    """One row per predicted epoch of P (through the paper event) that has any data."""
    k = np.arange(int(np.floor((t.min() - h.t1) / P)) - 1, int(np.ceil((t.max() - h.t1) / P)) + 2)
    ep = h.t1 + k * P
    v = core.epoch_verdicts(ep, t, f, dur_d, depth, depth_err, sig, edges)
    rows = []
    for i in np.where(v["frac"] > 0)[0]:
        sec = [(s, a) for s, a, lo, hi in spans if lo - 0.5 <= ep[i] <= hi + 0.5]
        verdict, old = v["verdict"][i], v["verdict_any_coverage"][i]
        if abs(ep[i] - h.t1) < dur_d / 2:
            verdict = old = "paper event"
        elif abs(ep[i] - h.t2) < dur_d / 2:
            verdict = old = "alias event (T2)"
        rows.append(dict(cand_id=h.cand_id, period_d=P, period_source=source, k=int(k[i]),
                         epoch_btjd=ep[i], sector=sec[0][0] if sec else np.nan,
                         pipeline=sec[0][1] if sec else "", coverage=v["frac"][i],
                         edge_dist_d=v["edge_d"][i], box_depth_ppm=v["depth"][i] * 1e6,
                         box_depth_err_ppm=v["depth_err"][i] * 1e6, verdict=verdict,
                         verdict_any_coverage=old))
    return rows


def epoch_status(rows):
    """Alias status from epoch rows, under the clean-epoch rule and the old any-coverage rule."""
    def st(col):
        v = {r[col] for r in rows}
        return "excluded" if "excludes" in v else "supported" if "supports" in v else "allowed"
    return st("verdict"), st("verdict_any_coverage")


def plot_epochs(h, P, rows, t, f, outdir):
    """Paper event next to every epoch that excludes P under either rule."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bad = [r for r in rows if "excludes" in (r["verdict"], r["verdict_any_coverage"])]
    if not bad:
        return ""
    show = [None] + bad[:11]
    nc = min(len(show), 4)
    nr = int(np.ceil(len(show) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(4 * nc, 3 * nr), squeeze=False, sharey=True)
    dur = h.dur_ref_h / 24
    w = max(3 * dur, 0.75)
    for a, r in zip(axs.flat, show):
        tc = h.t1 if r is None else r["epoch_btjd"]
        m = np.abs(t - tc) < w
        x, y = (t[m] - tc) * 24, (f[m] - 1) * 1e6
        a.plot(x, y, ".", ms=2, color="0.6")
        if m.sum():
            edges = np.arange(-w, w + 1 / 48, 1 / 48) * 24
            idx = np.digitize(x, edges)
            xb = [x[idx == j].mean() for j in np.unique(idx)]
            yb = [y[idx == j].mean() for j in np.unique(idx)]
            a.plot(xb, yb, "o", ms=3, color="C0")
        a.axhline(-h.depth_ref_ppm, color="C3", lw=0.8, ls="--")
        a.axvspan(-dur * 12, dur * 12, color="C3", alpha=0.08)
        if r is None:
            a.set_title(f"paper event T={h.t1:.3f}", fontsize=8)
        else:
            a.set_title(f"S{r['sector']:.0f} {r['pipeline']} T={tc:.3f}\ncov {r['coverage']:.0%}, edge "
                        f"{r['edge_dist_d']:.2f} d, {r['box_depth_ppm']:.0f}+/-{r['box_depth_err_ppm']:.0f} ppm"
                        f"\nnow: {r['verdict']}", fontsize=8)
        a.set_xlabel("hours")
    for a in list(axs.flat)[len(show):]:
        a.axis("off")
    fig.suptitle(f"TIC {h.tic} P = {P:.5f} d: epochs that exclude it ({len(bad)}"
                 + (", first 11 shown" if len(bad) > 11 else "") + ")", fontsize=10)
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    name = f"TIC{h.tic}_{h.cand_id}_P{P:.4f}.png"
    fig.savefig(outdir / name, dpi=100)
    plt.close(fig)
    return name


def distinct_events(g, merge_d=0.5):
    """Passing event times with events within merge_d of a higher-SNR one dropped."""
    kept = []
    for t2 in g.sort_values("snr", ascending=False)["t2"]:
        if all(abs(t2 - k) > merge_d for k in kept):
            kept.append(t2)
    return sorted(kept)


def common_periods(passing, work, p_floor):
    """Best common period for every candidate with 2+ distinct passing events.

    flag: 'periodic' when k >= 4, k/N >= 0.8 and P >= 5 d; 'chance' when the
    best P is under 3 d (dips that only line up at short periods). One stray
    dip can let a short chance period beat the true one, so the best P >= 5 d
    is recorded too (cp_long_*).
    """
    rows = []
    for cid, g in passing.groupby("cand_id"):
        times = distinct_events(g)
        if len(times) < 2:
            continue
        h = g.iloc[0]
        cp = core.common_period(h.t1, times, h.dur_ref_h / 24, p_floor)
        lp = core.common_period(h.t1, times, h.dur_ref_h / 24, max(p_floor, 5.0))
        flag = ""
        if np.isfinite(cp["P"]):
            if cp["k"] >= 4 and cp["k"] / cp["N"] >= 0.8 and cp["P"] >= 5:
                flag = "periodic"
            elif cp["P"] < 3:
                flag = "chance"
        rows.append(dict(cand_id=cid, tic=h.tic, cp_n_distinct_events=len(times), cp_N=cp["N"],
                         cp_k=cp["k"], cp_period_d=cp["P"], cp_flag=flag,
                         cp_matched=";".join(f"{x:.3f}" for x in cp.get("matched", [])),
                         cp_long_period_d=lp["P"], cp_long_k=lp["k"]))
    df = pd.DataFrame(rows, columns=["cand_id", "tic", "cp_n_distinct_events", "cp_N", "cp_k",
                                     "cp_period_d", "cp_flag", "cp_matched", "cp_long_period_d", "cp_long_k"])
    df.to_csv(work / "common_period.csv", index=False)
    print(f"wrote {work/'common_period.csv'}: {len(df)} candidates with 2+ distinct events, "
          f"{int((df['cp_flag'] == 'periodic').sum())} periodic, {int((df['cp_flag'] == 'chance').sum())} chance")
    return df


def step_aliases(cands, work, args):
    ev = pd.read_csv(work / "events_all.csv")
    passing = ev[ev["passes"] == True]
    hits = alias_events(passing)
    hyp = load_hypotheses(getattr(args, "hypotheses", None))
    cps = common_periods(passing, work, args.p_floor).set_index("cand_id")
    hyp_rows = []
    cov = pd.read_csv(work / "coverage.csv").set_index("tic")
    (work / "aliases").mkdir(exist_ok=True)
    summ = []
    for h in hits.itertuples():
        info = cov.loc[h.tic] if h.tic in cov.index else {}
        mstar = float(info.get("mstar", np.nan)) if len(info) else np.nan
        rstar = float(info.get("rstar", np.nan)) if len(info) else np.nan
        dur_d = h.dur_ref_h / 24
        depth = h.depth_ref_ppm * 1e-6
        t, f, secs, spans = all_lc_for_tic(h.tic, work, dur_d)
        pmin = core.pmin_from_duration(dur_d, mstar, rstar, depth)
        derr = h.depth_err_ppm * 1e-6
        rows = core.alias_table(h.t1, h.t2, dur_d, depth, derr, t, f,
                                p_floor=args.p_floor, pmin_density=pmin)
        at = pd.DataFrame(rows)
        at.to_csv(work / "aliases" / f"{h.cand_id}.csv", index=False)
        ok = at[at["status"].isin(["allowed", "supported"])]

        # per-epoch table: period hypotheses, then every alias still allowed or supported
        sig, edges = core.duration_scale_noise(t, f, dur_d), core.gap_edges(t)
        times = passing.loc[passing["cand_id"] == h.cand_id, "t2"].to_numpy()
        periods = []
        for P0 in hyp.get(h.tic, []):
            P, k, N = core.refine_period(h.t1, times, P0, dur_d)
            periods.append((P, f"hypothesis {P0:g} d ({k}/{N} events match)", P0))
        if h.cand_id in cps.index and np.isfinite(cps.loc[h.cand_id, "cp_period_d"]):
            c = cps.loc[h.cand_id]
            Pc, kc = (c.cp_period_d, c.cp_k) if c.cp_period_d >= 5 else (c.cp_long_period_d, c.cp_long_k)
            if np.isfinite(Pc) and all(abs(Pc - q[0]) > 1e-4 for q in periods):
                periods.append((Pc, f"common period ({kc}/{c.cp_N} events match)", Pc))
        periods += [(r.period_d, f"alias n={r.n} ({r.status})", None)
                    for r in ok.itertuples()]
        erows = []
        for P, source, P0 in periods:
            er = epoch_rows(h, P, source, t, f, spans, dur_d, depth, derr, sig, edges)
            erows += er
            if P0 is None:
                continue
            now, old = epoch_status(er)
            hyp_rows.append(dict(
                cand_id=h.cand_id, tic=h.tic, period_quoted_d=P0, period_used_d=P, source=source,
                status_old_rule=old, status_now=now,
                n_epochs_with_data=sum(r["verdict"] not in ("no data", "paper event", "alias event (T2)") for r in er),
                n_supporting=sum(r["verdict"] == "supports" for r in er),
                n_inconclusive=sum(r["verdict"] == "inconclusive" for r in er),
                excluding_old_rule=";".join(f"S{r['sector']:.0f}:{r['epoch_btjd']:.3f}" for r in er
                                            if r["verdict_any_coverage"] == "excludes"),
                excluding_now=";".join(f"S{r['sector']:.0f}:{r['epoch_btjd']:.3f}" for r in er
                                       if r["verdict"] == "excludes"),
                plot=plot_epochs(h, P, er, t, f, Path(args.out) / "epoch_plots")))
        pd.DataFrame(erows).to_csv(work / "aliases" / f"{h.cand_id}_epochs.csv", index=False)
        summ.append(dict(cand_id=h.cand_id, tic=h.tic, alias_event_t2=h.t2,
                         alias_event_sector=int(h.sector), alias_event_author=h.author,
                         alias_event_snr=h.snr, alias_event_reason=h.alias_event_reason,
                         n_clean_events=h.n_clean_events, dT_d=h.t2 - h.t1,
                         pmin_circ_d=pmin, n_aliases=len(at),
                         n_allowed=int((at["status"] == "allowed").sum()),
                         n_supported=int((at["status"] == "supported").sum()),
                         n_excluded=int((at["status"] == "excluded").sum()),
                         n_below_pmin=int((at["status"] == "allowed_below_pmin").sum()),
                         longest_allowed_d=float(ok["period_d"].max()) if len(ok) else np.nan,
                         shortest_allowed_d=float(ok["period_d"].min()) if len(ok) else np.nan,
                         supported_periods=";".join(f"{p:.4f}" for p in
                                                    at.loc[at["status"] == "supported", "period_d"]),
                         sectors_all=" ".join(map(str, secs))))
        print(f"[aliases] {h.cand_id}: dT={h.t2-h.t1:.2f} d, {summ[-1]['n_allowed']} allowed, "
              f"{summ[-1]['n_supported']} supported")
    df = pd.DataFrame(summ)
    df.to_csv(work / "aliases_summary.csv", index=False)
    pd.DataFrame(hyp_rows).to_csv(work / "aliases" / "hypotheses_summary.csv", index=False)
    return df


def plot_hit(h, work, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def load(s):
        p = work / "lc" / f"{h.tic}_s{int(s):03d}.npz"
        if not p.exists():
            return None, None
        d = np.load(p)
        return d["t"], d["f"]

    r = json.loads((work / "per_cand" / f"{h.cand_id}.json").read_text())
    ref_sector = (r.get("ref") or {}).get("sector")
    dur = h.dur_ref_h / 24
    fig, ax = plt.subplots(3, 1, figsize=(10, 9))
    for a, s, tc, lab in ((ax[0], ref_sector, h.t1, "paper event"),
                          (ax[1], h.sector, h.t2, "new event")):
        if s is None:
            a.set_title(f"{lab}: sector not found")
            continue
        t, f = load(s)
        if t is None:
            continue
        m = np.abs(t - tc) < max(4 * dur, 1.0)
        a.plot((t[m] - tc) * 24, (f[m] - 1) * 1e6, ".", ms=2, color="0.6")
        # 30-min bins
        tb = np.arange(-max(4 * dur, 1.0), max(4 * dur, 1.0), 1 / 48)
        idx = np.digitize(t[m] - tc, tb)
        fb = [np.mean(f[m][idx == k]) for k in np.unique(idx)]
        xb = [np.mean(t[m][idx == k] - tc) for k in np.unique(idx)]
        a.plot(np.array(xb) * 24, (np.array(fb) - 1) * 1e6, "o", ms=4, color="C0")
        a.axhline(-h.depth_ref_ppm, color="C3", lw=0.8, ls="--")
        a.axvspan(-dur * 12, dur * 12, color="C3", alpha=0.08)
        a.set_title(f"TIC {h.tic}  S{int(s)}  {lab}  T={tc:.4f} BTJD")
        a.set_xlabel("hours from mid-transit")
        a.set_ylabel("ppm")
    t, f = load(h.sector)
    if t is not None:
        ax[2].plot(t, (f - 1) * 1e6, ".", ms=1, color="0.5")
        ax[2].axvline(h.t2, color="C3", alpha=0.5)
        ax[2].set_title(f"S{int(h.sector)} full sector ({h.author}); SNR {h.snr:.1f}, "
                        f"depth ratio {h.depth_ratio:.2f}, duration ratio {h.dur_ratio:.2f}")
        ax[2].set_xlabel("BTJD")
        ax[2].set_ylabel("ppm")
    fig.tight_layout()
    p = outdir / f"TIC{h.tic}_{h.t2:.2f}.png"
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return p.name


def step_report(cands, work, args):
    out = Path(args.out)
    (out / "plots").mkdir(parents=True, exist_ok=True)
    ev = pd.read_csv(work / "events_all.csv")
    passing = ev[ev["passes"] == True].sort_values("snr", ascending=False)
    hits = alias_events(passing)
    # the other passing events, so no dip is hidden
    chosen = dict(zip(hits["cand_id"], hits["t2"]))
    others = {cid: g[(g["t2"] - chosen[cid]).abs() > 1e-6] for cid, g in passing.groupby("cand_id")}
    hits["n_passing_events"] = [len(others[c]) + 1 for c in hits["cand_id"]]
    hits["other_passing"] = [";".join(f"S{int(o.sector)}:{o.t2:.4f}:snr{o.snr:.1f}"
                                      for o in others[c].itertuples()) for c in hits["cand_id"]]
    toi = pd.read_csv(work / "toi_match.csv", dtype=str).fillna("") if (work / "toi_match.csv").exists() else None
    al = pd.read_csv(work / "aliases_summary.csv") if (work / "aliases_summary.csv").exists() else None
    cp = pd.read_csv(work / "common_period.csv") if (work / "common_period.csv").exists() else None
    if toi is not None:
        hits = hits.merge(toi.drop(columns=["tic"]), on="cand_id", how="left")
    if al is not None and len(al):
        hits = hits.merge(al.drop(columns=["tic", "n_clean_events", "alias_event_reason"], errors="ignore"),
                          on="cand_id", how="left")
    if cp is not None and len(cp):
        hits = hits.merge(cp.drop(columns=["tic"]), on="cand_id", how="left")
    paper = cands.set_index("cand_id")[["depth_ppm", "duration_h"]].rename(
        columns={"depth_ppm": "paper_depth_ppm", "duration_h": "paper_dur_h"})
    hits = hits.merge(paper, left_on="cand_id", right_index=True, how="left")
    hits["plot"] = [plot_hit(h, work, out / "plots") for h in hits.itertuples()]
    hits.to_csv(out / "shortlist.csv", index=False)

    cov = pd.read_csv(work / "coverage.csv")
    # tess-point lists scheduled sectors too; count only those up to the latest
    # sector that has any light curve on MAST, i.e. already observed
    lc_secs = [int(x.split(":")[0]) for v in cov["lc_sectors"].dropna() for x in str(v).split()]
    last_sec = max([x for x in lc_secs if x < 1000], default=0)
    missing = cov.get("pointed_no_lc_ge61", pd.Series(dtype=str)).fillna("").astype(str)
    n_missing = int(missing.apply(lambda v: any(int(float(x)) <= last_sec for x in v.split())).sum())
    lines = ["# Second-transit shortlist", "",
             f"Candidates in: {len(cands)}. TICs with S>={FIRST_NEW_SECTOR} light curves: "
             f"{int((cov.get('n_new_lc', 0) > 0).sum())} of {len(cov)}. "
             f"Events passing cuts: {len(hits)}. TICs that tess-point says were observed in "
             f"S{FIRST_NEW_SECTOR}-S{last_sec} with no light curve on MAST yet: {n_missing} "
             f"(see pointed_no_lc_ge61 in coverage.csv, which also lists scheduled sectors; "
             f"these need FFI extraction, e.g. TESSCut).", "",
             f"Cuts: SNR >= {args.snr_min}, depth ratio {args.depth_lo}-{args.depth_hi}, "
             f"duration ratio {args.dur_lo}-{args.dur_hi}, SNR >= {args.anti_ratio} x strongest brightening.",
             ""]
    for h in hits.itertuples():
        notes = []
        if getattr(h, "toi", ""):
            notes.append(f"already TOI {h.toi} (P={h.toi_period}, same signal: {h.toi_same_signal}, disp {h.toi_disp})")
        if getattr(h, "ctoi", ""):
            notes.append(f"already CTOI {h.ctoi} (P={h.ctoi_period})")
        if h.tic in PAPER_NOTES:
            notes.append(PAPER_NOTES[h.tic])
        notes.append(f"alias event: {h.alias_event_reason}")
        if pd.notna(getattr(h, "cp_period_d", np.nan)):
            notes.append(f"common period {h.cp_period_d:.4f} d matches {int(h.cp_k)} of {int(h.cp_N)} "
                         f"events incl. the paper's" + (f" [{h.cp_flag}]" if h.cp_flag == h.cp_flag and h.cp_flag else "")
                         + (f"; best P >= 5 d: {h.cp_long_period_d:.4f} d, {int(h.cp_long_k)} of {int(h.cp_N)}"
                            if pd.notna(h.cp_long_period_d) and abs(h.cp_long_period_d - h.cp_period_d) > 1e-4 else ""))
        if h.edge_dist_d < 0.5:
            notes.append(f"event {h.edge_dist_d:.2f} d from a gap or sector edge, check systematics")
        if h.frac < 0.8:
            notes.append(f"only {h.frac:.0%} of the transit window has data")
        for o in others[h.cand_id].itertuples():
            flags = [x for x, bad in (("near edge", o.edge_dist_d < 0.5),
                                      (f"{o.frac:.0%} covered", o.frac < 0.8)) if bad]
            notes.append(f"also passes: S{int(o.sector)} T={o.t2:.4f} SNR {o.snr:.1f}, depth ratio "
                         f"{o.depth_ratio:.2f}, duration ratio {o.dur_ratio:.2f}"
                         + (f" ({', '.join(flags)})" if flags else ""))
        a = ""
        if al is not None and pd.notna(getattr(h, "n_allowed", np.nan)):
            a = (f" dT = {h.dT_d:.2f} d, so P <= {h.dT_d:.2f} d. Aliases: {int(h.n_allowed)} allowed, "
                 f"{int(h.n_excluded)} excluded by data, {int(h.n_below_pmin)} below the circular-orbit "
                 f"minimum ({h.pmin_circ_d:.1f} d). Longest allowed {h.longest_allowed_d:.2f} d.")
            if isinstance(h.supported_periods, str) and h.supported_periods:
                a += f" Possible third transit at P = {h.supported_periods} d."
        lines += [f"## TIC {h.tic} (paper T1 {h.t1:.4f}, new T2 {h.t2:.4f}, S{int(h.sector)} {h.author})",
                  f"Depth {h.depth_ppm:.0f} +/- {h.depth_err_ppm:.0f} ppm vs {h.depth_ref_ppm:.0f} ppm; "
                  f"duration {h.dur_h:.2f} h vs {h.dur_ref_h:.2f} h (box fits; paper CNN: "
                  f"{h.paper_depth_ppm:.0f} ppm, {h.paper_dur_h:.2f} h); SNR {h.snr:.1f}.{a}",
                  ("Notes: " + "; ".join(notes)) if notes else "Notes: none flagged.",
                  f"![](plots/{h.plot})", ""]
    (out / "shortlist.md").write_text("\n".join(lines))
    print(f"wrote {out/'shortlist.csv'} and {out/'shortlist.md'} ({len(hits)} hits)")


def fresh_lc(lc):
    """The same light curve rebuilt from writable copies of time, flux, flux_err and quality.

    QLP files from S94 on carry a read-only column that makes lightkurve's
    bin() (inside the CNN's preprocess_lc) fail with "output array is read-only".
    """
    import lightkurve as lk
    new = lk.LightCurve(time=lc.time.copy(), flux=np.array(lc.flux.value, copy=True) * lc.flux.unit,
                        flux_err=np.array(lc.flux_err.value, copy=True) * lc.flux_err.unit)
    q = getattr(lc, "quality", None)
    if q is not None:
        new["quality"] = np.array(getattr(q, "value", q), copy=True)
    return new


def cnn_targets(work, extra_path):
    """Every passing event, the paper event of each such candidate, and extra transits."""
    ev = pd.read_csv(work / "events_all.csv")
    passing = ev[ev["passes"] == True]
    rows = [dict(cand_id=h.cand_id, tic=h.tic, sector=int(h.sector), author=h.author, tc=h.t2, which="new")
            for h in passing.itertuples()]
    for cid, g in passing.groupby("cand_id"):
        ref = json.loads((work / "per_cand" / f"{cid}.json").read_text()).get("ref") or {}
        if ref.get("sector") is not None:
            rows.append(dict(cand_id=cid, tic=int(g.tic.iloc[0]), sector=int(ref["sector"]),
                             author=ref["author"], tc=float(g.t1.iloc[0]), which="paper"))
    if extra_path and Path(extra_path).exists():
        rows += pd.read_csv(Path(extra_path).resolve(), comment="#").to_dict("records")
    return pd.DataFrame(rows).drop_duplicates(["tic", "sector", "author", "tc"])


def step_cnn(cands, work, args):
    """Refit events with the Salinas CNN (their repo, any pipeline)."""
    repo = Path(args.cnn_repo).resolve()
    work = work.resolve()  # before chdir, or a relative --work points into the CNN repo
    todo = cnn_targets(work, getattr(args, "cnn_extra", None))
    sys.path.insert(0, str(repo))
    import os
    os.chdir(repo)  # their config uses repo-relative paths
    import lightkurve as lk
    from computeparams.param_inference_core import pick_lightcurve, preprocess_lc, _run_inference_on_lc
    rows = []
    for k, h in enumerate(todo.itertuples(), 1):
        base = dict(cand_id=h.cand_id, tic=h.tic, sector=h.sector, author=h.author, which=h.which, t_input=h.tc)
        try:
            sr = lk.search_lightcurve(f"TIC {h.tic}", sector=int(h.sector), author=h.author)
            if len(sr) == 0:
                rows.append(dict(base, error="no light curve on MAST"))
                continue
            lc = preprocess_lc(fresh_lc(pick_lightcurve(sr[0].download(download_dir=str(work / "cnn_cache")))))
            r = _run_inference_on_lc(lc, h.tc, plot=False)
            rows.append(dict(base, **{k_: r[k_] for k_ in ("depth_ppm", "sigma_depth", "duration_h",
                                                           "sigma_duration", "t_mid_pred_btjd")}))
        except Exception as e:
            rows.append(dict(base, error=str(e)))
        print(f"[cnn {k}/{len(todo)}] {h.cand_id} S{h.sector} {h.which}: {rows[-1].get('error', 'ok')}")
    df = pd.DataFrame(rows)
    df.to_csv(work / "cnn_new_events.csv", index=False)
    n_err = int(df["error"].notna().sum()) if "error" in df else 0
    print(f"wrote cnn_new_events.csv: {len(df) - n_err} fitted, {n_err} failed")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["coverage", "search", "toi", "aliases", "report", "cnn", "all"])
    ap.add_argument("--candidates", default="candidates.csv")
    ap.add_argument("--work", default="work")
    ap.add_argument("--out", default="results")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--snr-min", type=float, default=7.0)
    ap.add_argument("--depth-lo", type=float, default=0.5)
    ap.add_argument("--depth-hi", type=float, default=2.0)
    ap.add_argument("--dur-lo", type=float, default=0.5)
    ap.add_argument("--dur-hi", type=float, default=2.0)
    ap.add_argument("--anti-ratio", type=float, default=1.5,
                    help="event SNR must beat the strongest brightening by this factor")
    ap.add_argument("--p-floor", type=float, default=1.0, help="shortest alias period to list (d)")
    ap.add_argument("--hypotheses", default="period_hypotheses.csv",
                    help="tic,period_d list tabulated epoch by epoch in the alias step")
    ap.add_argument("--types", default="single,multiplanet",
                    help="comma list of candidate types to search (single, multi, multiplanet, all)")
    ap.add_argument("--cnn-repo", default="singletransiters")
    ap.add_argument("--cnn-extra", default="cnn_extra.csv",
                    help="extra transits for the CNN step: cand_id,tic,sector,author,tc,which")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    cands = load_candidates(args.candidates)
    if args.types != "all" and cands["type"].notna().any():
        keep = [x.strip().lower() for x in args.types.split(",")]
        cands = cands[cands["type"].astype(str).str.lower().isin(keep)]
    print(f"{len(cands)} candidates, {cands['tic'].nunique()} TICs")
    steps = ["coverage", "search", "toi", "aliases", "report"] if args.step == "all" else [args.step]
    for s in steps:
        print(f"\n=== {s} ===")
        if s == "coverage":
            step_coverage(cands, work, args.workers)
        else:
            globals()[f"step_{s}"](cands, work, args)


if __name__ == "__main__":
    main()
