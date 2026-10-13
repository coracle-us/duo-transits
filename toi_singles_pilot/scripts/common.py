"""
Shared code for the TOI-singles pilot (toi_singles_pilot/). Reuses the validated k2singles v3 search
(k2singles/scripts/k2lib.py: grid / find_repeats / robust_periods) on TESS-only light curves.

All downloads are redirected into toi_singles_pilot/data/ (gitignored):
  data/lightcurves/<TIC>/s<NNN>_<AUTHOR>_<EXP>s.csv.gz   pipeline light curves + own TESScut photometry
  data/tesscut/<TIC>/                                   TESScut FFI cutouts (own photometry, difference images)
  data/lk_cache/                                        lightkurve FITS download cache
Nothing under data/raw or other project folders is written.

Light curves: ONE product per sector is downloaded, the one alias_pipeline.load_lightcurves would pick
(SPOC 120 s > TESS-SPOC > QLP > SPOC 20 s), which avoids pulling GBs of 20-s and TGLC files. A sector
that TESScut lists but that has none of these products gets our own aperture photometry
(the same algorithm as singles/scripts/tesscut_photometry.py, writing into data/lightcurves).

Per TOI (run_toi):
  1. light curves; transits of OTHER TOIs on the same star that have a real period (TOI table,
     data/raw/exofop_toi.csv; not single/duo/max-period rows) are masked (+-0.75 duration + 1 h)
  2. catalogued event check (alias_pipeline.test_window at the TOI epoch, +-0.01 d scan)
  3. k2lib.grid (box test of the TOI depth and duration every D/4), k2lib.find_repeats (SNR >= 8,
     depth 0.5-2x, coverage >= 0.8, baseline both sides), k2lib.robust_periods (P >= 5 d)
  4. period constraints on a FINE period grid (dP small enough that every predicted epoch moves by
     <= D/8 between grid points, so the allowed set is not a grid artefact; review point M3):
     frac_logP_excluded, minimum allowed P, number of allowed ranges
  5. windows plot; optional difference images at the event and repeats (pilot.compute_diffs)
"""
import json
import re
import sys
import time
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
PILOT = HERE.parent
ROOT = PILOT.parent
DATA = PILOT / "data"
LCD = DATA / "lightcurves"
sys.path.insert(0, str(ROOT / "k2singles" / "scripts"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
import k2lib  # noqa: E402
import alias_pipeline as ap  # noqa: E402
import pilot  # noqa: E402
import diff_image as di  # noqa: E402

# redirect every cache/output path of the reused modules into toi_singles_pilot/data
ap.LC = LCD
di.RAW = DATA / "tesscut"
k2lib.CACHE = DATA / "k2cache"
PRODUCTS = ap.PRODUCT_RANK + ["OWN"]
AUTHORS = ("SPOC", "TESS-SPOC", "QLP")
EXTRA_COLS = ("sap_flux", "pdcsap_flux", "kspsap_flux", "sap_bkg", "mom_centr1", "mom_centr2", "pos_corr1",
              "pos_corr2", "orbitid")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _col(x):
    v = x.value if hasattr(x, "value") else x
    mask = getattr(v, "mask", None)
    arr = np.asarray(getattr(v, "unmasked", v), dtype=float)
    if mask is not None:
        arr = np.where(np.asarray(mask), np.nan, arr)
    return arr


def _rank(author, exp):
    key = f"SPOC_{int(round(exp))}" if author == "SPOC" else author
    return ap.PRODUCT_RANK.index(key) if key in ap.PRODUCT_RANK else 99


def fetch_lightcurves(tic, extra_authors=(), only_sectors=None):
    """Best product per sector (plus every product of `extra_authors`, for vetting). Resumable."""
    import lightkurve as lk
    tdir = LCD / str(tic)
    tdir.mkdir(parents=True, exist_ok=True)
    sr = None
    for attempt in range(4):
        try:
            sr = lk.search_lightcurve(f"TIC {tic}")
            break
        except Exception as e:  # noqa: BLE001
            log(f"  TIC {tic} search retry {attempt + 1}: {e}")
            time.sleep(10 * (attempt + 1))
    if sr is None:
        raise RuntimeError("lightkurve search failed")
    tab = sr.table.to_pandas()
    tab["idx"] = np.arange(len(tab))
    tab["sector"] = pd.to_numeric(tab.mission.astype(str).str.extract(r"(\d+)")[0], errors="coerce")
    tab = tab[tab.author.isin(AUTHORS) & tab.sector.notna() & (tab.sector <= 200)]
    if "target_name" in tab:
        tn = tab.target_name.astype(str).str.lstrip("0")
        tab = tab[tn == str(tic)]
    tab["rank"] = [_rank(a, e) for a, e in zip(tab.author, tab.exptime)]
    tab[["mission", "author", "exptime", "target_name", "sector", "rank"]].to_csv(tdir / "search_results.csv",
                                                                                   index=False)
    pick = tab.sort_values("rank").groupby("sector").head(1)
    if extra_authors:
        ex = tab[tab.author.isin(extra_authors)]
        if only_sectors is not None:
            ex = ex[ex.sector.isin(list(only_sectors))]
        pick = pd.concat([pick, ex]).drop_duplicates("idx")
    n_fail = 0
    for r in pick.itertuples():
        fname = tdir / f"s{int(r.sector):03d}_{r.author}_{int(round(r.exptime))}s.csv.gz"
        if fname.exists():
            continue
        lc = None
        for attempt in range(3):
            try:
                lc = sr[int(r.idx)].download(quality_bitmask="none", download_dir=str(DATA / "lk_cache"))
                break
            except Exception as e:  # noqa: BLE001
                log(f"  TIC {tic} download retry {attempt + 1} ({r.author} s{int(r.sector)}): {str(e)[:120]}")
                m = re.search(r"Data product (.+?\.fits)", str(e))
                if m and Path(m.group(1)).exists():
                    Path(m.group(1)).unlink()
                time.sleep(10)
        if lc is None:
            n_fail += 1
            continue
        df = pd.DataFrame({"time": _col(lc.time), "flux": _col(lc.flux), "flux_err": _col(lc.flux_err),
                           "quality": np.nan_to_num(_col(lc.quality)).astype(int)})
        for extra in EXTRA_COLS:
            if extra in lc.columns:
                df[extra] = _col(lc[extra])
        tmp = fname.with_suffix(".tmp")
        df.to_csv(tmp, index=False, compression="gzip")
        tmp.replace(fname)
    return dict(n_products=len(pick), n_failed=n_fail)


def own_photometry(tic, sector, ra, dec):
    """singles/scripts/tesscut_photometry.photometry, writing into data/lightcurves."""
    from astropy.io import fits
    from astropy.wcs import WCS
    path = di.get_cutout(tic, ra, dec, sector)
    if path is None:
        raise RuntimeError(f"no TESScut cutout for S{sector}")
    with fits.open(path) as h:
        tab, wcs = h[1].data, WCS(h[2].header)
    t, flux, q = tab["TIME"], tab["FLUX"].astype(float), tab["QUALITY"]
    tx, ty = wcs.all_world2pix(ra, dec, 0)
    yy, xx = np.mgrid[0:flux.shape[1], 0:flux.shape[2]]
    med = np.nanmedian(flux, axis=0)
    bkg_pix = med <= np.nanpercentile(med, 40)
    bkg = np.nanmedian(flux[:, bkg_pix], axis=1)
    apm = (np.hypot(xx - tx, yy - ty) <= 1.6) & (med > np.nanpercentile(med, 60))
    f = np.nansum(flux[:, apm] - bkg[:, None], axis=1)
    good = np.isfinite(t) & np.isfinite(f) & (q == 0) & (f > 0)
    if good.sum() < 100:
        raise RuntimeError(f"S{sector}: only {good.sum()} good cadences")
    cad = int(round(np.nanmedian(np.diff(t[good])) * 86400))
    fe = np.sqrt(np.abs(f)) / np.nanmedian(f[good])
    df = pd.DataFrame({"time": t[good], "flux": f[good], "flux_err": fe[good] * np.nanmedian(f[good]), "quality": 0})
    out = LCD / str(tic) / f"s{sector:03d}_OWN_{cad}s.csv.gz"
    df.to_csv(out, index=False, compression="gzip")
    return out


def add_own_photometry(tic, ra, dec, min_sector=1):
    import lightkurve as lk
    have = set()
    for f in (LCD / str(tic)).glob("s*.csv.gz"):
        m = re.fullmatch(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", f.name)
        if m and m.group(2) in ("SPOC", "TESS-SPOC", "QLP", "OWN"):
            have.add(int(m.group(1)))
    secs = None
    for attempt in range(3):
        try:
            secs = sorted({int(x.split()[-1]) for x in lk.search_tesscut(f"{ra} {dec}").table["mission"]})
            break
        except Exception:  # noqa: BLE001
            time.sleep(10)
    if secs is None:
        return [], "tesscut search failed"
    added, failed = [], []
    for sec in secs:
        if min_sector <= sec <= 200 and sec not in have:
            try:
                own_photometry(tic, sec, ra, dec)
                added.append(sec)
            except Exception as e:  # noqa: BLE001
                failed.append(f"S{sec}: {str(e)[:80]}")
    return added, "; ".join(failed)


# ------------------------------------------------------------------ other planets on the star
_TOI = None


def toi_table():
    global _TOI
    if _TOI is None:
        t = pd.read_csv(ROOT / "data" / "raw" / "exofop_toi.csv", comment="#", low_memory=False)
        _TOI = t
    return _TOI


SINGLE_RE = re.compile(r"single|duo|mono|maximum period|max period|max\. period|period could be|actual period",
                       re.I)


def sibling_ephemerides(tic, toi):
    """Other TOIs on the star with a usable period (P > 0 and comment not single/duo/max-period)."""
    t = toi_table()
    s = t[(t["TIC ID"] == tic) & (t["TOI"].astype(str) != str(toi))]
    out = []
    for _, r in s.iterrows():
        P, ep, dur = (pd.to_numeric(r.get(c), errors="coerce") for c in
                      ("Period (days)", "Epoch (BJD)", "Duration (hours)"))
        com = str(r.get("Comments", ""))
        if not (np.isfinite(P) and P > 0 and np.isfinite(ep) and np.isfinite(dur)):
            continue
        if SINGLE_RE.search(com):
            continue
        out.append(dict(TOI=str(r["TOI"]), P=float(P), T0=float(ep) - 2457000, dur_h=float(dur),
                        disp=str(r.get("TFOPWG Disposition"))))
    return out


def mask_siblings(lc, sibs):
    if not sibs or len(lc) == 0:
        return lc, 0
    t = lc.time.to_numpy()
    bad = np.zeros(len(t), bool)
    for s in sibs:
        ph = (t - s["T0"] + 0.5 * s["P"]) % s["P"] - 0.5 * s["P"]
        bad |= np.abs(ph) < 0.75 * s["dur_h"] / 24 + 1 / 24
    return lc[~bad].reset_index(drop=True), int(bad.sum())


# ------------------------------------------------------------------ fine period constraints
def fine_period_constraints(g, epoch, D, pmin=5.0, pmax=1500.0):
    """Exclusion of trial periods by flat windows on a period grid fine enough that consecutive
    trial periods move the furthest predicted epoch by <= D/8."""
    Tg = g.T_pred.to_numpy()
    flat = (g.verdict == "excluded").to_numpy()
    span = max(abs(Tg.max() - epoch), abs(Tg.min() - epoch))
    Ps = [pmin]
    while Ps[-1] < pmax:
        P = Ps[-1]
        kmax = max(1.0, span / P)
        Ps.append(P + D / 8 / kmax)
    Ps = np.array(Ps)
    tlo, thi = Tg.min(), Tg.max()
    status = np.zeros(len(Ps), bool)
    for i, P in enumerate(Ps):
        ks = np.arange(np.ceil((tlo - epoch) / P), np.floor((thi - epoch) / P) + 1)
        ks = ks[ks != 0]
        if not len(ks):
            continue
        T = epoch + ks * P
        j = np.clip(np.searchsorted(Tg, T), 1, len(Tg) - 1)
        jn = np.where(np.abs(Tg[j - 1] - T) < np.abs(Tg[j] - T), j - 1, j)
        status[i] = np.any(flat[jn] & (np.abs(Tg[jn] - T) <= D / 8 + 1e-9))
    w = np.gradient(np.log(Ps))
    allowed = ~status
    ranges, start = [], None
    for i, ok in enumerate(allowed):
        if ok and start is None:
            start = Ps[i]
        if (not ok or i == len(Ps) - 1) and start is not None:
            ranges.append((float(start), float(Ps[i - 1] if not ok else Ps[i])))
            start = None
    return dict(frac_logP_excluded_fine=float(np.sum(w[status]) / np.sum(w)),
                min_allowed_P_fine=float(Ps[allowed].min()) if allowed.any() else None,
                n_allowed_ranges=len(ranges), allowed_ranges=ranges)


def summarize_ranges(ranges, maxn=6):
    """Compact text: the widest allowed ranges in log P."""
    if not ranges:
        return "none"
    r = sorted(ranges, key=lambda x: -np.log(x[1] / x[0]))[:maxn]
    r = sorted(r)
    return "; ".join(f"{a:.1f}-{b:.1f}" for a, b in r) + (f" (+{len(ranges) - maxn} more)" if len(ranges) > maxn
                                                          else "")


# ------------------------------------------------------------------ one TOI
def run_toi(row, out, do_diff=True, own_min_sector=1, fetch=True):
    """row: dict with TOI, TIC, ra, dec, Tmag, epoch_btjd, dur_h, depth_ppm. Returns a JSON-able dict."""
    tic = int(row["TIC"])
    D, delta, epoch = row["dur_h"] / 24, row["depth_ppm"] * 1e-6, float(row["epoch_btjd"])
    rep = dict(TOI=str(row["TOI"]), TIC=tic, Tmag=row.get("Tmag"), depth_ppm=row["depth_ppm"], dur_h=row["dur_h"],
               epoch_btjd=epoch)
    if fetch:
        rep["fetch"] = fetch_lightcurves(tic)
        rep["own_photometry_sectors"], rep["own_photometry_failed"] = add_own_photometry(
            tic, row["ra"], row["dec"], own_min_sector)
    lc = ap.load_lightcurves(tic, products=PRODUCTS)
    rep["sectors_loaded"] = " ".join(f"{s}:{p}" for s, p in lc.groupby("sector")["product"].first().items())
    rep["n_sectors"] = int(lc.sector.nunique())
    sibs = sibling_ephemerides(tic, row["TOI"])
    lc, nmask = mask_siblings(lc, sibs)
    rep["siblings_masked"] = [s["TOI"] for s in sibs]
    rep["n_masked_cadences"] = nmask
    if len(lc) == 0:
        rep["status"] = "no light curves"
        return rep
    t, f, sg = lc.time.to_numpy(), lc.flux.to_numpy(), lc.segment.to_numpy()
    ev = ap.test_window(t, f, epoch, pilot.EPOCH_ERR, D, delta, sg)
    rep["event_verdict"] = ev.get("verdict")
    rep["event_depth_ppm"] = round(ev["depth"] * 1e6, 1) if "depth" in ev else None
    rep["event_snr"] = round(ev["snr"], 1) if "snr" in ev else None
    ev_T = epoch + (ev.get("offset", 0.0) or 0.0)
    rep["event_T"] = float(ev_T)
    g = k2lib.grid(lc, D, delta, ev_T)
    reps = k2lib.find_repeats(g, delta, D)
    # Re-anchoring (added before the validation answers were frozen, prompted by our own event check, not by
    # any reference period): a TOI epoch can be a propagated reference epoch that is not itself a transit in
    # the data. If the catalogued event is not 'supported' but robust dips exist, the highest-SNR dip becomes
    # the anchor event and the repeat search is redone relative to it.
    rep["re_anchored"] = False
    if rep["event_verdict"] != "supported" and reps:
        a = max(reps, key=lambda x: x["snr"])
        rep.update(re_anchored=True, toi_epoch_offset_d=round(a["T"] - epoch, 4))
        ev_T = float(a["T"])
        g["near_event"] = (g.T_pred - ev_T).abs() < D
        reps = k2lib.find_repeats(g, delta, D)
        rep["event_T"] = ev_T
    rep["repeats"] = reps
    rep["n_repeats"] = len(reps)
    rep["periods"] = k2lib.robust_periods(g, ev_T, reps, D)
    rep["n_periods"] = len(rep["periods"])
    rep["max_support"] = max([p["n_support"] for p in rep["periods"]], default=0)
    rep["frac_logP_excluded"] = k2lib.frac_excluded(g, ev_T, D)
    rep.update(fine_period_constraints(g, ev_T, D))
    rep["n_windows"] = int(len(g))
    rep["n_excluded_windows"] = int((g.verdict == "excluded").sum())
    g.to_csv(out / f"{tic}_{rep['TOI']}_grid.csv.gz", index=False)
    if do_diff:
        prow = pd.Series(dict(RA=row["ra"], Dec=row["dec"]))
        rep["diff_images"] = pilot.compute_diffs(tic, prow, lc, ev_T, reps, D, out)
    pilot.plot_windows(tic, lc, ev_T, reps, D, delta, out)
    rep["status"] = "ok"
    return rep


def run_safe(row, out, **kw):
    f = out / f"{row['TOI']}.json"
    if f.exists():
        return json.loads(f.read_text())
    try:
        rep = run_toi(row, out, **kw)
    except Exception:  # noqa: BLE001
        rep = dict(TOI=str(row["TOI"]), TIC=int(row["TIC"]), status="error", error=traceback.format_exc()[-800:])
    f.write_text(json.dumps(rep, indent=1, default=float))
    return rep
