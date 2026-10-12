"""
K2 singles x TESS, full run: every LaCourse & Jacobs (2018) event with TESS coverage and expected
TESS SNR >= 7. That covers both the planet-like ones (re-run) and the eclipses. K2 light curves from
ALL campaigns of the star are added to the search.

Why K2 too: in the pilot, EPIC 211924561's K2 C16 transit was hidden from the standard products
(EVEREST flagged it as outliers, K2SFF had a gap). Re-observed fields (C5/C16/C18, C6/C17, ...)
give extra repeats and flat-window period exclusions.

K2 data: the K2 mission pipeline light curves (author 'K2', PDCSAP flux, lightkurve 'default'
quality bitmask, which removes thruster-firing cadences) for every campaign. They are cached as
k2singles/cache/<EPIC>/k2_c<NN>.npz and added to the TESS light curve with sector = 1000 + campaign,
product 'K2PDC' and one segment per gap > 0.5 d. pilot.run_target then searches K2 and TESS
together. Difference images are attempted only for TESS sectors; K2 ones fail gracefully.

Selection (fixed in advance): n_tess_observed >= 1, tess_snr_per_transit >= 7, and catalogued
depth < 1e6 ppm (some rows carry impossible depths > 100%).
K2-campaign period floor: as in run_pilot.py (only the campaign containing the epoch).
Output: k2singles/results/full/<TIC>.json, summary.csv. Resumable.
Usage:  python k2singles/scripts/run_full.py [n_workers]
"""
import json
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
sys.path.insert(0, str(ROOT / "scripts"))
OUT = ROOT / "k2singles" / "results" / "full"
CACHE = ROOT / "k2singles" / "cache"


def select():
    d = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv").drop_duplicates("EPIC")
    return d[(d.n_tess_observed > 0) & (d.tess_snr_per_transit >= 7) & (d.depth_ppm < 1e6)].copy()


K2_FOLDERS = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "102", "111", "112", "12", "13", "14", "15", "16",
              "17", "18", "19"]   # archive folder names (C9 absent; C10 only as c102; C11 split)


def _get(url, dest, tries=3):
    import time
    import urllib.error
    import urllib.request
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                dest.write_bytes(r.read())
            return True
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return False          # this EPIC was not observed in that campaign
            time.sleep(5 * (k + 1))
        except Exception:  # noqa: BLE001
            time.sleep(5 * (k + 1))
    raise RuntimeError(f"download failed after {tries} tries: {url}")


def k2_lightcurves(epic):
    """All K2 campaigns of this EPIC: K2 mission pipeline PDCSAP, default quality mask, cached.
    Files come straight from the MAST archive by path. The MAST search portal was timing out, and
    the first version wrongly cached a failed search as 'no K2 data' (bug fixed 2026-10-12).
    Network failures raise, so nothing is cached for that star."""
    from astropy.io import fits
    from lightkurve.utils import KeplerQualityFlags
    cdir = CACHE / str(epic)
    (cdir / "dl").mkdir(parents=True, exist_ok=True)
    (cdir / "none.txt").unlink(missing_ok=True)
    s = str(epic)
    if not (cdir / "checked.txt").exists():
        for c in K2_FOLDERS:
            fn = cdir / "dl" / f"ktwo{s}-c{int(c):02d}_llc.fits"
            if fn.exists():
                continue
            url = (f"https://archive.stsci.edu/missions/k2/lightcurves/c{c}/{s[:4]}00000/{s[4:6]}000/"
                   f"ktwo{s}-c{int(c):02d}_llc.fits")
            _get(url, fn)
        for fn in sorted((cdir / "dl").glob("ktwo*_llc.fits")):
            camp = int(fn.stem.split("-c")[1].split("_")[0])
            with fits.open(fn) as h:
                d = h[1].data
                t = np.asarray(d["TIME"], float) + 2454833 - 2457000
                f = np.asarray(d["PDCSAP_FLUX"], float)
                q = np.asarray(d["SAP_QUALITY"], int)
            ok = np.isfinite(t) & np.isfinite(f) & ((q & KeplerQualityFlags.DEFAULT_BITMASK) == 0)
            if ok.sum() > 50:
                np.savez(cdir / f"k2_c{camp:03d}.npz", t=t[ok], f=f[ok] / np.nanmedian(f[ok]))
        (cdir / "checked.txt").write_text("all K2 campaign folders checked")
    done = sorted(cdir.glob("k2_c*.npz"))
    parts = []
    for fn in done:
        camp = int(fn.stem.split("_c")[1])
        z = np.load(fn)
        t, f = z["t"], z["f"]
        o = np.argsort(t)
        t, f = t[o], f[o]
        sec = 1000 + camp
        seg = sec * 1000 + np.concatenate([[0], np.cumsum(np.diff(t) > 0.5)])
        parts.append(pd.DataFrame({"time": t, "flux": f, "err": np.nan, "sector": sec, "product": "K2PDC",
                                   "segment": seg}))
    return pd.concat(parts, ignore_index=True) if parts else None


def work(row_dict):
    import alias_pipeline as ap
    import pilot
    import run_pilot as rp
    pilot.OWN_MIN_SECTOR = 80
    row = pd.Series(row_dict)
    tic, epic = int(row.TIC), int(row.EPIC)
    f = OUT / f"{tic}.json"
    if f.exists():
        return tic, "exists"
    try:
        k2 = k2_lightcurves(epic)
        orig = ap.load_lightcurves

        def load_with_k2(t, max_sector=None, products=ap.PRODUCT_RANK):
            lc = orig(t, max_sector=max_sector, products=products)
            if k2 is None:
                return lc
            return pd.concat([lc, k2], ignore_index=True).sort_values("time").reset_index(drop=True)

        ap.load_lightcurves = load_with_k2
        pilot.ap.load_lightcurves = load_with_k2
        rep = pilot.run_target(row, OUT, True)
        # K2 coverage from the downloaded files themselves (no MAST portal call)
        iv = [] if k2 is None else [(float(g.time.min()), float(g.time.max())) for _, g in k2.groupby("sector")]
        cont = [(a, b) for a, b in iv if a - 0.5 <= float(row.epoch_btjd) <= b + 0.5]
        tmin, tmax = (cont[0] if cont else (None, None))
        rep.update(EPIC=epic, k2_campaign=str(row.campaign), k2_tmin_btjd=tmin, k2_tmax_btjd=tmax,
                   k2_intervals_btjd=iv, k2_campaigns_loaded=[] if k2 is None else sorted({int(s) - 1000 for s in k2.sector}),
                   catalog_comment=str(row.comment))
        if tmin is not None:
            pmin = max(float(row.epoch_btjd) - tmin, tmax - float(row.epoch_btjd))
            rep["P_min_from_K2"] = pmin
            cp = rep.get("candidate_periods") or []
            rep["candidate_periods_all"] = cp
            rep["candidate_periods"] = [x for x in cp if x[1] >= pmin]
            rep["n_candidate_periods"] = len(rep["candidate_periods"])
    except Exception:  # noqa: BLE001
        rep = dict(TIC=tic, EPIC=epic, status="error", error=traceback.format_exc()[-600:])
    f.write_text(json.dumps(rep, indent=1, default=float))
    reps = rep.get("repeats") or []
    return tic, (f"repeats {len(reps)} (K2 {sum(1 for x in reps if x['sector'] >= 1000)}) periods "
                 f"{[round(p, 2) for _, p in (rep.get('candidate_periods') or [])][:4]}")


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    OUT.mkdir(parents=True, exist_ok=True)
    t = select()
    rows = [dict(TIC=int(r.TIC), EPIC=int(r.EPIC), campaign=r.campaign, RA=float(r.ra), Dec=float(r.dec),
                 Tmag=float(r.Tmag), depth_ppm=float(r.depth_ppm), dur_h=float(r.dur_h),
                 epoch_btjd=float(r.epoch_btjd), CTOI=None, comment=str(r.comment)) for r in t.itertuples()]
    todo = [r for r in rows if not (OUT / f"{r['TIC']}.json").exists()]
    print(f"targets: {len(rows)}; to do: {len(todo)}", flush=True)
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = [ex.submit(work, r) for r in todo]
        for i, fu in enumerate(as_completed(futs), 1):
            tic, msg = fu.result()
            print(f"[{i}/{len(todo)}] TIC {tic}: {msg}", flush=True)
    summ = []
    for r in rows:
        fn = OUT / f"{r['TIC']}.json"
        if fn.exists():
            j = json.loads(fn.read_text())
            reps = j.get("repeats") or []
            summ.append(dict(TIC=r["TIC"], EPIC=r["EPIC"], Tmag=r["Tmag"], depth_ppm=r["depth_ppm"], dur_h=r["dur_h"],
                             comment=r["comment"], status=j.get("status", "ok"),
                             k2_campaigns=" ".join(map(str, j.get("k2_campaigns_loaded") or [])),
                             n_repeats=len(reps), n_repeats_k2=sum(1 for x in reps if x["sector"] >= 1000),
                             n_periods=j.get("n_candidate_periods"), P_min_K2=j.get("P_min_from_K2"),
                             frac_logP_excluded=j.get("frac_logP_excluded"),
                             periods=" ".join(f"{p:.3f}" for _, p in (j.get("candidate_periods") or [])[:6])))
    pd.DataFrame(summ).to_csv(OUT / "summary.csv", index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
