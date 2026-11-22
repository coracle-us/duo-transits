"""Stage 1: sample of non-EB BeT stars and TESS light-curve download.

Pilot (BE_RUN=pilot): draw in seed order until 100 stars pass (serial, as first run).
Full (BE_RUN=full): every non-EB star, in the same seed order, with a thread pool; resumable
(re-running skips stars already in sample.csv unless their draw ended in an error).

Product rule (one per star): rank SPOC 120 s > TESS-SPOC > QLP > GSFC-ELEANOR-LITE.
The star's preferred product is the highest-ranked one that covers >= 2 qualifying sectors
(else the highest-ranked present). A sector without the preferred product uses the next
available product and is flagged mixed_pipeline. Qualifying = sector mid-time before the
last NEOWISE epoch (read from the China-VO WISE curves). Skips CDIPS/PATHOS/TASOC/TGLC.
"""
import glob, time, warnings, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np, pandas as pd
from astropy.io import fits
import lightkurve as lk
from astroquery.mast import Catalogs
from common import CAT_DIR, PILOT, RESULTS, TESS_DIR, SEED, RUN

warnings.filterwarnings("ignore")
N_TARGET = 100 if RUN == "pilot" else None
N_WORKERS = 1 if RUN == "pilot" else 6
RANK = ["SPOC120", "TESS-SPOC", "QLP", "GSFC-ELEANOR-LITE"]
FLUX_COL = {"SPOC120": "PDCSAP_FLUX", "TESS-SPOC": "PDCSAP_FLUX", "QLP": "SAP_FLUX",
            "GSFC-ELEANOR-LITE": "CORR_FLUX"}


def last_neowise_mjd():
    return max(pd.read_csv(f, usecols=["mjd"]).mjd.max()
               for f in glob.glob(str(CAT_DIR / "wise_lc" / "wise_lc" / "*.csv")))


def product_of(row):
    a, e = row["author"], float(row["exptime"])
    if a == "SPOC":
        return "SPOC120" if abs(e - 120) < 1 else None
    return a if a in RANK else None


def contratio(tic):
    try:
        t = Catalogs.query_criteria(catalog="Tic", ID=int(tic))
        v = t["contratio"][0]
        return float(v) if np.isfinite(v) else np.nan
    except Exception:
        return np.nan


def retry(fn, n=3):
    for attempt in range(n):
        try:
            return fn()
        except Exception:
            if attempt == n - 1:
                raise
            time.sleep(5 * (attempt + 1))


def process(k, s, cutoff):
    rec = dict(draw=k, tic=int(s.tic), gaia_dr3_id=int(s.gaia_dr3_id), status="fail",
               reason="", n_sectors_any=0, n_sectors_pre=0, preferred="", n_mixed=0)
    sectors = []
    try:
        sr = retry(lambda: lk.search_lightcurve(f"TIC {int(s.tic)}", mission="TESS"))
        if len(sr) == 0:
            rec["reason"] = "no TESS light curves"
            return rec, sectors
        tb = sr.table.to_pandas()
        tb["idx"] = np.arange(len(tb))
        tb["product"] = [product_of(r) for _, r in tb.iterrows()]
        tb = tb[tb["product"].notna()].copy()
        # only the target itself (lightkurve can return neighbours within its radius)
        tb = tb[tb["target_name"].astype(str).str.lstrip("0") == str(int(s.tic))]
        tb["sector"] = tb["sequence_number"].astype(int)
        tb["mid_mjd"] = 0.5 * (tb.t_min + tb.t_max)
        tb = tb.drop_duplicates(["product", "sector"])
        rec["n_sectors_any"] = tb.sector.nunique()
        pre = tb[tb.mid_mjd < cutoff]
        rec["n_sectors_pre"] = pre.sector.nunique()
        if pre.sector.nunique() < 2:
            rec["reason"] = (f"<2 sectors before NEOWISE end in the 4 products "
                             f"({pre.sector.nunique()} pre, {tb.sector.nunique()} total)")
            return rec, sectors
        cover = {p: pre[pre["product"] == p].sector.nunique() for p in RANK}
        pref = next((p for p in RANK if cover[p] >= 2), next(p for p in RANK if cover[p] > 0))
        rec["preferred"] = pref
        rows = []
        for sec, g in pre.groupby("sector"):
            avail = set(g["product"])
            use = pref if pref in avail else next(p for p in RANK if p in avail)
            rows.append((sec, use, use != pref, g[g["product"] == use].iloc[0]))
        ctr = None
        for sec, use, mixed, r in rows:
            lc = retry(lambda: sr[int(r.idx)].download(download_dir=str(TESS_DIR),
                                                       quality_bitmask="none"))
            path = lc.meta.get("FILENAME") or ""
            hdr1 = fits.getheader(path, 1)
            crowd = hdr1.get("CROWDSAP", np.nan) if use in ("SPOC120", "TESS-SPOC") else np.nan
            if use in ("QLP", "GSFC-ELEANOR-LITE"):
                if ctr is None:
                    ctr = contratio(s.tic)
                crowd = 1.0 / (1.0 + ctr) if np.isfinite(ctr) else np.nan
            sectors.append(dict(tic=int(s.tic), sector=int(sec), product=use,
                                flux_col=FLUX_COL[use], cadence_s=float(r.exptime),
                                mixed_pipeline=bool(mixed), crowdsap=hdr1.get("CROWDSAP", np.nan),
                                contratio=ctr if ctr is not None else np.nan,
                                crowd=float(crowd) if crowd is not None else np.nan,
                                mid_mjd_mast=float(r.mid_mjd),
                                file=str(path).replace(str(TESS_DIR) + "/", "")))
        rec["n_mixed"] = sum(m for _, _, m, _ in rows)
        rec["status"] = "pass"
    except Exception as e:
        rec["reason"] = f"error: {type(e).__name__}: {e}"[:200]
        sectors = []
    return rec, sectors


def main():
    cat = pd.read_parquet(RESULTS / "catalogue.parquet")
    pool = cat[~cat.eb_type.isin([1, 2, 3])].reset_index(drop=True)
    order = np.random.default_rng(SEED).permutation(len(pool))
    cutoff = last_neowise_mjd()
    print(f"run {RUN}: pool {len(pool)} non-EB stars; last NEOWISE epoch MJD {cutoff:.2f}", flush=True)
    draws, sectors = [], []
    if RUN != "pilot" and (PILOT / "sample.csv").exists():
        old = pd.read_csv(PILOT / "sample.csv")
        old = old[~old.reason.fillna("").str.startswith("error")]
        draws = old.to_dict("records")
        if (PILOT / "sectors.csv").exists() and len(old):
            so = pd.read_csv(PILOT / "sectors.csv")
            sectors = so[so.tic.isin(old.tic)].to_dict("records")
    done = {d["tic"] for d in draws}
    lock = threading.Lock()
    t0 = time.time()

    def save():
        pd.DataFrame(draws).sort_values("draw").to_csv(PILOT / "sample.csv", index=False)
        pd.DataFrame(sectors).to_csv(PILOT / "sectors.csv", index=False)

    todo = [(k, pool.iloc[i]) for k, i in enumerate(order) if int(pool.iloc[i].tic) not in done]
    if N_TARGET is not None:      # pilot: serial, stop at N_TARGET passes
        for k, s in todo:
            if sum(d["status"] == "pass" for d in draws) >= N_TARGET:
                break
            rec, sec = process(k, s, cutoff)
            draws.append(rec); sectors += sec
            save()
    else:
        with ThreadPoolExecutor(N_WORKERS) as ex:
            futs = [ex.submit(process, k, s, cutoff) for k, s in todo]
            for n, f in enumerate(as_completed(futs)):
                rec, sec = f.result()
                with lock:
                    draws.append(rec); sectors += sec
                    if n % 25 == 0:
                        save()
                        npass = sum(d["status"] == "pass" for d in draws)
                        nerr = sum(str(d["reason"]).startswith("error") for d in draws)
                        print(f"[{time.time()-t0:6.0f}s] {len(draws)}/{len(pool)} drawn, "
                              f"{npass} pass, {nerr} errors", flush=True)
    save()
    print("done", time.time() - t0, "s", flush=True)


if __name__ == "__main__":
    main()
