"""Stage 1: pilot sample (100 non-EB BeT stars) and TESS light-curve download.

Product rule (one per star): rank SPOC 120 s > TESS-SPOC > QLP > GSFC-ELEANOR-LITE.
The star's preferred product is the highest-ranked one that covers >= 2 qualifying sectors
(else the highest-ranked present). A sector without the preferred product uses the next
available product and is flagged mixed_pipeline. Qualifying = sector mid-time before the
last NEOWISE epoch (read from the China-VO WISE curves). Skips CDIPS/PATHOS/TASOC/TGLC.
"""
import glob, sys, time, traceback, warnings
import numpy as np, pandas as pd
from astropy.io import fits
import lightkurve as lk
from astroquery.mast import Catalogs
from common import CAT_DIR, PILOT, RESULTS, TESS_DIR, SEED, BTJD_TO_MJD

warnings.filterwarnings("ignore")
N_TARGET = 100
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


def main():
    cat = pd.read_parquet(RESULTS / "catalogue.parquet")
    pool = cat[~cat.eb_type.isin([1, 2, 3])].reset_index(drop=True)
    order = np.random.default_rng(SEED).permutation(len(pool))
    cutoff = last_neowise_mjd()
    print(f"pool {len(pool)} non-EB stars; last NEOWISE epoch MJD {cutoff:.2f}", flush=True)
    draws, sectors = [], []
    printed_qlp = False
    t0 = time.time()
    for k, i in enumerate(order):
        if sum(d["status"] == "pass" for d in draws) >= N_TARGET:
            break
        s = pool.iloc[i]
        rec = dict(draw=k, tic=int(s.tic), gaia_dr3_id=int(s.gaia_dr3_id), status="fail",
                   reason="", n_sectors_any=0, n_sectors_pre=0, preferred="", n_mixed=0)
        try:
            sr = lk.search_lightcurve(f"TIC {int(s.tic)}", mission="TESS")
            if len(sr) == 0:
                rec["reason"] = "no TESS light curves"; draws.append(rec); continue
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
                draws.append(rec); continue
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
                lc = sr[int(r.idx)].download(download_dir=str(TESS_DIR), quality_bitmask="none")
                path = lc.meta.get("FILENAME") or ""
                hdr1 = fits.getheader(path, 1)
                if use == "QLP" and not printed_qlp:
                    print("QLP columns:", fits.getdata(path, 1).columns.names, flush=True)
                    printed_qlp = True
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
            traceback.print_exc()
        draws.append(rec)
        npass = sum(d["status"] == "pass" for d in draws)
        print(f"[{time.time()-t0:6.0f}s] draw {k} TIC {int(s.tic)} {rec['status']} "
              f"{rec['preferred']} pre={rec['n_sectors_pre']} {rec['reason']} | pass {npass}", flush=True)
        pd.DataFrame(draws).to_csv(PILOT / "sample.csv", index=False)
        pd.DataFrame(sectors).to_csv(PILOT / "sectors.csv", index=False)
    pd.DataFrame(draws).to_csv(PILOT / "sample.csv", index=False)
    pd.DataFrame(sectors).to_csv(PILOT / "sectors.csv", index=False)
    print("done", time.time() - t0, "s", flush=True)


if __name__ == "__main__":
    main()
