"""
Step 2a: which of the target stars were observed by K2?

For every target with |ecliptic latitude| < 25 deg (all K2 fields lie within ~10 deg of the ecliptic;
the margin is generous), two independent checks:
  1. VizieR IV/34 (EPIC) cone search, 3 arcsec -> EPIC ID(s), Kp
  2. MAST Observations.query_region (3 arcsec), obs_collection K2 / HLSP -> K2 light-curve products
     actually archived (EPIC catalogue membership alone does not mean the star was observed)
Then for each EPIC with K2 light curves, the K2 PDCSAP files are downloaded (archive URLs, as in
k2singles/scripts/k2lib.py, cached under toi_singles_k2/data/k2cache) and the campaign time ranges listed.

Alias / window test (per TOI on a K2 star):
  - with a max period Pmax (duotransit-like): aliases P_n = Pmax/n for n=1..n_max (P >= 5 d);
    for each alias the predicted transit times T = epoch + k*P_n are checked against K2 cadences
    (in-window if any cadence lies within +-T14/2 of T). Reports how many aliases put >=1 transit
    in K2 data.
  - single (no period): any period shorter than the epoch - K2 baseline gap could put a transit in
    K2; the K2 search is then blind (whole light curve). Reported as 'blind'.

Output: toi_singles_k2/k2_overlap.csv, toi_singles_k2/derived/epic_matches.csv
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "toi_singles_k2"
DATA = OUT / "data"
sys.path.insert(0, str(ROOT / "k2singles" / "scripts"))


def epic_cone(ra, dec, r_as=3.0):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier
    v = Vizier(columns=["ID", "RAJ2000", "DEJ2000", "Kpmag", "_r"], row_limit=-1)
    for k in range(3):
        try:
            r = v.query_region(SkyCoord(ra, dec, unit="deg"), radius=r_as * u.arcsec, catalog="IV/34/epic")
            return r[0].to_pandas() if len(r) else pd.DataFrame()
        except Exception:  # noqa: BLE001
            time.sleep(3 * (k + 1))
    return None


def mast_k2(ra, dec, r_as=3.0):
    from astroquery.mast import Observations
    for k in range(3):
        try:
            o = Observations.query_region(f"{ra} {dec}", radius=f"{r_as} arcsec").to_pandas()
            break
        except Exception:  # noqa: BLE001
            time.sleep(5 * (k + 1))
    else:
        return None
    if not len(o):
        return o
    keep = (o.obs_collection == "K2") | ((o.obs_collection == "HLSP") & o.target_name.astype(str).str.contains(
        "EPIC|ktwo|k2", case=False)) | o.provenance_name.astype(str).str.contains("K2|EVEREST|K2SFF|K2SC|POLAR|KEPLER",
                                                                                  case=False)
    keep &= o.dataproduct_type.isin(["timeseries", "image"]) | (o.obs_collection == "HLSP")
    return o[keep]


def main():
    import k2lib
    k2lib.CACHE = DATA / "k2cache"
    d = pd.read_csv(OUT / "targets.csv")
    near = d[d.ecl_lat.abs() < 25]
    print(f"{len(near)} TOIs within 25 deg of the ecliptic (of {len(d)})")
    rows = []
    for r in near.itertuples():
        e = epic_cone(r.ra_deg, r.dec_deg)
        m = mast_k2(r.ra_deg, r.dec_deg)
        epics = [] if e is None or not len(e) else e.sort_values("_r").ID.astype(int).tolist()
        mk = m if m is not None else pd.DataFrame()
        k2obs = (mk[(mk.obs_collection == "K2") & (mk.dataproduct_type == "timeseries")
                    & mk.target_name.astype(str).str.startswith("ktwo")] if len(mk) else mk)
        ffi = (sorted(set(mk[(mk.obs_collection == "K2") & (mk.dataproduct_type == "image")
                             ].sequence_number.dropna().astype(int))) if len(mk) else [])
        hl = sorted(set(mk[mk.obs_collection == "HLSP"].provenance_name.astype(str))) if len(mk) else []
        camps = sorted(set(k2obs.get("sequence_number", pd.Series(dtype=float)).dropna().astype(int))) if len(k2obs) else []
        names = sorted(set(k2obs.target_name.astype(str))) if len(k2obs) else []
        rows.append(dict(TOI=r.TOI, TIC=r.TIC, ra=r.ra_deg, dec=r.dec_deg, ecl_lat=r.ecl_lat,
                         epic_ids=" ".join(map(str, epics)),
                         epic_sep_as=(round(float(e.sort_values("_r")._r.iloc[0]), 2) if epics else np.nan),
                         Kp=(float(e.sort_values("_r").Kpmag.iloc[0]) if epics else np.nan),
                         mast_k2_targets=" ".join(names), mast_k2_campaigns=" ".join(map(str, camps)),
                         mast_hlsp=" ".join(hl), k2_ffi_campaigns=" ".join(map(str, ffi)), mast_query_ok=m is not None, vizier_ok=e is not None))
        print(r.TOI, rows[-1]["epic_ids"], rows[-1]["mast_k2_targets"], rows[-1]["mast_k2_campaigns"], hl,
              "FFI:", ffi, flush=True)
    em = pd.DataFrame(rows)
    em.to_csv(OUT / "derived" / "epic_matches.csv", index=False)

    # K2 light curves and campaign ranges, alias windows
    ov = []
    for r in em.itertuples():
        epics = [int(x.replace("ktwo", "").replace("EPIC", "").strip()) for x in r.mast_k2_targets.split()
                 if x.replace("ktwo", "").replace("EPIC", "").strip().isdigit()]
        epics = sorted(set(epics) | ({int(x) for x in r.epic_ids.split()} if r.mast_k2_targets else set()))
        if not epics:
            continue
        t = d[d.TOI == r.TOI].iloc[0]
        for epic in epics:
            lc = k2lib.k2_lightcurves(epic)
            if lc is None:
                ov.append(dict(TOI=r.TOI, TIC=r.TIC, EPIC=epic, campaigns="", note="EPIC in MAST/VizieR but no K2 llc"))
                continue
            camps, rng = [], []
            for sec, g in lc.groupby("sector"):
                camps.append(str(sec - 1000))
                rng.append(f"C{sec - 1000}:{g.time.min():.2f}-{g.time.max():.2f}")
            tt = np.sort(lc.time.to_numpy())
            D = t.dur_h / 24
            Pmax = t.max_period_d
            n_alias = n_hit = 0
            hit_list = []
            if np.isfinite(Pmax) and Pmax > 0:
                nmax = int(Pmax // 5)
                for n in range(1, nmax + 1):
                    P = Pmax / n
                    k = np.arange(np.floor((tt.min() - t.epoch_btjd) / P) - 1, np.ceil((tt.max() - t.epoch_btjd) / P) + 2)
                    T = t.epoch_btjd + k * P
                    T = T[(T > tt.min() - D) & (T < tt.max() + D)]
                    j = np.searchsorted(tt, T)
                    inw = [(abs(tt[min(i, len(tt) - 1)] - x) < D / 2) or (i > 0 and abs(tt[i - 1] - x) < D / 2)
                           for x, i in zip(T, j)]
                    n_alias += 1
                    if any(inw):
                        n_hit += 1
                        if len(hit_list) < 12:
                            hit_list.append(f"n={n}:P={P:.3f}")
                mode = "aliases"
            else:
                mode = "blind (no period)"
            ov.append(dict(TOI=r.TOI, TIC=r.TIC, EPIC=epic, campaigns=" ".join(camps), k2_ranges_btjd=" ".join(rng),
                           n_k2_cadences=len(lc), k2_span_d=round(float(sum(g.time.max() - g.time.min()
                                                                             for _, g in lc.groupby("sector"))), 1),
                           epoch_btjd=t.epoch_btjd, period_d=t.period_d, max_period_d=Pmax, dur_h=t.dur_h,
                           depth_ppm=t.depth_ppm, mode=mode, n_aliases=n_alias, n_aliases_with_k2_transit=n_hit,
                           example_aliases_in_k2=" ".join(hit_list), mast_hlsp=r.mast_hlsp, note=""))
    pd.DataFrame(ov).to_csv(OUT / "k2_overlap.csv", index=False)
    print(pd.DataFrame(ov).to_string())


if __name__ == "__main__":
    main()
