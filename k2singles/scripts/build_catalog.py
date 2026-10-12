"""
K2 single transits re-observed by TESS: build the candidate list and check feasibility.

Sources (arXiv LaTeX in data/raw/k2singles/):
  LaCourse & Jacobs 2018, RNAAS 2, 28 (arXiv:1802.05823): table of single transits and eclipses
    in K2 C0-C17. Columns: EPIC, campaign, Kp, Teff, R*, BJD0 (BJD-2454833), depth (ppm),
    duration (h), comments.
  (Osborn et al. 2016, arXiv:1512.03722, adds 7 C1-C3 candidates; not parsed here yet.)

Steps:
  1. parse the LaCourse & Jacobs table -> epoch in BTJD (= BJD - 2457000 = BJD0 - 2167)
  2. EPIC coordinates (VizieR IV/34) -> nearest TIC within 3" (MAST TIC v8.2)
  3. TESS coverage S1-S134 (tess-point); sectors already observed (end date before today)
  4. expected per-transit TESS SNR = depth / sigma_1h * sqrt(duration_h), with
     sigma_1h = 220 ppm * 10^(0.2 (Tmag - 9.1)). Calibrated on our own TIC 229055790 recovery
     (Tmag 9.1, 1100 ppm, 6.5 h, SNR ~12.6 in an independent reduction).
  5. already-known status: TOI / CTOI / confirmed planet on the TIC (ExoFOP/pscomppars
     tables in data/raw/), noted, not dropped

Output: k2singles/data/k2_singles.csv, printed summary.
Run:    python k2singles/scripts/build_catalog.py
"""
import re
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "k2singles" / "data"


def parse_lacourse():
    tex = (RAW / "k2singles" / "1802.05823" / "ms.tex").read_text(encoding="utf-8", errors="ignore")
    body = tex.split("\\startdata", 1)[1].split("\\enddata", 1)[0]
    rows = []
    for line in body.splitlines():
        line = line.strip()
        if not line or line.startswith("%") or "&" not in line:
            continue
        f = [x.strip() for x in line.rstrip("\\").split("&")]
        if len(f) < 8 or not re.match(r"^\d{9}$", f[0]):
            continue
        num = lambda x: float(x) if re.match(r"^-?[\d.]+$", x) else np.nan  # noqa: E731
        rows.append(dict(EPIC=int(f[0]), campaign=f[1], Kp=num(f[2]), Teff=num(f[3]), Rstar=num(f[4]),
                         bkjd=num(f[5]), depth_ppm=num(f[6]), dur_h=num(f[7]),
                         comment=f[8] if len(f) > 8 else ""))
    d = pd.DataFrame(rows)
    d["epoch_btjd"] = d.bkjd + 2454833 - 2457000
    d["source"] = "LaCourse & Jacobs 2018"
    return d


def epic_coords(epics):
    from astroquery.vizier import Vizier
    v = Vizier(columns=["ID", "RAJ2000", "DEJ2000", "Kpmag"], row_limit=-1)
    out = []
    for e in epics:  # VizieR's batched OR syntax returns nothing for this catalogue; query one by one
        for attempt in range(3):
            try:
                r = v.query_constraints(catalog="IV/34/epic", ID=str(e))
                break
            except Exception:  # noqa: BLE001
                time.sleep(3)
        else:
            continue
        if r:
            out.append(r[0].to_pandas())
    c = pd.concat(out).rename(columns={"ID": "EPIC", "RAJ2000": "ra", "DEJ2000": "dec"})
    c["EPIC"] = c.EPIC.astype(int)
    return c[["EPIC", "ra", "dec", "Kpmag"]].drop_duplicates("EPIC")


def tic_match(ra, dec):
    from astroquery.mast import Catalogs
    for attempt in range(3):
        try:
            t = Catalogs.query_region(f"{ra} {dec}", radius=3 / 3600, catalog="TIC").to_pandas()
            break
        except Exception:  # noqa: BLE001
            time.sleep(5)
    else:
        return None
    if not len(t):
        return None
    t = t.sort_values("dstArcSec")
    r = t.iloc[0]
    return dict(TIC=int(r["ID"]), Tmag=float(r["Tmag"]), sep_arcsec=float(r["dstArcSec"]),
                n_tic_3as=int(len(t)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = parse_lacourse()
    print(f"parsed {len(d)} rows from LaCourse & Jacobs 2018")
    c = epic_coords(d.EPIC.tolist())
    d = d.merge(c, on="EPIC", how="left")
    print(f"EPIC coordinates: {d.ra.notna().sum()}/{len(d)}")
    ms = []
    for r in d.itertuples():
        ms.append(tic_match(r.ra, r.dec) if np.isfinite(r.ra) else None)
    m = pd.DataFrame([x if x else {} for x in ms])
    d = pd.concat([d.reset_index(drop=True), m], axis=1)
    print(f"TIC matches within 3\": {d.TIC.notna().sum()}/{len(d)}")

    import tess_stars2px as tp
    ends = pd.read_csv(ROOT / "singles" / "data" / "tess_orbit_times_latest.csv", comment="#")
    ends["end"] = pd.to_datetime(ends["End of Orbit"], utc=True)
    sec_end = ends.groupby("Sector").end.max()
    today = pd.Timestamp(date.today(), tz="UTC")
    done = {int(s) for s, e in sec_end.items() if e < today}
    obs, done_secs = [], []
    for r in d.itertuples():
        if not np.isfinite(r.ra):
            obs.append(""), done_secs.append("")
            continue
        o = tp.tess_stars2px_function_entry(int(r.TIC) if np.isfinite(r.TIC) else 0, r.ra, r.dec)
        s = sorted({int(x) for x in o[3]})
        obs.append(" ".join(map(str, s)))
        done_secs.append(" ".join(str(x) for x in s if x in done))
    d["tess_sectors_all"] = obs
    d["tess_sectors_observed"] = done_secs
    d["n_tess_observed"] = [len(x.split()) if x else 0 for x in done_secs]
    sig1h = 220 * 10 ** (0.2 * (d.Tmag - 9.1))
    d["tess_snr_per_transit"] = d.depth_ppm / sig1h * np.sqrt(d.dur_h)

    toi = pd.read_csv(RAW / "exofop_toi.csv")
    ctoi = pd.read_csv(RAW / "exofop_ctoi.csv")
    ps = pd.read_csv(RAW / "nasa_pscomppars.csv")
    ps_tic = set(pd.to_numeric(ps.tic_id.str.replace("TIC ", ""), errors="coerce").dropna().astype(int))
    toi_map = toi.groupby("TIC ID")["TOI"].apply(lambda x: " ".join(str(v) for v in x)).to_dict()
    ctoi_col = "TIC ID" if "TIC ID" in ctoi.columns else ctoi.columns[0]
    ctoi_set = set(pd.to_numeric(ctoi[ctoi_col], errors="coerce").dropna().astype(int))
    d["TOI"] = [toi_map.get(int(t), "") if np.isfinite(t) else "" for t in d.TIC]
    d["is_CTOI"] = [bool(np.isfinite(t) and int(t) in ctoi_set) for t in d.TIC]
    d["confirmed_host"] = [bool(np.isfinite(t) and int(t) in ps_tic) for t in d.TIC]
    d.to_csv(OUT / "k2_singles.csv", index=False)

    feas = d[(d.n_tess_observed > 0) & (d.tess_snr_per_transit >= 10)]
    print(f"\nwith >=1 TESS sector observed: {(d.n_tess_observed > 0).sum()}")
    print(f"... and per-transit TESS SNR >= 10: {len(feas)}")
    print(f"... of which flagged 'Deep' / eclipse-like in the comments: "
          f"{feas.comment.str.contains('Deep|EB|eclips', case=False).sum()}")
    print(f"... TOI on star: {(feas.TOI != '').sum()}, CTOI: {feas.is_CTOI.sum()}, confirmed host: {feas.confirmed_host.sum()}")
    print("\nTESS sectors observed (feasible ones): median", feas.n_tess_observed.median())
    pd.set_option("display.width", 220)
    print(feas.sort_values("tess_snr_per_transit", ascending=False)[
        ["EPIC", "TIC", "campaign", "Tmag", "depth_ppm", "dur_h", "Rstar", "comment", "n_tess_observed",
         "tess_snr_per_transit", "TOI", "is_CTOI", "confirmed_host"]].head(60).round(1).to_string(index=False))


if __name__ == "__main__":
    main()
