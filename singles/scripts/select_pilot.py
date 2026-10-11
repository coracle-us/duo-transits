"""
Single-transit project, step 2: choose pilot targets.

For candidates in singles/data/ctoi_singles.csv (brightest first, Tmag >= TMAG_MIN to avoid
saturated FFIs):
  - sectors in which TESS observed the star (lightkurve.search_tesscut)
  - sector time ranges from the TESS orbit-times table (data/raw/tess_orbit_times.csv)
  - sectors observed AFTER the CTOI was last modified ("new" data the uploader never had)
  - ExoFOP per-target JSON (cached in data/raw/exofop_targets/): drop any star with a
    community planet-parameter upload carrying a period, or a period note
Stop when N_PILOT targets qualify. Output: singles/data/pilot_targets.csv (+ a log of skips).

Options:  --by-snr      order by expected event SNR, depth_ppm * 10^(-0.2 (Tmag - 10)), instead of
                        brightness (the brightest-first pilot picked mostly shallow events)
          --exclude F   skip TICs listed in CSV file F (e.g. an earlier batch)
          --out NAME    output name stem (default: pilot_targets)
          --depth-min X / --depth-max Y   keep depth_ppm within [X, Y] (default 200-50000)
          --gaia-clean  keep only stars with Gaia DR3 flags computed (data/gaia_binarity_flags.csv)
                        and none raised: no NSS orbit, non_single_star == 0, not RV-variable,
                        RUWE <= 1.4. Batch 2 showed deep events on such stars are mostly EBs.

Usage:  python singles/scripts/select_pilot.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import fetch_exofop as fx  # noqa: E402

DATA = ROOT / "singles" / "data"
TMAG_MIN = 7.0
N_PILOT = 20
MAX_SECTOR = 120


def sector_times():
    """Sector start/end in BTJD from the TESS orbit-times table (tess.mit.edu, saved in
    data/raw/tess_orbit_times.csv; UTC, so BTJD is approximate to ~10 min, fine here).
    (Querying MAST per sector took ~3 min/query, so it is not used.)"""
    o = pd.read_csv(ROOT / "data" / "raw" / "tess_orbit_times.csv", comment="#")
    o["t0"] = pd.to_datetime(o["Start of Orbit"]).map(pd.Timestamp.to_julian_date) - 2457000
    o["t1"] = pd.to_datetime(o["End of Orbit"]).map(pd.Timestamp.to_julian_date) - 2457000
    g = o.groupby("Sector").agg(t_min_btjd=("t0", "min"), t_max_btjd=("t1", "max"))
    g.index.name = "sector"
    g.to_csv(DATA / "sector_times.csv")
    return g


def community_period(tic):
    """True if ExoFOP has a non-TOI, non-archive planet-parameter row with a period, or the
    observing notes mention a period/alias solution."""
    d = json.loads((fx.OUT / f"{tic}.json").read_text(encoding="utf-8"))
    prov, hits = None, []
    for p in d.get("planet_parameters", []):
        if "prov" in p:
            prov = p.get("prov_title", "")
            continue
        per = pd.to_numeric(p.get("per"), errors="coerce")
        if prov and not prov.startswith("TOIs") and per and per > 0:
            hits.append(f"{per:.4f} d by {p.get('puser')} ({str(p.get('pdate'))[:10]}): {p.get('pnotes')}")
    notes = fx.obsnotes_text(tic).lower()
    if any(k in notes for k in ("period", "alias", "second transit", "additional transit")):
        hits.append("observing notes mention period/alias/additional transit")
    return hits


def main():
    import lightkurve as lk
    args = sys.argv[1:]
    by_snr = "--by-snr" in args
    exclude = set()
    if "--exclude" in args:
        exclude = set(pd.read_csv(args[args.index("--exclude") + 1]).TIC.astype(int))
    stem = args[args.index("--out") + 1] if "--out" in args else "pilot_targets"
    dmin = float(args[args.index("--depth-min") + 1]) if "--depth-min" in args else 200.0
    dmax = float(args[args.index("--depth-max") + 1]) if "--depth-max" in args else 50000.0
    st = sector_times()
    cand = pd.read_csv(DATA / "ctoi_singles.csv")
    cand = cand[(cand.Tmag >= TMAG_MIN) & ~cand.TIC.isin(exclude)].drop_duplicates("TIC")
    if by_snr:
        # depths > 5% are implausible for a planet (several CTOIs list 95-99%, a units/entry
        # error) and < 200 ppm are hopeless for difference imaging in FFIs
        bad = (cand.depth_ppm > dmax) | (cand.depth_ppm < dmin)
        cand[bad][["TIC", "CTOI", "User", "Tmag", "depth_ppm"]].to_csv(
            DATA / f"{stem}_excluded_depth.csv", index=False)
        print(f"excluded {bad.sum()} targets with depth outside {dmin:.0f}-{dmax:.0f} ppm")
        cand = cand[~bad].copy()
        if "--gaia-clean" in args:
            g = pd.read_csv(DATA / "gaia_binarity_flags.csv")
            g = g[g.source_id.notna()]
            nss = g["nss_nss_solution_type"].notna() if "nss_nss_solution_type" in g else False
            clean = g[~nss & (g.non_single_star.fillna(1) == 0) & ~g.rv_variable.fillna(False).astype(bool)
                      & (g.ruwe.fillna(99) <= 1.4)]
            n0 = len(cand)
            cand = cand[cand.TIC.isin(set(clean.tic.astype(int)))].copy()
            print(f"Gaia-clean filter: {len(cand)} of {n0} remain (flags computed for {len(g)} stars)")
        cand["snr_score"] = cand.depth_ppm * 10 ** (-0.2 * (cand.Tmag - 10))
        cand = cand.sort_values("snr_score", ascending=False)
    else:
        cand = cand.sort_values("Tmag")
    fx.OUT.mkdir(parents=True, exist_ok=True)
    picked, log = [], []
    for _, r in cand.iterrows():
        tic = int(r.TIC)
        try:
            sr = lk.search_tesscut(f"TIC {tic}")
            secs = sorted({int(m.split()[-1]) for m in sr.table["mission"]})
        except Exception as e:  # noqa: BLE001
            log.append(dict(TIC=tic, skip=f"tesscut search failed: {e}"))
            continue
        lastmod_btjd = (pd.Timestamp(r["CTOI lastmod"]).to_julian_date() - 2457000)
        ev_sec = [s for s in secs if s in st.index and st.t_min_btjd[s] - 1 <= r.epoch_btjd <= st.t_max_btjd[s] + 1]
        new = [s for s in secs if s in st.index and st.t_min_btjd[s] > lastmod_btjd]
        after_event = [s for s in secs if s in st.index and st.t_min_btjd[s] > r.epoch_btjd + 1]
        if not new:
            log.append(dict(TIC=tic, skip="no sectors after CTOI lastmod", sectors=" ".join(map(str, secs))))
            continue
        try:
            fx.fetch(tic)
            hits = community_period(tic)
        except Exception as e:  # noqa: BLE001
            log.append(dict(TIC=tic, skip=f"ExoFOP fetch failed: {e}"))
            continue
        if hits:
            log.append(dict(TIC=tic, skip="community period on ExoFOP: " + " | ".join(hits)[:300]))
            continue
        picked.append(dict(TIC=tic, CTOI=r.CTOI, User=r.User, Tmag=r.Tmag, RA=r.RA, Dec=r.Dec,
                           epoch_btjd=r.epoch_btjd, depth_ppm=r.depth_ppm, dur_h=r.dur_h,
                           event_sector=" ".join(map(str, ev_sec)), all_sectors=" ".join(map(str, secs)),
                           sectors_after_lastmod=" ".join(map(str, new)),
                           n_sectors_after_event=len(after_event), lastmod=r["CTOI lastmod"]))
        print(f"  pick {len(picked):2d}: TIC {tic} Tmag {r.Tmag:.2f} depth {r.depth_ppm:.0f} ppm "
              f"dur {r.dur_h:.1f} h event S{ev_sec} new {new}")
        if len(picked) >= N_PILOT:
            break
    pd.DataFrame(picked).to_csv(DATA / f"{stem}.csv", index=False)
    pd.DataFrame(log).to_csv(DATA / f"{stem}_selection_log.csv", index=False)
    print(f"picked {len(picked)}; skipped {len(log)}")
    if log:
        print(pd.DataFrame(log).skip.str.split(":").str[0].value_counts().to_string())


if __name__ == "__main__":
    main()
