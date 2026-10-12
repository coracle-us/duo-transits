"""
Pick validation asteroids: LCDB U = 3 (secure) periods, crossing a 200-s-cadence sector (S56+).

Rule (fixed before any photometry was run):
  - LCDB 2023-Oct summary, U == '3', numbered, not flagged binary ('BIN' empty), AMAX >= 0.15 mag.
  - Period bins (h): [2,4), [4,8), [8,16), [16,30), [30,100), [100,400); 6 random objects per bin
    (numpy seed 1) among numbers < 20000 (bright enough), plus 3 Jupiter Trojans and 3 Hildas with
    U = 3 (from cache/population_lcdb_merge.csv.gz), also seed 1.
  - For each, TESScut moving-target sector list; keep sectors >= 56. Then JPL Horizons V at the
    sector mid-time from TESS; keep V <= 17.
  - Validation set = per period bin the first 2 objects (in random order) passing the cuts, plus
    the first 1 Trojan and 1 Hilda passing; the earliest qualifying sector is used.
Output: data/validation_candidates.csv (all queried, with sectors and V) and the chosen set
marked chosen=True.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lcdb_coverage as lc  # noqa: E402
import select_targets as st  # noqa: E402

BINS = [(2, 4), (4, 8), (8, 16), (16, 30), (30, 100), (100, 400)]


def sectors(num):
    from astroquery.mast import Tesscut
    try:
        s = Tesscut.get_sectors(objectname=f"{num};", moving_target=True)
        return sorted(set(int(x) for x in s["sector"]))
    except Exception as e:  # noqa: BLE001
        return f"ERR {e}"[:80]


def v_mid(num, sector):
    from astroquery.jplhorizons import Horizons
    t0, t1 = st.sector_span(sector)
    jd = 0.5 * (t0 + t1)
    try:
        e = Horizons(id=str(num), id_type="smallbody", location="500@-95", epochs=[jd]).ephemerides()
        return float(e["V"][0]), float(e["alpha"][0])
    except Exception:  # noqa: BLE001
        return np.nan, np.nan


def main():
    d = lc.read_lcdb()
    d = d[(d.U == "3") & d.number.notna() & (d.BIN == "") & (d.amax >= 0.15) & (d.number < 20000)]
    rng = np.random.default_rng(1)
    picks = []
    for lo, hi in BINS:
        g = d[(d.period_h >= lo) & (d.period_h < hi)]
        idx = rng.permutation(len(g))[:6]
        for j in idx:
            picks.append(dict(number=int(g.number.iloc[j]), name=g.NAME.iloc[j], P_lcdb_h=g.period_h.iloc[j],
                              amax=g.amax.iloc[j], group=f"P{lo}-{hi}h"))
    m = pd.read_csv(HERE / "cache" / "population_lcdb_merge.csv.gz")
    for pop in ["Jupiter Trojan", "Hilda"]:
        g = m[(m.population == pop) & (m.U == "3") & m.period_h.notna()]
        idx = rng.permutation(len(g))[:3]
        for j in idx:
            picks.append(dict(number=int(g.number.iloc[j]), name=g.name.iloc[j], P_lcdb_h=g.period_h.iloc[j],
                              amax=g.amax.iloc[j], group=pop))
    p = pd.DataFrame(picks)
    with ThreadPoolExecutor(6) as ex:
        p["sectors"] = list(ex.map(sectors, p.number))
    p["sector"] = [next((s for s in ss if s >= 56), np.nan) if isinstance(ss, list) else np.nan
                   for ss in p.sectors]
    vs = []
    for n, s in zip(p.number, p.sector):
        vs.append(v_mid(n, int(s)) if np.isfinite(s) else (np.nan, np.nan))
    p["V_mid"] = [v[0] for v in vs]
    p["alpha_mid"] = [v[1] for v in vs]
    ok = p.sector.notna() & (p.V_mid <= 17)
    p["chosen"] = False
    for grp, g in p[ok].groupby("group", sort=False):
        k = 1 if grp in ("Jupiter Trojan", "Hilda") else 2
        p.loc[g.index[:k], "chosen"] = True
    p.to_csv(HERE / "data" / "validation_candidates.csv", index=False)
    print(p.to_string(index=False))


if __name__ == "__main__":
    main()
