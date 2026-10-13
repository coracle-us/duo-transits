"""
Step 3: new TESS data since each TOI was last updated (quantification only, no search).

For every target:
  - all TESS sectors on silicon (tess-point, S1-S111 pointing table in data/raw/tess_orbit_times.csv)
  - observed = sector ended before RUN_DATE (2026-10-12)
  - new_since_update = observed sectors whose START is after 'Date TOI Updated' (data the last
    update cannot have used)
  - beyond_listed = observed sectors greater than the largest sector in the TOI 'Sectors' column
  - comment_last_sector: largest sector number named in the comment ("single transit in s46",
    "still single as of Sector 72"); flag_stale_comment = observed sectors > that number exist
Output: toi_singles_k2/new_tess.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "toi_singles_k2"
RUN_DATE = pd.Timestamp("2026-10-12", tz="UTC")


def main():
    import tess_stars2px as tp
    d = pd.read_csv(OUT / "targets.csv")
    ot = pd.read_csv(ROOT / "data" / "raw" / "tess_orbit_times.csv", comment="#")
    ot["start"] = pd.to_datetime(ot["Start of Orbit"], utc=True)
    ot["end"] = pd.to_datetime(ot["End of Orbit"], utc=True)
    sec = ot.groupby("Sector").agg(start=("start", "min"), end=("end", "max"))
    rows = []
    for r in d.itertuples():
        o = tp.tess_stars2px_function_entry(int(r.TIC), r.ra_deg, r.dec_deg)
        allsec = sorted({int(x) for x in o[3]})
        obs = [s for s in allsec if s in sec.index and sec.loc[s, "end"] < RUN_DATE]
        upd = pd.Timestamp(r.updated, tz="UTC")
        new = [s for s in obs if sec.loc[s, "start"] > upd]
        listed = [int(x) for x in str(r.Sectors).split(",") if x.strip().isdigit()]
        beyond = [s for s in obs if listed and s > max(listed)]
        cl = r.last_sector_in_comment
        stale = [s for s in obs if np.isfinite(cl) and s > cl]
        upcoming = [s for s in allsec if s in sec.index and sec.loc[s, "end"] >= RUN_DATE]
        rows.append(dict(TOI=r.TOI, TIC=r.TIC, Tmag=r.Tmag, period_d=r.period_d, max_period_d=r.max_period_d,
                         updated=r.updated, listed_sectors=" ".join(map(str, listed)),
                         observed_sectors=" ".join(map(str, obs)), n_observed=len(obs),
                         new_since_update=" ".join(map(str, new)), n_new_since_update=len(new),
                         beyond_listed=" ".join(map(str, beyond)), n_beyond_listed=len(beyond),
                         comment_last_sector=cl, n_after_comment_sector=len(stale),
                         flag_stale_comment=bool(len(stale)), upcoming_to_S111=" ".join(map(str, upcoming)),
                         Comments=r.Comments))
    n = pd.DataFrame(rows)
    n.to_csv(OUT / "new_tess.csv", index=False)
    print(f"targets: {len(n)}")
    print(f"with >=1 observed sector after 'Date TOI Updated': {(n.n_new_since_update > 0).sum()} "
          f"(total {n.n_new_since_update.sum()} sector-observations; median among those "
          f"{n[n.n_new_since_update > 0].n_new_since_update.median()})")
    print(f"with >=1 observed sector beyond the TOI 'Sectors' list: {(n.n_beyond_listed > 0).sum()} "
          f"(total {n.n_beyond_listed.sum()})")
    print(f"comment names a sector: {n.comment_last_sector.notna().sum()}; of those, later sectors observed: "
          f"{n.flag_stale_comment.sum()}")
    print("distribution of n_new_since_update:")
    print(n.n_new_since_update.value_counts().sort_index().to_string())
    print(f"no new data at all: {(n.n_new_since_update == 0).sum()}")
    yrs = pd.to_datetime(n.updated).dt.year.value_counts().sort_index()
    print("Date TOI Updated by year:\n" + yrs.to_string())


if __name__ == "__main__":
    main()
