"""
Frozen target list for the southern run (pre-registration prereg/southern_run_v1.md).

From data/sample_ranked_precision.csv (the 1,191 EBs passing the sample cuts, ranked by expected
LTTE amplitude error) keep:
  - Tmag <= 13.5: measured MAST light-curve availability (data/lc_availability_bins.csv) is
    ~100% of on-silicon sectors up to Tmag 13.5 and drops to 64-79% beyond 14;
  - >= 3 observing seasons (groups of on-silicon sectors separated by > 100 d, tess-point
    coverage up to S107), the minimum for the season-jackknife rule.
Columns added: n_seasons, est_files, in_pilot, in_validation.
Cost model: est_files = on-silicon sectors (<= S107); 0.31 MB per sector file (mean over the 4,636
files cached so far); timing ~8.3 s CPU per sector file + noise model ~45 s per star (measured).
Output: data/southern_run_targets.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
D = HERE / "data"
TMAG_MAX = 13.5
MIN_SEASONS = 3
MB_PER_FILE = 0.31
CPU_S_PER_FILE = 8.3   # run_etv, measured: 27 pilot/validation systems, 1,072 files, ~150 process-min
CPU_S_PER_STAR = 45.0  # noise_fit (bootstrap), measured median on the NCVZ sets


def main():
    r = pd.read_csv(D / "sample_ranked_precision.csv", comment="#")
    cov = pd.read_csv(D / "coverage_all.csv").fillna("").set_index("TIC").sectors
    st = pd.read_csv(D / "sector_table.csv").set_index("sector")

    def seasons(s):
        ss = sorted(int(x) for x in s.split(";") if x and int(x) in st.index)
        t = st.loc[ss, "t_start"].values
        return (1 + int(np.sum(np.diff(t) > 100)) if len(t) else 0), len(ss)

    sv = r.TIC.map(cov).apply(seasons)
    r["n_seasons"] = [x[0] for x in sv]
    r["est_files"] = [x[1] for x in sv]
    pilot = set(pd.read_csv(D / "pilot_list.csv", comment="#").TIC)
    val = set(pd.read_csv(D / "validation_list.csv").TIC)
    r["in_pilot"] = r.TIC.isin(pilot)
    r["in_validation"] = r.TIC.isin(val)
    t = r[(r.Tmag <= TMAG_MAX) & (r.n_seasons >= MIN_SEASONS)].reset_index(drop=True)
    with open(D / "southern_run_targets.csv", "w") as fh:
        fh.write(f"# southern run v1 targets: Tmag <= {TMAG_MAX}, >= {MIN_SEASONS} seasons; {len(t)} of {len(r)}; "
                 f"est {t.est_files.sum()} sector files, {t.est_files.sum() * MB_PER_FILE / 1024:.1f} GB, "
                 f"{(t.est_files.sum() * CPU_S_PER_FILE + len(t) * CPU_S_PER_STAR) / 3600:.0f} CPU-h\n")
        t.to_csv(fh, index=False, float_format="%.6g")
    print(f"{len(t)} targets of {len(r)}; files {t.est_files.sum()}; "
          f"{t.est_files.sum() * MB_PER_FILE / 1024:.1f} GB; {(t.est_files.sum() * CPU_S_PER_FILE + len(t) * CPU_S_PER_STAR) / 3600:.0f} CPU-h; "
          f"SCVZ {int((t.region == 'SCVZ').sum())}; pilot {int(t.in_pilot.sum())}; known multiples {int(t.known_multiple.notna().sum())}")


if __name__ == "__main__":
    main()
