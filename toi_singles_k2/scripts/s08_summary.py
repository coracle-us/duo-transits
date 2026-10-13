"""
Step 5: merge targets + new-TESS counts + ExoFOP uploads into one per-TOI table and print the
summary counts used in REPORT.txt.

An upload 'belongs' to the TOI when its toi field equals 'TOI <number>' (exact planet). A
'post-update community period' is such an upload with period > 0, dated after 'Date TOI Updated',
by a user other than the TESS project, and differing from the TOI-table period by > 0.1%.
Output: toi_singles_k2/targets_merged.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[2] / "toi_singles_k2"


def main():
    t = pd.read_csv(OUT / "targets.csv")
    n = pd.read_csv(OUT / "new_tess.csv")[["TOI", "observed_sectors", "n_observed", "new_since_update",
                                           "n_new_since_update", "n_beyond_listed", "flag_stale_comment"]]
    s = pd.read_csv(OUT / "targets_status.csv")[["TOI", "n_uploads", "in_pscomppars", "n_spectra", "n_imaging",
                                                 "n_time_series"]]
    u = pd.read_csv(OUT / "exofop_uploads.csv")
    u["date_ts"] = pd.to_datetime(u.date, errors="coerce", format="mixed")
    rows = []
    for r in t.itertuples():
        own = u[(u.TOI == r.TOI) & (u.toi_field.astype(str).str.strip() == f"TOI {r.TOI:.2f}")]
        upd = pd.Timestamp(r.updated)
        post = own[(own.period > 0) & (own.date_ts > upd) & (own.user.astype(str) != "TESS project")]
        if r.period_d > 0:
            post = post[(post.period / r.period_d - 1).abs() > 1e-3]
        rows.append(dict(TOI=r.TOI, n_own_uploads=len(own), n_post_update_periods=len(post),
                         post_update_periods=" ".join(f"{p:.4f}" for p in post.period),
                         post_update_users=" ".join(sorted(set(post.user.astype(str)))),
                         post_update_notes=" | ".join(str(x)[:120] for x in post.notes.fillna(""))))
    m = t.merge(n, on="TOI").merge(s, on="TOI").merge(pd.DataFrame(rows), on="TOI")
    m["group"] = np.select(
        [m.comment_says_resolved | (m.n_post_update_periods > 0),
         m.n_new_since_update > 0],
        ["community/comment period after last TOI update", "new TESS sectors since update, no newer period"],
        default="no new TESS data since update")
    m.to_csv(OUT / "targets_merged.csv", index=False)
    print(m.group.value_counts().to_string())
    g = m[m.group == "new TESS sectors since update, no newer period"]
    print(f"\n'new data, no newer period' group: {len(g)} TOIs; sectors-after-update total {g.n_new_since_update.sum()}; "
          f">=3 new sectors: {(g.n_new_since_update >= 3).sum()}")
    print(f"  of which period column = 0 (true singles): {(g.period_d == 0).sum()}; max-period/duo: "
          f"{(g.max_period_d.notna()).sum()}")
    print(f"  host in pscomppars: {g.in_pscomppars.sum()}")
    print(f"stale-comment flag (comment names sector N, later sectors observed): {m.flag_stale_comment.sum()}")
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 70)
    print(g.sort_values("n_new_since_update", ascending=False)[
        ["TOI", "TIC", "Tmag", "period_d", "max_period_d", "depth_ppm", "dur_h", "updated", "new_since_update",
         "Comments"]].head(25).to_string(index=False))
    print(m[m.group.str.startswith("community")][["TOI", "period_d", "post_update_periods", "post_update_users",
                                                 "post_update_notes"]].to_string(index=False))


if __name__ == "__main__":
    main()
