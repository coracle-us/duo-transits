"""
Build the frozen candidate and host lists for the Gaia cross-match from tables already in this repo.

One row per candidate (or per period we derived ourselves):
    tic, catalogue, cand_id, t1_btjd, t2_btjd, period_d
  t1/t2: the catalogued transit times (t2 only for duotransits)
  period_d: a period the catalogue or this project states (blank if none)

Catalogues:
  hawthorn     data/catalog/candidates.csv                 Hawthorn et al. 2024 duotransits
  ctoi         singles/data/ctoi_singles.csv               stale single/duo CTOIs (ExoFOP 2026-10-11)
  salinas      salinas-second-transits/candidates.csv      Salinas et al. (arXiv 2610.24882)
  k2           k2singles/data/k2_singles.csv               LaCourse & Jacobs 2018 events with a TIC
  toi_single   toi_singles_k2/targets.csv                  open single/duo TOIs
  own_ctoi     singles/results/final_classification.csv    rows where we found exactly one period
  own_salinas  salinas-second-transits/results/round3/oddeven_summary.csv

Outputs (gaia_dr4/): candidates_frozen.csv, hosts_frozen.csv, candidates_dropped.csv
Usage: python gaia_dr4/build_hosts.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "gaia_dr4"
COLS = ["tic", "catalogue", "cand_id", "t1_btjd", "t2_btjd", "period_d"]


def _frame(cat, tic, cand_id, t1=None, t2=None, period=None):
    n = len(tic)
    blank = pd.Series([float("nan")] * n)
    return pd.DataFrame(dict(
        tic=pd.to_numeric(pd.Series(list(tic)), errors="coerce"),
        catalogue=cat,
        cand_id=pd.Series(list(cand_id)).astype(str),
        t1_btjd=pd.to_numeric(pd.Series(list(t1)), errors="coerce") if t1 is not None else blank,
        t2_btjd=pd.to_numeric(pd.Series(list(t2)), errors="coerce") if t2 is not None else blank,
        period_d=pd.to_numeric(pd.Series(list(period)), errors="coerce") if period is not None else blank,
    ))


def load_all(root=ROOT):
    parts = []
    h = pd.read_csv(root / "data/catalog/candidates.csv")
    parts.append(_frame("hawthorn", h.TIC, "H" + h.ID.astype(str), h.Tc1, h.Tc2))

    c = pd.read_csv(root / "singles/data/ctoi_singles.csv")
    parts.append(_frame("ctoi", c.TIC, c.CTOI.map(lambda v: f"{v:.2f}"), c.epoch_btjd,
                        period=c.period_listed.where(c.period_listed > 0)))

    s = pd.read_csv(root / "salinas-second-transits/candidates.csv")
    sid = s.tic.astype(str) + "_" + s.planet.fillna("").astype(str) + s.t0.map(lambda v: f"{v:.3f}")
    parts.append(_frame("salinas", s.tic, sid, s.t0, period=s.period_d))

    k = pd.read_csv(root / "k2singles/data/k2_singles.csv")
    kid = "EPIC" + k.EPIC.astype(str) + "_" + k.epoch_btjd.map(lambda v: f"{v:.3f}")
    parts.append(_frame("k2", k.TIC, kid, k.epoch_btjd))

    t = pd.read_csv(root / "toi_singles_k2/targets.csv")
    parts.append(_frame("toi_single", t.TIC, "TOI" + t.TOI.map(lambda v: f"{v:.2f}"), t.epoch_btjd,
                        period=t.period_d.where(t.period_d > 0)))

    f = pd.read_csv(root / "singles/results/final_classification.csv")
    f = f[f.n_candidate_periods == 1]
    parts.append(_frame("own_ctoi", f.TIC, "own" + f.TIC.astype(str), period=f.periods))

    o = pd.read_csv(root / "salinas-second-transits/results/round3/oddeven_summary.csv")
    parts.append(_frame("own_salinas", o.tic, "own" + o.cand_id.astype(str), o.T0_btjd, period=o.P_d))

    return pd.concat(parts, ignore_index=True)[COLS]


def split(df):
    """Return (kept, dropped). A row is dropped, with a reason, if it has no usable TIC, or
    duplicates an earlier (catalogue, cand_id)."""
    reason = pd.Series("", index=df.index)
    reason[df.tic.isna() | (df.tic <= 0)] = "no TIC"
    dup = df.duplicated(["catalogue", "cand_id"]) & (reason == "")
    reason[dup] = "duplicate (catalogue, cand_id)"
    kept = df[reason == ""].copy()
    kept["tic"] = kept.tic.astype("int64")
    dropped = df[reason != ""].assign(reason=reason[reason != ""])
    return kept, dropped


def hosts(kept):
    g = kept.groupby("tic")
    return pd.DataFrame(dict(
        catalogues=g.catalogue.agg(lambda s: ";".join(sorted(set(s)))),
        n_cand=g.size(),
    )).reset_index()


def main():
    kept, dropped = split(load_all())
    kept.to_csv(OUT / "candidates_frozen.csv", index=False)
    dropped.to_csv(OUT / "candidates_dropped.csv", index=False)
    h = hosts(kept)
    h.to_csv(OUT / "hosts_frozen.csv", index=False)
    print(kept.catalogue.value_counts().to_string())
    print(f"candidates kept {len(kept)}, dropped {len(dropped)}, hosts {len(h)}")
    print("hosts in 2+ external catalogues:",
          int(h.catalogues.map(lambda c: len([x for x in c.split(';') if not x.startswith('own')]) > 1).sum()))


if __name__ == "__main__":
    main()
