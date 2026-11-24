"""
Trial definitions for the catalogue run, from catalogues and the TIC only (no light curve is read).

One row per scored row of gaia_dr4/candidates_frozen.csv (prereg_catalogue.yaml, scope.scored):
  row_id, tic, catalogue, cand_id, t_cat (first catalogued transit, BTJD), t2_cat (Hawthorn only),
  dur_cat_d, depth_cat (fraction), depth_cat_err, mstar, rstar, tmag, multi, stated_period_d,
  flags (t16_coincident, t16_gaia_id_disagrees, depth_cat_out_of_range), n_rows_on_star,
  n_toi_on_star, n_ctoi_on_star
Rows not scored (k2, own_ctoi, own_salinas) go to defs_not_scored.csv with the reason.

"multi": more than one scored row on the star, or more than one TOI, or more than one CTOI, in the
ExoFOP tables of 2026-11-22. A star with one TOI and one CTOI is not multi on that ground: the two
are usually the same signal (counting them together would make 96 of the 156 TOI rows multi).
Flag depth_cat_out_of_range: catalogued depth under 50 ppm or over 25%; the row is run as it is (data/raw/exofop_toi.csv, exofop_ctoi.csv as downloaded for gaia_dr4/t16_overlap.py).

Inputs  : gaia_dr4/candidates_frozen.csv and the source catalogues named in gaia_dr4/build_hosts.py;
          the TIC (MAST, network) for mass, radius and Tmag, cached in catalogue/tic_info.csv
Outputs : catalogue/defs.csv, catalogue/defs_not_scored.csv, catalogue/tic_info.csv
Usage   : python catalogue/build_defs.py [--toi path --ctoi path]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "catalogue"
SCORED = ("hawthorn", "ctoi", "salinas", "toi_single", "t16")


def sources(root=ROOT):
    """cand_id -> catalogued duration (d), depth (fraction), depth error, per catalogue; ids built
    exactly as gaia_dr4/build_hosts.py builds them."""
    parts = []
    h = pd.read_csv(root / "data/catalog/candidates.csv")
    parts.append(pd.DataFrame(dict(catalogue="hawthorn", cand_id="H" + h.ID.astype(str),
                                   dur_cat_d=h.Dur_h / 24, depth_cat=h.Depth1_ppt * 1e-3, depth_cat_err=np.nan)))
    c = pd.read_csv(root / "singles/data/ctoi_singles.csv")
    parts.append(pd.DataFrame(dict(catalogue="ctoi", cand_id=c.CTOI.map(lambda v: f"{v:.2f}"),
                                   dur_cat_d=c.dur_h / 24, depth_cat=c.depth_ppm * 1e-6, depth_cat_err=np.nan)))
    s = pd.read_csv(root / "salinas-second-transits/candidates.csv")
    sid = s.tic.astype(str) + "_" + s.planet.fillna("").astype(str) + s.t0.map(lambda v: f"{v:.3f}")
    parts.append(pd.DataFrame(dict(catalogue="salinas", cand_id=sid, dur_cat_d=s.duration_h / 24,
                                   depth_cat=s.depth_ppm * 1e-6, depth_cat_err=s.depth_err * 1e-6)))
    t = pd.read_csv(root / "toi_singles_k2/targets.csv")
    parts.append(pd.DataFrame(dict(catalogue="toi_single", cand_id="TOI" + t.TOI.map(lambda v: f"{v:.2f}"),
                                   dur_cat_d=t.dur_h / 24, depth_cat=t.depth_ppm * 1e-6, depth_cat_err=np.nan)))
    x = pd.read_csv(root / "gaia_dr4/data_t16/t16_overlap.csv")
    x = x[x.drop_reason.isna()]
    parts.append(pd.DataFrame(dict(catalogue="t16", cand_id="T16_" + x.tic.astype(str), dur_cat_d=x.dur_d,
                                   depth_cat=1 - 10 ** (-0.4 * x.depth_mag), depth_cat_err=np.nan)))
    return pd.concat(parts, ignore_index=True)


def t16_flags(root=ROOT):
    x = pd.read_csv(root / "gaia_dr4/data_t16/t16_overlap.csv")
    dis = set(pd.read_csv(root / "gaia_dr4/data_t16/t16_gaia_id_disagreements.csv").tic)
    f = {}
    for r in x.itertuples():
        fl = (["t16_coincident"] if r.n_coincident > 0 else []) + (["t16_gaia_id_disagrees"] if r.tic in dis else [])
        f[int(r.tic)] = ";".join(fl)
    return f


def tic_table(tics, path):
    """Mass, radius, Tmag from the TIC; cached. Queried in chunks of 200 IDs."""
    have = pd.read_csv(path) if path.exists() else pd.DataFrame(columns=["tic", "mstar", "rstar", "tmag", "teff"])
    need = sorted(set(map(int, tics)) - set(have.tic.astype(int)))
    if need:
        from astroquery.mast import Catalogs
        rows = []
        for i in range(0, len(need), 200):
            r = Catalogs.query_criteria(catalog="Tic", ID=need[i:i + 200]).to_pandas()
            rows.append(pd.DataFrame(dict(tic=r.ID.astype(int), mstar=pd.to_numeric(r["mass"], errors="coerce"),
                                          rstar=pd.to_numeric(r["rad"], errors="coerce"),
                                          tmag=pd.to_numeric(r["Tmag"], errors="coerce"),
                                          teff=pd.to_numeric(r["Teff"], errors="coerce"))))
        have = pd.concat([have] + rows, ignore_index=True).drop_duplicates("tic").sort_values("tic")
        have.to_csv(path, index=False)
    return have


def n_toi_ctoi(toi_path, ctoi_path):
    def count(p):
        d = pd.read_csv(p, comment="#", low_memory=False)
        col = next(c for c in d.columns if c.strip().lower() in ("tic id", "tic", "tic_id"))
        return pd.to_numeric(d[col], errors="coerce").dropna().astype(int).value_counts()
    return count(toi_path), count(ctoi_path)


def build(toi_path, ctoi_path, root=ROOT, tic_path=OUT / "tic_info.csv"):
    c = pd.read_csv(root / "gaia_dr4/candidates_frozen.csv")
    reason = {"k2": "reference transit is in K2 data, not TESS",
              "own_ctoi": "a period this project derived, not a candidate",
              "own_salinas": "a period this project derived, not a candidate"}
    not_scored = c[~c.catalogue.isin(SCORED)].assign(reason=lambda d: d.catalogue.map(reason))
    d = c[c.catalogue.isin(SCORED)].merge(sources(root), on=["catalogue", "cand_id"], how="left", validate="1:1")
    assert d[["dur_cat_d", "depth_cat"]].notna().all().all(), "a scored row lacks depth or duration"
    d = d.merge(tic_table(d.tic, Path(tic_path)), on="tic", how="left")
    nrow = d.groupby("tic").tic.transform("size")
    ct, cc = n_toi_ctoi(toi_path, ctoi_path)
    ntoi, nctoi = d.tic.map(ct).fillna(0).astype(int), d.tic.map(cc).fillna(0).astype(int)
    odd = (d.depth_cat < 50e-6) | (d.depth_cat > 0.25)
    fl = t16_flags(root)
    out = pd.DataFrame(dict(
        row_id=d.catalogue + ":" + d.cand_id, tic=d.tic.astype(int), catalogue=d.catalogue, cand_id=d.cand_id,
        t_cat=d.t1_btjd, t2_cat=d.t2_btjd, dur_cat_d=d.dur_cat_d, depth_cat=d.depth_cat,
        depth_cat_err=d.depth_cat_err, mstar=d.mstar, rstar=d.rstar, tmag=d.tmag,
        multi=(nrow > 1) | (ntoi > 1) | (nctoi > 1), stated_period_d=d.period_d,
        flags=[";".join(x for x in ((fl.get(int(t), "") if cat == "t16" else ""),
                                    "depth_cat_out_of_range" if o else "") if x)
               for t, cat, o in zip(d.tic, d.catalogue, odd)],
        n_rows_on_star=nrow, n_toi_on_star=ntoi, n_ctoi_on_star=nctoi))
    assert out.row_id.is_unique
    return out, not_scored


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--toi", default=str(ROOT / "data/raw/exofop_toi.csv"))
    ap.add_argument("--ctoi", default=str(ROOT / "data/raw/exofop_ctoi.csv"))
    a = ap.parse_args()
    out, ns = build(a.toi, a.ctoi)
    out.to_csv(OUT / "defs.csv", index=False)
    ns.to_csv(OUT / "defs_not_scored.csv", index=False)
    print(len(out), "scored rows on", out.tic.nunique(), "stars;", len(ns), "not scored;",
          int(out.multi.sum()), "rows multi;", int(out.t_cat.isna().sum()), "without a time")
