"""
Rank the Prsa et al. 2022 TESS EBs by expected eclipse count and time span over TESS
Sectors 1-107, and write the top 200.

Per EB: sectors on silicon (tess-point, data/coverage_all.csv); good_days per sector
(data/sector_table.csv, measured from SPOC 2-min files); N_ecl = sum(good_days) / P
(expected primary eclipses); span = t_end(last sector) - t_start(first sector).

Selection rule (fixed before any timing was done):
  - catalog period present, 0.4 d <= P <= 30 d (shorter: contact/ellipsoidal-dominated and
    spot-dominated timing; longer: too few eclipses);
  - morphology coefficient 0 <= morph < 0.6 (detached + semi-detached; Prsa et al. 2022);
  - primary depth (2-gaussian, Dp-2g) >= 0.01, so eclipses are individually timeable;
  - span >= 2000 d (>= 5.5 yr, so P3 of a few years can show at least one cycle);
  - ranked by N_ecl; top 200.
Flags: hemisphere/CVZ (|ecliptic latitude| > 78 deg), known_multiple (data/known_multiples.csv).
Output: data/sample_all.csv (every EB) and data/sample_top200.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic
import astropy.units as u

D = Path(__file__).resolve().parents[1] / "data"


def main():
    cat = pd.read_csv(D / "prsa2022_tess_ebs.csv").drop_duplicates("TIC")
    cov = pd.read_csv(D / "coverage_all.csv").fillna("")
    sec = pd.read_csv(D / "sector_table.csv").set_index("sector")
    known = pd.read_csv(D / "known_multiples.csv").set_index("TIC").papers
    df = cat.merge(cov, on="TIC")
    rows = []
    for r in df.itertuples():
        s = [int(x) for x in r.sectors.split(";") if x] if r.sectors else []
        s = [x for x in s if x in sec.index]
        gd = sec.loc[s, "good_days"].sum() if s else 0.0
        span = (sec.loc[s, "t_end"].max() - sec.loc[s, "t_start"].min()) if s else 0.0
        rows.append(dict(n_sectors=len(s), good_days=gd, span_d=span,
                         first_sector=min(s) if s else np.nan, last_sector=max(s) if s else np.nan))
    df = pd.concat([df.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    df["N_ecl"] = df.good_days / df.Per
    ecl = SkyCoord(df.RAJ2000.values * u.deg, df.DEJ2000.values * u.deg).transform_to(GeocentricTrueEcliptic())
    df["ecl_lat"] = ecl.lat.deg
    df["region"] = np.where(df.ecl_lat > 78, "NCVZ", np.where(df.ecl_lat < -78, "SCVZ",
                            np.where(df.ecl_lat > 0, "north", "south")))
    df["known_multiple"] = df.TIC.map(known).fillna("")
    sel = (df.Per.between(0.4, 30) & (df.Morph >= 0) & (df.Morph < 0.6) & (df["Dp-2g"] >= 0.01)
           & (df.span_d >= 2000))
    df["selected"] = sel
    df.to_csv(D / "sample_all.csv", index=False)
    top = df[sel].sort_values("N_ecl", ascending=False).head(200).reset_index(drop=True)
    top.insert(0, "rank", np.arange(1, len(top) + 1))
    cols = ["rank", "TIC", "RAJ2000", "DEJ2000", "Tmag", "Per", "BJD0", "Morph", "Dp-2g", "Ds-2g",
            "Phis-2g", "Wp-2g", "Ws-2g", "n_sectors", "first_sector", "last_sector", "good_days",
            "span_d", "N_ecl", "ecl_lat", "region", "known_multiple"]
    top[cols].to_csv(D / "sample_top200.csv", index=False, float_format="%.6g")
    print(f"catalog {len(df)}; passing cuts {sel.sum()}")
    print("top200 by region:", top.region.value_counts().to_dict())
    print(f"top200 known multiples: {(top.known_multiple != '').sum()}")
    print(top[cols[:3] + ["Per", "Morph", "n_sectors", "span_d", "N_ecl", "region", "known_multiple"]].head(30).to_string())
    print(top[["n_sectors", "span_d", "N_ecl", "Per", "Tmag"]].describe().T)


if __name__ == "__main__":
    main()
