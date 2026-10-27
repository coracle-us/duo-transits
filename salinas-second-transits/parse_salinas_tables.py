#!/usr/bin/env python3
"""Parse Tables 4-6 of Salinas et al. (arXiv 2610.24882) from the saved full text
into candidates.csv. T0 is given in the paper as 'BJD' but the values are BTJD."""
import re, sys
import pandas as pd

src = sys.argv[1] if len(sys.argv) > 1 else "salinas_2610.24882_fulltext.txt"
s = open(src, encoding="utf-8").read()
i4, i5, i6 = (s.find(f"Table {n}:") for n in (4, 5, 6))
num = r"(-?\d+\.?\d*)"
pm = rf"{num}\s*±\s*{num}"
single = re.compile(rf"^(\d{{5,10}})\s+{pm}\s+{pm}\s+{pm}\s+(-?\d+\.?\d*|—)\s*$", re.M)
multi = re.compile(rf"^(\d{{5,10}})\s+{num}\s+{pm}\s+{pm}\s+{pm}\s+(-?\d+\.?\d*|—)\s*$", re.M)
mplan = re.compile(rf"^(\d{{5,10}})\s*\((\w)\)\s+{pm}\s+{pm}\s+{pm}\s+(-?\d+\.?\d*|—)\s*$", re.M)
rows = []
for m in single.finditer(s[i4:i5]):
    tic, t0, et0, d, ed, dep, edep, rp = m.groups()
    rows.append(dict(tic=int(tic), type="single", planet="", period_d=None, t0=float(t0), t0_err=float(et0),
                     duration_h=float(d), duration_err=float(ed), depth_ppm=float(dep), depth_err=float(edep), rp_rj=(float(rp) if rp != "—" else None)))
for m in multi.finditer(s[i5:i6]):
    tic, p, t0, et0, d, ed, dep, edep, rp = m.groups()
    rows.append(dict(tic=int(tic), type="multi", planet="", period_d=float(p), t0=float(t0), t0_err=float(et0),
                     duration_h=float(d), duration_err=float(ed), depth_ppm=float(dep), depth_err=float(edep), rp_rj=(float(rp) if rp != "—" else None)))
end6 = s.find("\n", s.find("Table 6:") + 3000) if False else len(s)
for m in mplan.finditer(s[i6:]):
    tic, pl, t0, et0, d, ed, dep, edep, rp = m.groups()
    rows.append(dict(tic=int(tic), type="multiplanet", planet=pl, period_d=None, t0=float(t0), t0_err=float(et0),
                     duration_h=float(d), duration_err=float(ed), depth_ppm=float(dep), depth_err=float(edep), rp_rj=(float(rp) if rp != "—" else None)))
df = pd.DataFrame(rows)
df.to_csv("candidates.csv", index=False)
print(df["type"].value_counts().to_string(), "\nrows:", len(df), "unique TICs:", df.tic.nunique())
