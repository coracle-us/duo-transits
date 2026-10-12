"""
Parse the third-body solution tables of Borkovits, Mitnyan, Czavalinga & Rappaport 2026
(Universe 12, 107; arXiv:2604.04004; LaTeX source cached in cache/arxiv/2604.04004.bin by
known_triples.py) into data/b26_triples.csv.

Tables (labels in the LaTeX) -> certainty subgroup:
  Tab:Orbelem_LTTE1 -> L1 (certain / very likely pure LTTE)        35 systems
  Tab:Orbelem_LTTE2 -> L2 (less certain but likely)                 28
  Tab:Orbelem_LTTE3 -> L3 (very uncertain, questionable)            64
  Tab:Orbelemdyn1   -> D1 (LTTE + dynamical, P2 < half the data)    31
  Tab:Orbelemdyn2   -> D2 (less certain LTTE + dynamical)           10
L tables: P1, dP1 (1e-10 d/c), P2 (d), a_AB sin i2 (Rsun), e2, w2 (deg), tau2 (MBJD), f(mC),
          (mC)min, Adyn/ALTTE, m_AB.
D tables: P1, P2, a2 (Rsun, relative orbit), e2, w2, tau2, f(mC), mC/mABC, m_AB, mC, Adyn/ALTTE.
Values like '1013.7 (6)' are value (uncertainty in the last digits) -> 1013.7 +- 0.6.
A row whose TIC cell is empty is a second (outer, fourth-body) solution of the row above; it
is kept with solution=2.
A_LTTE (s) for L rows = a sin i / c * sqrt(1 - e^2 cos^2 w), Rsun/c = 2.3206 s, error from
the a sin i error only (e and w errors are small for L1).
"""
import io
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
SRC = HERE / "cache" / "arxiv" / "2604.04004.bin"
OUT = HERE / "data" / "b26_triples.csv"
RSUN_S = 2.3206
GROUPS = {"Tab:Orbelem_LTTE1": "L1", "Tab:Orbelem_LTTE2": "L2", "Tab:Orbelem_LTTE3": "L3",
          "Tab:Orbelemdyn1": "D1", "Tab:Orbelemdyn2": "D2"}
LCOLS = ["P1", "dP1", "P2", "asini_rsun", "e2", "w2", "tau2", "f_m", "mC_min", "Adyn_ALTTE", "m_AB"]
DCOLS = ["P1", "P2", "a2_rsun", "e2", "w2", "tau2", "f_m", "mC_over_mABC", "m_AB", "mC", "Adyn_ALTTE"]


def val_err(cell):
    """'1013.7 (6)' -> (1013.7, 0.6); '58,284 (2)' -> (58284, 2); '$-$' -> (nan, nan)."""
    c = re.sub(r"\\textbf\{|\\boldmath\{|\$|\{|\}|\\,|\*|:", "", cell).replace("−", "-").strip()
    c = c.replace(",", "")
    m = re.match(r"^\s*(-?[\d.]+)(?:e(-?\d+))?\s*(?:\((\d+)\))?", c)
    if not m:
        return np.nan, np.nan
    s = m.group(1)
    v = float(s)
    err = np.nan
    if m.group(3):
        dec = len(s.split(".")[1]) if "." in s else 0
        err = int(m.group(3)) * 10 ** (-dec)
    return v, err


def table_body(t, label):
    i = t.index("\\label{" + label + "}")
    j = t.index("\\midrule", i)
    k = t.index("\\bottomrule", j)
    return t[j + len("\\midrule"):k]


def main():
    with tarfile.open(fileobj=io.BytesIO(SRC.read_bytes())) as tf:
        t = tf.extractfile("universe-4176579-arXiv.tex").read().decode("utf-8", "replace")
    rows = []
    for label, grp in GROUPS.items():
        body = table_body(t, label)
        cols = LCOLS if grp.startswith("L") else DCOLS
        last_tic = None
        for line in body.split("\\\\"):
            line = re.sub(r"%.*", "", line).strip()
            if not line or "&" not in line:
                continue
            cells = [c.strip() for c in line.split("&")]
            m = re.match(r"^(\d{6,10})", cells[0])
            note = re.findall(r"\^\{?([a-z])\}?", cells[0])
            if m:
                tic, sol = int(m.group(1)), 1
                last_tic = tic
            elif last_tic is not None and cells[0] == "":
                tic, sol = last_tic, 2
            else:
                continue
            r = dict(tic=tic, group=grp, solution=sol, notes="".join(note))
            for name, cell in zip(cols, cells[1:]):
                v, e = val_err(cell)
                r[name], r[name + "_err"] = v, e
            rows.append(r)
    df = pd.DataFrame(rows)
    L = df.group.str.startswith("L")
    e = df.e2.fillna(0)
    w = np.radians(df.w2.fillna(0))
    df["A_ltte_s"] = np.where(L, df.asini_rsun * RSUN_S * np.sqrt(1 - e ** 2 * np.cos(w) ** 2), np.nan)
    df["A_ltte_err_s"] = np.where(L, df.asini_rsun_err * RSUN_S, np.nan)
    df.to_csv(OUT, index=False)
    s1 = df[df.solution == 1]
    print(s1.group.value_counts().sort_index().to_dict(), "unique TIC:", s1.tic.nunique(), "-> ", OUT)
    print(s1[s1.group == "L1"][["tic", "P1", "P2", "P2_err", "A_ltte_s", "A_ltte_err_s", "e2", "e2_err"]].head(40).to_string())


if __name__ == "__main__":
    main()
