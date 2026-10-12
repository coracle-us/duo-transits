"""
Validation: compare our O-C solutions with published third-body orbits.

Published values, transcribed from the arXiv LaTeX sources (cached in cache/arxiv/):
  M24 = Marcadon & Prsa 2024 (arXiv 2403.07694): 2-min data S14-60, LTTE + linear
        ephemeris; A_LTTE = a sin i / c * sqrt(1 - e^2 cos^2 w).
  B26 = Borkovits, Mitnyan, Czavalinga & Rappaport 2026 (arXiv 2604.04004): FFI S14-86 (the
        same baseline as ours for the NCVZ). Table columns P2, a_AB sin i2 (Rsun), e2, w2;
        A computed here from a sin i with the same formula (Rsun/c = 2.3206 s). For
        the dynamical systems (259006185, 356896561) only P2 and e2 are comparable.
        B26 lists 229771234 as TIC 229771231 (M24 notes the two TIC IDs are the same star).
  B20 = Borkovits et al. 2020 (arXiv 2002.06902): photodynamical; dynamically dominated
        ETVs, only P_out and e_out comparable.
Output: sidequests/eb_timing/data/validation_table.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
RSUN_S = 2.3206


def A_from(asini_rsun, e, w_deg):
    return asini_rsun * RSUN_S * np.sqrt(1 - e ** 2 * np.cos(np.radians(w_deg)) ** 2)


PUB = [
    # tic, ref, P3, A_s, e3, note
    (219900027, "M24", 1159.4, 113.0, 0.506, "LTTE"),
    (219900027, "B26", 1631.0, A_from(58.8, 0.49, 66), 0.49, "B26: 'very uncertain, questionable'"),
    (229771234, "M24", 359.8, 179.7, 0.341, "LTTE"),
    (229771234, "B26", 359.1, A_from(85.5, 0.50, 146), 0.50, "B26 'certain'; listed as TIC 229771231"),
    (259006185, "M24", 199.6, 38.6, 0.079, "LTTE+dyn (A_dyn/A_LTTE = 0.6)"),
    (259006185, "B26", 200.0, None, 0.12, "dynamical+LTTE; apsidal motion"),
    (356896561, "M24", 477.5, 51.2, 0.005, "GZ Dra"),
    (356896561, "B26", 475.3, None, 0.52, "dynamical+LTTE (A_dyn/A_LTTE = 0.34)"),
    (424461577, "M24", 995.6, 316.5, 0.175, "LTTE"),
    (424461577, "B26", 1013.7, A_from(205.6, 0.480, 121.3), 0.480, "B26 'less certain but likely'"),
    (167692429, "B20", 331.5, None, 0.557, "photodynamical; dynamical ETV; eclipse depth variations"),
    (220397947, "B20", 77.08, None, 0.225, "photodynamical; dynamical ETV (2+1+1)"),
]


def main():
    rows = []
    for tic, ref, P3, A, e3, note in PUB:
        f = HERE / "results" / "oc" / f"{tic}_fit.json"
        if not f.exists():
            continue
        r = json.loads(f.read_text())
        L = r["ltte"]
        rows.append(dict(tic=tic, ref=ref, P3_pub=P3, P3_ours=round(L["P3"], 1),
                         P3_err_ours=round(L["P3_err"], 1), P3_ratio=round(L["P3"] / P3, 3),
                         A_pub_s=round(A, 1) if A else None, A_ours_s=round(L["A_ltte_s"], 1),
                         A_ratio=round(L["A_ltte_s"] / A, 3) if A else None,
                         e_pub=e3, e_ours=round(L["e3"], 3), our_base=r["ltte_base"],
                         P3_lin=round(r["ltte_variants"]["lin"]["P3"], 1),
                         A_lin_s=round(r["ltte_variants"]["lin"]["A_s"], 1),
                         P3_quad=round(r["ltte_variants"]["quad"]["P3"], 1),
                         A_quad_s=round(r["ltte_variants"]["quad"]["A_s"], 1),
                         quad_at_bound=r["ltte_variants"]["quad"]["at_bound"],
                         dBICc_ltte_vs_quad=round(r["dBICc_ltte_vs_quad"], 1),
                         rho_pri_sec=round(r["rho_pri_sec"], 2), span_d=round(r["span_d"]),
                         n_pri=r["n_pri"], n_sec=r["n_sec"], sig_pri_s=round(r["median_sig_pri_s"], 1),
                         note=note))
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "data" / "validation_table.csv", index=False)
    print(df.drop(columns=["note"]).to_string())


if __name__ == "__main__":
    main()
