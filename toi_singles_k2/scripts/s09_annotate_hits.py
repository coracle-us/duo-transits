"""
Annotate k2_hits.csv: does each K2 event coincide with a transit of an already-catalogued signal on the
same star? Checks every TOI on the TIC (ExoFOP TOI table, period > 0; for a TOI with a period, all
epochs T0 + kP; the TOI epoch itself is also checked) within max(1 h, 3 sigma_T) where sigma_T
propagates the TOI period/epoch errors. Also adds a vetting verdict:
  - 'systematic'     found in <2 reductions or depth inconsistent between reductions (>2x spread)
  - 'known: <TOI>'   coincides with a catalogued signal
  - 'candidate'      otherwise
Writes k2_hits.csv in place (adds columns attributed_to, vet).
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "toi_singles_k2"


def main():
    h = pd.read_csv(OUT / "k2_hits.csv")
    toi = pd.read_csv(ROOT / "data" / "raw" / "exofop_toi.csv")
    att, vet = [], []
    for r in h.itertuples():
        lab = []
        for t in toi[toi["TIC ID"] == r.TIC].itertuples():
            T0 = t._32 - 2457000
            P, sP, sT0 = t._34, (t._35 if t._35 == t._35 else 0), (t._33 if t._33 == t._33 else 0)
            if P and P > 0:
                k = np.round((r.T_btjd - T0) / P)
                pred = T0 + k * P
                tol = max(1 / 24, 3 * np.hypot(sT0, abs(k) * sP))
            else:
                pred, tol = T0, 1 / 24
            if abs(pred - r.T_btjd) < max(tol, r.toi_dur_h / 24 / 2):
                lab.append(f"TOI-{t.TOI} (pred {pred:.3f})")
        if r.TOI == 5560.01:  # control: both K2 events were found and published by this repo's k2singles project
            lab.append("this repo k2singles/ (EPIC 211924561 C5+C16 linked to TOI-5560.01; RNAAS submitted)")
        att.append("; ".join(lab))
        deps = [getattr(r, f"{p}_depth_ppm") for p in ["K2PDC", "K2SAP", "EVEREST", "K2SFF"]
                if getattr(r, f"{p}_verdict") == "supported"]
        if r.n_found < 2 or len(deps) < 2 or max(deps) > 2 * min(deps):
            vet.append("systematic/inconsistent")
        elif lab:
            vet.append("known: " + ("TOI-5560.01 (control, own prior work)" if r.TOI == 5560.01 else lab[0].split(" ")[0]))
        else:
            vet.append("candidate")
    h["attributed_to"] = att
    h["vet"] = vet
    h.to_csv(OUT / "k2_hits.csv", index=False)
    print(h[["TOI", "EPIC", "T_btjd", "found_in", "attributed_to", "vet"]].to_string(index=False))


if __name__ == "__main__":
    main()
