"""
Automatic acceptance rule + comparison with known periods.

Rule (set on the validation set, data/validation_results.csv, and then frozen for the pilot):
  auto_pass = boot_fap <= 0.005  (observed LS peak above all 200 block-shuffled maxima)
              & frac_contam < 0.4 (fewer than 40% of good frames star-contaminated)
              & n_cycles_2P >= 2  (at least two rotations of the adopted period in the data span)
Suggested LCDB-style code for auto_pass objects (the final code is decided by eye, in
data/<set>_manual.csv):
  3- : odd-harmonic test significant (p_odd_red < 1e-3), >= 5 cycles, amp_snr >= 20
  2  : otherwise (typically P vs 2P not settled by the data, or few cycles)
'match' (when a known period is given): |P2 - P_known| / P_known < 1% -> 'P'; within 1% of
P_known/2 or 2 P_known -> 'alias x0.5 / x2'; else 'wrong'.
Usage: python classify.py RESULTS.csv [KNOWN.csv with number,P_known_h]
"""
import sys

import numpy as np
import pandas as pd


def classify(d):
    d = d.copy()
    d["auto_pass"] = (d.boot_fap <= 0.005) & (d.frac_contam < 0.4) & (d.n_cycles_2P >= 2)
    d["U_suggest"] = np.where(~d.auto_pass, "", np.where(
        (d.p_odd_red < 1e-3) & (d.n_cycles_2P >= 5) & (d.amp_snr >= 20), "3-", "2"))
    if "P_known_h" in d.columns:
        r = d.P2_h / d.P_known_h
        d["match"] = np.select([np.abs(r - 1) < 0.01, np.abs(r - 0.5) < 0.005, np.abs(r - 2) < 0.02],
                               ["P", "alias x0.5", "alias x2"], "wrong")
        d.loc[d.P_known_h.isna(), "match"] = ""
        d["dP_pct"] = (100 * (r - 1)).round(3)
    return d


if __name__ == "__main__":
    d = pd.read_csv(sys.argv[1])
    if len(sys.argv) > 2:
        k = pd.read_csv(sys.argv[2])[["number", "P_known_h"]]
        d = d.drop(columns=[c for c in ["P_known_h"] if c in d.columns]).merge(
            k.rename(columns={"number": "id"}), on="id", how="left")
    c = classify(d)
    out = sys.argv[1].replace(".csv", "_classified.csv")
    c.to_csv(out, index=False)
    cols = [x for x in ["id", "label", "set", "V_med", "P2_h", "sP2_h", "amp2", "P_known_h", "dP_pct", "match",
                        "boot_fap", "frac_contam", "n_cycles_2P", "p_odd_red", "auto_pass", "U_suggest"] if x in c]
    print(c[cols].to_string(index=False))
    if "match" in c:
        a = c[c.auto_pass]
        print(f"\nauto_pass {len(a)}/{len(c)}; of these match P: {(a.match == 'P').sum()}, "
              f"alias: {a.match.str.startswith('alias').sum()}, wrong: {(a.match == 'wrong').sum()}")
        f = c[~c.auto_pass]
        print(f"rejected {len(f)}; of these match P anyway: {(f.match == 'P').sum()}")
