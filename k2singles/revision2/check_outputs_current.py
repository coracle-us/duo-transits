"""
AAS82225 revision: check that the saved stage-2 outputs (periodic/*.json, thresholds.json) match what the current
code gives, without rerunning the search. Recomputes the thresholds, the common-mode flags (current Poisson rule)
and each star's significant events from per_star/*.json, and compares them with the saved files. Also counts what
the old common-mode rule (>= 2 other stars within 0.5 d) would have flagged, to explain counts written under it.
Usage:  python k2singles/revision2/check_outputs_current.py
Output: stdout only.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import search_periodic as sp  # noqa: E402

DATA = HERE / "data"


def old_common_mode(reps):
    rows = [(r["EPIC"], e["sector"], e["T"]) for r in reps for e in r["events"]]
    df = pd.DataFrame(rows, columns=["EPIC", "sector", "T"])
    flag = set()
    for _, g in df.groupby("sector"):
        for epic, T in zip(g.EPIC, g["T"]):
            if g[(g.EPIC != epic) & (np.abs(g["T"] - T) < sp.CM_WINDOW)].EPIC.nunique() >= 2:
                flag.add((epic, round(T, 5)))
    return flag


def significant(reps, thr, cm):
    out = {}
    for r in reps:
        out[int(r["EPIC"])] = sorted(round(e["T_fit"], 4) for e in r["events"]
                                     if sp.eligible(e) and e["snr"] >= thr[e["mission"]]["S_star"]
                                     and (int(r["EPIC"]), round(e["T"], 5)) not in cm)
    return out


def catalogued_significant(reps, sig):
    n = 0
    for r in reps:
        T = sig[int(r["EPIC"])]
        n += all(any(abs(t - e["epoch"]) < max(0.1, (e["blind_D_h"] or 0) / 48) for t in T) for e in r["catalogued_events"])
    return n


def main():
    reps = [json.loads(f.read_text()) for f in sorted((DATA / "per_star").glob("*.json"))]
    reps = [r for r in reps if r.get("status") == "ok"]
    saved_thr = json.loads((DATA / "thresholds.json").read_text())
    thr = sp.thresholds(reps)
    for m in ("K2", "TESS"):
        print(f"{m}: S* saved {saved_thr[m]['S_star']}, recomputed {thr[m]['S_star']}")
    cm_new, cm_old = sp.common_mode(reps), old_common_mode(reps)
    print(f"common-mode events: saved {saved_thr['n_common_mode_events']}, current rule {len(cm_new)}, old rule {len(cm_old)}")
    sig_new, sig_old = significant(reps, thr, cm_new), significant(reps, thr, cm_old)
    mismatch = []
    for r in reps:
        p = json.loads((DATA / "periodic" / f"{r['EPIC']}.json").read_text())
        if sorted(s["T"] for s in p["significant"]) != sig_new[int(r["EPIC"])]:
            mismatch.append(int(r["EPIC"]))
    print(f"stars whose saved significant events differ from the current rule: {len(mismatch)} {mismatch[:10]}")
    print(f"catalogued events recovered at a significant level: current rule {catalogued_significant(reps, sig_new)}, "
          f"old rule {catalogued_significant(reps, sig_old)}")


if __name__ == "__main__":
    main()
