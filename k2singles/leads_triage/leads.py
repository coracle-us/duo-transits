"""
The leads to triage: every signal on the 41 stars of k2singles/revision2/results/new_leads.csv (pre-registered blind
search: linked families of >= 2 significant events, and BLS detections) plus the post-hoc linked-period detections not in
the submitted Table 1 (results/link_calibrated.csv; 8 stars, 2 of them also in new_leads). 47 stars in all.

A "signal" is one row of new_leads.csv, or the top-support period of a post-hoc detection. Each carries the event
times it rests on and, where it has exactly one, a period.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REV = HERE.parent / "revision2"
DATA = REV / "data"
POSTHOC = ["202085278", "204776782", "247795097", "215067200", "205966706", "247756662", "211598816", "212805678"]


def lead_epics():
    L = pd.read_csv(REV / "results" / "new_leads.csv", dtype={"EPIC": str})
    return sorted(set(L.EPIC) | set(POSTHOC))


def signals():
    """List of dicts: EPIC, TIC, sig_id, kind, P (or None), T0, T_events (times the signal rests on), D (d),
    depth (fraction), stat (support or SDE), block (t0, t1) for BLS."""
    L = pd.read_csv(REV / "results" / "new_leads.csv", dtype={"EPIC": str})
    lc = json.loads((DATA / "link_calibrated.json").read_text())
    lcs = {str(s["EPIC"]): s for s in lc["stars"]}
    out = []
    for e in lead_epics():
        per = json.loads((DATA / "periodic" / f"{e}.json").read_text())
        tic = int(per["TIC"])
        k = 0
        for fm in per["families"]:
            if fm["n_members"] < 2:  # as compare.py
                continue
            Ps = [q["P"] for q in fm["periods"]]
            out.append(dict(EPIC=e, TIC=tic, sig_id=f"{e}-F{k}", kind="family", n_allowed=fm["n_allowed"],
                            P=Ps[0] if len(Ps) == 1 else None, P_list=Ps, T_events=list(fm["T"]),
                            T0=min(fm["T"]), D=fm["D_h"] / 24, depth=fm["depth_ppm"] * 1e-6,
                            stat=fm["n_members"], stat_name="members"))
            k += 1
        j = 0
        for b in per["bls"]:
            for s in b["signals"]:
                if not s.get("detected"):
                    continue
                n0, n1 = np.ceil((b["t0"] - s["T0"]) / s["P"]), np.floor((b["t1"] - s["T0"]) / s["P"])
                Tp = list(s["T0"] + np.arange(n0, n1 + 1) * s["P"])
                out.append(dict(EPIC=e, TIC=tic, sig_id=f"{e}-B{j}", kind="BLS", n_allowed=1, P=s["P"], P_list=[s["P"]],
                                T_events=Tp, T0=s["T0"], D=s["D_h"] / 24, depth=s["depth_ppm"] * 1e-6, stat=s["sde"],
                                stat_name="SDE", block=(b["t0"], b["t1"])))
                j += 1
        if e in POSTHOC:
            st = lcs[e]
            best = st["real"][0]  # sorted by (-support, -P)
            top = [p for p in st["real"] if p["support"] == best["support"]]
            out.append(dict(EPIC=e, TIC=tic, sig_id=f"{e}-H", kind="posthoc", n_allowed=len(top), P=best["P"],
                            P_list=[p["P"] for p in top], T_events=list(best["T"]), T0=min(best["T"]),
                            D=best["D_h"] / 24, depth=best["depth_ppm"] * 1e-6, stat=st["max_support"],
                            stat_name="support", stat_control=st["max_support_control"]))
    return out
