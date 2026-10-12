"""Write results.csv from data/final_details.json plus the vetted per-star conclusions."""
import os, json, sys
sys.path.insert(0, os.path.dirname(__file__))
import pandas as pd, numpy as np, lib
from targets import T
D = json.load(open(os.path.join(lib.ROOT, "data", "final_details.json")))
EXTRA_EVENTS = {  # supporting events measured outside the automatic pipeline (see REPORT)
    202071902: "3273.0023+-0.0003 S72 partial ingress (TARS + own raw TESScut, scattered-light cadences)",
    202072917: "3261.1888+-0.0003 S72 sector start (own raw TESScut, scattered-light cadences; depth 19%)",
    246089124: "SECONDARY 810.986+-0.001 K2 C12 depth ~6400 ppm, T14~8 h, phase 0.494",
}
CONCL = {
 248811085: ("49.25012 +- 0.00002", "unique", "only if twin eclipses: depths K2 89.8k vs TESS 89.0-89.3k ppm (equal); no dip at phase 0.5 (K2 -57+-26 ppm)", "high"),
 246696483: ("182.2452 +- 0.0003", "unique (P/2 excluded by one S98 epoch only)", "only if twin eclipses (K2 epoch even, TESS odd); depths K2 5.7% vs QLP 5.4-5.6% (equal)", "medium-high"),
 247967714: ("185.62668 | 92.81334 | 61.87556 (each +-0.00003)", "3 aliases", "any of 2x the aliases if twin eclipses; all events 24-26% deep, V-shaped -> probably EB", "high (events); period 3-fold ambiguous"),
 246089124: ("73.96241 +- 0.00002", "unique", "no: secondary eclipse found at phase 0.494 (6400 ppm) -> EB, P_orb = 73.962", "high"),
 211892898: ("72.5100 +- 0.0014", "unique", "only if twin eclipses (C16, C18 odd cycles): C5 0.69% vs C16/C18 0.61-0.62% (similar)", "high"),
 211503363: ("122.9049 +- 0.0007 | 61.4525 +- 0.0003", "2 aliases", "possible if twin eclipses", "medium (TESS events ~10 sigma each)"),
 211390677: ("184.3103 +- 0.0010 | 92.1552 +- 0.0005", "2 aliases", "possible if twin eclipses", "high (events); period 2-fold ambiguous"),
 211411112: ("62.59479 +- 0.00005", "unique", "only if twin eclipses (C16 odd): depths 2.6-3.0% all similar; no dip at phase 0.5", "high"),
 211489484: ("82.25843 +- 0.00007", "unique", "no: distinct secondary eclipses at phase 0.194 (eccentric EB)", "high"),
 202071902: ("22.904130 +- 0.000002", "unique", "only if twin eclipses: S71 (E=151) and S72 partial (E=152) equal depth ~9%; no dip at phase 0.5", "high"),
 202126877: ("21.77521 (mean; ETVs +-0.017 d)", "unique", "no: distinct secondary at phase 0.33 (EB with growing eclipse depths)", "high"),
 202072917: ("34.745606 +- 0.000007", "unique", "only if twin eclipses: consecutive S43/S45 events 23.7% vs 23.8% (identical)", "high"),
 230887315: ("> 70.0 d (windows in data/allowed_windows_230887315.csv)", "single event", "n/a", "high that no further event is in the data"),
 212152316: ("> 68.5 d (windows in data/allowed_windows_212152316.csv)", "single event", "n/a", "high that no further event is in the data"),
}
rows = []
for _, r in T.iterrows():
    e = int(r.epic); v = D[str(e)]
    evs = []
    for x in v["events"]:
        if "t0" not in x: continue
        if e == 202071902 and x["label"] == "p": continue  # replaced by the TARS/raw-own ingress fit in EXTRA_EVENTS
        lab = {"P": "", "S": "SECONDARY ", "p": "partial "}[x["label"]]
        evs.append(f"{lab}{x['t0']:.4f}+-{x['e_t0']:.4f} {x['best_src']} depth {x['depth']*1e6:.0f}ppm T14 {x['T14_h']:.1f}h")
    if e in EXTRA_EVENTS: evs.append(EXTRA_EVENTS[e])
    P, nal, twoP, conf = CONCL[e]
    nprim = sum(1 for x in v["events"] if x.get("label") == "P")
    rows.append(dict(epic=e, tic=int(r.tic), k2_campaign_listed=int(r.camp), cat_epoch=r.t0, cat_depth_ppm=r.depth, cat_dur_h=r.dur,
                     n_events_found=len(evs), n_primary_events=nprim, allowed_periods_d=P, period_status=nal,
                     could_be_2x=twoP, confidence=conf, events=" ; ".join(evs)))
pd.DataFrame(rows).to_csv(os.path.join(lib.ROOT, "results.csv"), index=False)
print(pd.DataFrame(rows)[["epic", "n_events_found", "allowed_periods_d", "confidence"]].to_string())
