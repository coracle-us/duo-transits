"""
AAS82225 revision: recompute the counts quoted in the revised note (k2singles/rnaas/rnaas.tex) from saved outputs.
Output: k2singles/revision2/results/note_numbers.txt
Usage:  python -I k2singles/revision2/note_numbers.py
"""
import glob
import json
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "results")
D = os.path.join(HERE, "data")
out = []


def say(s):
    out.append(s)
    print(s)


cat = pd.read_csv(os.path.join(HERE, "..", "data", "k2_singles.csv"))
say(f"catalogue: {len(cat)} events on {cat.EPIC.nunique()} stars")
c = pd.read_csv(os.path.join(R, "compare.csv"))
miss = c[~c.catalogued_recovered]
say(f"stars with every catalogued event recovered: {c.catalogued_recovered.sum()} of {len(c)}; "
    f"not: {list(miss.EPIC)} (catalogued events on it: {int((cat.EPIC.isin(miss.EPIC)).sum())}, one recovered per PLAN.md)")
say(f"stars with catalogued events significant: {int(c.catalogued_significant.sum())}")
th = json.load(open(os.path.join(R, "thresholds.json")))
say(f"thresholds: {json.dumps(th)[:300]}")
t1 = pd.read_csv(os.path.join(R, "compare_table1.csv"))
rec = t1.groupby("EPIC").apply(lambda g: bool(g.blind_allowed.any() or g.blind_bls.any()), include_groups=False)
say(f"Table 1 systems: {len(rec)}; blind period recovered: {int(rec.sum())}; not: {sorted(rec.index[~rec])}")
lk = pd.read_csv(os.path.join(R, "link_calibrated.csv"))
say(f"post-hoc test: {int(lk.detected.sum())} stars detected; detected among Table 1 misses: "
    f"{sorted(set(lk.EPIC[lk.detected]) & set(rec.index[~rec]))}")
nl = pd.read_csv(os.path.join(R, "new_leads.csv"))
ph = set(lk.EPIC[lk.detected]) - set(rec.index)
say(f"new leads: {nl.EPIC.nunique()} stars ({len(nl)} rows); post-hoc detections not in Table 1: {len(ph)}; "
    f"union: {len(set(nl.EPIC) | ph)} stars")
nb = nbc = 0
for f in glob.glob(os.path.join(D, "periodic", "*.json")):
    d = json.load(open(f))
    for b in d["bls_control"]:
        nb += 1
        nbc += bool(b["signals"]) and b["signals"][0]["sde"] >= d["bls_sde_threshold"]
say(f"BLS blocks: {nb}; inverted blocks at or above the threshold: {nbc}")
lit = os.path.join(R, "literature")
a = pd.read_csv(os.path.join(lit, "table1_classes_final.csv"))
b = pd.read_csv(os.path.join(lit, "table1_vizier_classes_final.csv"))
say(f"Table 1 sources classified: ADS/SIMBAD {len(a)} + VizieR {len(b)} = {len(a) + len(b)}")
vc = pd.read_csv(os.path.join(lit, "vizier_catalogues.csv"))
say(f"VizieR catalogue IDs: {len(vc)}; exist on the mirror: {int(vc.exists.astype(bool).sum())}")
open(os.path.join(R, "note_numbers.txt"), "w").write("\n".join(out) + "\n")
