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

# real BLS blocks above the threshold (against the inverted count above)
nr = 0
for f in glob.glob(os.path.join(D, "periodic", "*.json")):
    d = json.load(open(f))
    nr += sum(any(s.get("detected") for s in b["signals"]) for b in d["bls"])
say(f"real BLS blocks with a detection: {nr}")

# TESS coverage of the blind search: observed target-sectors (catalogue) vs sectors searched (per_star json)
obs = got = 0
for r in cat.drop_duplicates("EPIC").itertuples():
    o = {int(x) for x in str(r.tess_sectors_observed).split() if x.isdigit()} if pd.notna(r.tess_sectors_observed) else set()
    g = set(json.load(open(os.path.join(D, "per_star", f"{r.EPIC}.json"))).get("tess_sectors") or [])
    obs += len(o)
    got += len(o & g)
say(f"TESS target-sectors observed: {obs}; searched: {got}; not searched (no pipeline light curve, sector < 80): {obs - got}")

# other catalogue periods of the 'agrees' C0 binaries against consecutive observed eclipses
# (eclipse times: independent_table/scripts/s7_final.py, results/ephemerides/202126877.json)
SPACING = {202072917: 2531.533 - 2496.787, 202126877: 2510.368 - 2488.583, 202071902: 3273.002 - 3250.097}  # 202071902: S71 eclipse to the S72 ingress (consecutive)
OTHER = {202072917: [3.842, 23.255, 34.742, 1007.6, 15.17, 14.20, 5.579], 202126877: [3.877], 202071902: [8.6973]}
for e, ps in OTHER.items():
    for p in ps:
        n = SPACING[e] / p
        say(f"EPIC {e}: catalogue period {p} d; observed spacing {SPACING[e]:.3f} d = {n:.3f} cycles; "
            f"offset from an integer {abs(n - round(n)) * p:.3f} d")
open(os.path.join(R, "note_numbers.txt"), "w").write("\n".join(out) + "\n")

# TOI-5560.01 predicted transits shared by both periods (results/ephemerides/211924561.json)
e = json.load(open(os.path.join(HERE, "..", "results", "ephemerides", "211924561.json")))
for n in (30, 31):
    t_ = e["T0"] + n * e["P"]
    sig = (e["sigma_T0"] ** 2 + (n * e["sigma_P"]) ** 2) ** 0.5
    say(f"TOI-5560.01 cycle {n} (140.49 d): BJD_TDB {t_ + 2457000:.3f} +- {sig:.4f}")
open(os.path.join(R, "note_numbers.txt"), "w").write("\n".join(out) + "\n")
