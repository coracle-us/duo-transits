"""
AAS82225 revision: compare the blind search (search_blind.py, search_periodic.py) with the submitted note.

For every star: was the catalogued event recovered blindly (any first-stage event; significant event), how many
significant events, the linked families with their allowed periods, and BLS detections. For the 15 stars in the
submitted Table 1 (table_submitted.csv): is each table period among the blind family's allowed periods
(within 0.05%), and how many periods does the blind search allow. Stars not in Table 1 with a blind family of
>= 2 significant events, or a BLS detection, are listed as new leads for vetting.
Output: k2singles/revision2/data/compare.csv, compare_table1.csv, new_leads.csv; summary on stdout.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"


def main():
    tab = pd.read_csv(HERE / "table_submitted.csv")
    tab_p = {int(r.EPIC): [float(x) for x in str(r.periods).split(";")] for r in tab.itertuples()}
    thr = json.loads((DATA / "thresholds.json").read_text())
    rows, t1, leads = [], [], []
    for f in sorted((DATA / "per_star").glob("*.json")):
        r = json.loads(f.read_text())
        if r.get("status") != "ok":
            rows.append(dict(EPIC=r["EPIC"], status=r.get("status")))
            continue
        p = json.loads((DATA / "periodic" / f.name).read_text())
        sig_T = [s["T"] for s in p["significant"]]
        cat = r["catalogued_events"]
        cat_hit = [e["blind_hit"] for e in cat]
        cat_sig = [any(abs(T - e["epoch"]) < max(0.1, (e["blind_D_h"] or 0) / 48) for T in sig_T) for e in cat]
        fams = [fm for fm in p["families"] if fm["n_members"] >= 2]
        bls = [s for b in p["bls"] for s in b["signals"] if s.get("detected")]
        row = dict(EPIC=r["EPIC"], TIC=r["TIC"], comment=r.get("comment"), k2_campaigns=" ".join(map(str, r["k2_campaigns"])),
                   n_tess_sectors=len(r["tess_sectors"]), n_events=p["n_events"], n_significant=p["n_significant"],
                   catalogued_recovered=all(cat_hit), catalogued_significant=all(cat_sig),
                   n_linked_families=len(fams), max_family=max((fm["n_members"] for fm in fams), default=1),
                   min_allowed=min((fm["n_allowed"] for fm in fams), default=None),
                   n_bls=len(bls), bls_periods=" ".join(f"{s['P']:.4f}(SDE {s['sde']:.1f})" for s in bls),
                   in_table1=r["EPIC"] in tab_p)
        rows.append(row)
        if r["EPIC"] in tab_p:
            allowed = [q["P"] for fm in fams for q in fm["periods"]]
            for P in tab_p[r["EPIC"]]:
                ok = any(abs(a / P - 1) < 5e-4 for a in allowed)
                bls_match = any(abs(s["P"] / P - 1) < 5e-3 or abs(s["P"] * 2 / P - 1) < 5e-3 for s in bls)
                t1.append(dict(EPIC=r["EPIC"], table_period=P, blind_allowed=ok, blind_bls=bls_match,
                               n_allowed_total=len(allowed),
                               families=" | ".join(f"{fm['n_members']}ev {fm['depth_ppm']:.0f}ppm {fm['D_h']}h: {fm['n_allowed']} allowed"
                                                   for fm in fams)))
        elif fams or bls:
            for fm in fams:
                leads.append(dict(EPIC=r["EPIC"], TIC=r["TIC"], kind="family", n_members=fm["n_members"],
                                  depth_ppm=fm["depth_ppm"], D_h=fm["D_h"], missions=" ".join(fm["missions"]),
                                  T=" ".join(map(str, fm["T"])), n_allowed=fm["n_allowed"],
                                  shortest=min((q["P"] for q in fm["periods"]), default=None),
                                  longest=max((q["P"] for q in fm["periods"]), default=None)))
            for s in bls:
                leads.append(dict(EPIC=r["EPIC"], TIC=r["TIC"], kind="BLS", depth_ppm=round(s["depth_ppm"]), D_h=round(s["D_h"], 1),
                                  T=round(s["T0"], 4), shortest=s["P"], longest=s["P"], n_allowed=1))
    c = pd.DataFrame(rows)
    c.to_csv(DATA / "compare.csv", index=False)
    pd.DataFrame(t1).to_csv(DATA / "compare_table1.csv", index=False)
    pd.DataFrame(leads).to_csv(DATA / "new_leads.csv", index=False)
    ok = c[c.get("status").isna()] if "status" in c else c
    print(f"thresholds: K2 S*={thr['K2']['S_star']}, TESS S*={thr['TESS']['S_star']}, BLS SDE={thr.get('BLS', {}).get('sde_threshold')}")
    print(f"stars: {len(c)}; searched ok: {len(ok)}")
    print(f"catalogued event recovered blindly: {int(ok.catalogued_recovered.sum())}; significant: {int(ok.catalogued_significant.sum())}")
    print(f"stars with a linked family (>= 2 significant events): {int((ok.n_linked_families > 0).sum())}; "
          f">= 3 events: {int((ok.max_family >= 3).sum())}; BLS detections: {int((ok.n_bls > 0).sum())}")
    if t1:
        t = pd.DataFrame(t1)
        print(f"Table 1 periods: {len(t)}; among blind allowed periods: {int(t.blind_allowed.sum())}; matched by BLS: {int(t.blind_bls.sum())}")
    print(f"new leads (not in Table 1): {len(set(l['EPIC'] for l in leads))} stars, {len(leads)} rows")


if __name__ == "__main__":
    main()
