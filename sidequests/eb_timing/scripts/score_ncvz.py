"""
Score the Northern CVZ test: recall on the Borkovits+2026 L1 / D1 systems and the false-alarm
rate on the seeded random NCVZ EB sample, for the v1 (pre-registered) and v2 classification.

Inputs: results/noise_ncvz.csv (noise_fit on results/oc_ncvz, CAL=1), results/calibration.csv
(calibrate.py), data/ncvz_fa_list.csv, data/ncvz_fa_vetting.csv (by-eye verdicts).
'Testable' = timed with >= 10 sectors (13 L1 / 12 D1 stars have <= 2 sectors of MAST light
curves; B26 used their own FFI photometry for them).
TIC 377307592 (B26 L3, in the FA list by a harvesting defect) is excluded from the FA rate.
Output: results/ncvz_summary.json, results/ncvz_fa_results.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
EXCLUDE_FA = {377307592}


def main():
    n = pd.read_csv(HERE / "results" / "noise_ncvz.csv").set_index("tic")
    cal = pd.read_csv(HERE / "results" / "calibration.csv").set_index("tic")
    fa = pd.read_csv(HERE / "data" / "ncvz_fa_list.csv", comment="#")
    vet = pd.read_csv(HERE / "data" / "ncvz_fa_vetting.csv").set_index("tic")
    out = {}
    for grp in ("L1", "D1"):
        c = cal[cal.group == grp].copy()
        c["testable"] = c.n_sectors.fillna(0) >= 10
        for v in ("cls_v1", "cls"):
            c[v] = [n.loc[t, v] if t in n.index else "no timing" for t in c.index]
        t = c[c.testable]
        ok = t.P3_recovered.fillna(False).astype(bool)
        out[grp] = dict(n=int(len(c)), testable=int(len(t)), P3_recovered=int(ok.sum()),
                        flagged_and_recovered_v1=int((ok & t.cls_v1.isin(["LTTE", "LTTE-long"])).sum()),
                        flagged_and_recovered_v2=int((ok & t.cls.isin(["LTTE", "LTTE-long"])).sum()),
                        not_testable=[int(x) for x in c.index[~c.testable]])
    f = fa[~fa.TIC.isin(EXCLUDE_FA)]
    rows = []
    for tic in f.TIC:
        r = n.loc[tic] if tic in n.index else None
        rows.append(dict(tic=tic, n_sectors=r.n_sectors if r is not None else 0,
                         cls_v1=r.cls_v1 if r is not None else "no timing", cls_v2=r.cls if r is not None else "no timing",
                         v2_veto=r.v2_veto if r is not None else "", P3=r.P3 if r is not None else np.nan,
                         A_s=r.A_s if r is not None else np.nan, dBIC=r.dBIC if r is not None else np.nan,
                         verdict=vet.loc[tic, "verdict"] if tic in vet.index else ""))
    fr = pd.DataFrame(rows)
    fr.to_csv(HERE / "results" / "ncvz_fa_results.csv", index=False)
    fl1 = fr.cls_v1.isin(["LTTE", "LTTE-long"])
    fl2 = fr.cls_v2.isin(["LTTE", "LTTE-long"])
    pfa = fr.verdict == "pipeline-false-alarm"
    real = fr.verdict.str.startswith("real-signal")
    out["FA"] = dict(n=int(len(fr)), timed=int((fr.n_sectors > 0).sum()),
                     flagged_v1=int(fl1.sum()), false_alarms_v1=int((fl1 & pfa).sum()), real_v1=int((fl1 & real).sum()),
                     flagged_v2=int(fl2.sum()), false_alarms_v2=int((fl2 & pfa).sum()), real_v2=int((fl2 & real).sum()),
                     real_unflagged_any=[int(x) for x in fr.tic[real & ~fl2]],
                     note="v2 rules were derived from these same FA systems, so the v2 FA rate is in-sample")
    (HERE / "results" / "ncvz_summary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    print(fr[fl1 | real].to_string())


if __name__ == "__main__":
    main()
