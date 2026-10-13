"""
Out-of-sample false-alarm rate from NCVZ false-alarm set 2 (data/ncvz_fa2_list.csv), using the
frozen classification (ocfit v1 and v2) on results/noise_ncvz_fa2.csv, and the by-eye verdicts
in data/ncvz_fa2_vetting.csv.
Rates: flagged / N, and pipeline false alarms / N (verdict 'pipeline-false-alarm'), with exact
Clopper-Pearson 95% intervals.
Output: results/ncvz_fa2_summary.json, results/ncvz_fa2_results.csv
"""
import json
from pathlib import Path

import pandas as pd
from scipy.stats import beta

HERE = Path(__file__).resolve().parents[1]


def cp(k, n, a=0.05):
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k)
    return [round(k / n, 3), round(float(lo), 3), round(float(hi), 3)]


def main():
    fa = pd.read_csv(HERE / "data" / "ncvz_fa2_list.csv", comment="#")
    n = pd.read_csv(HERE / "results" / "noise_ncvz_fa2.csv").set_index("tic")
    vf = HERE / "data" / "ncvz_fa2_vetting.csv"
    vet = pd.read_csv(vf).set_index("tic") if vf.exists() else pd.DataFrame(columns=["verdict"])
    rows = []
    for tic in fa.TIC:
        r = n.loc[tic] if tic in n.index else None
        rows.append(dict(tic=tic, timed=r is not None, n_sectors=r.n_sectors if r is not None else 0,
                         cls_v1=r.cls_v1 if r is not None else "no timing", cls_v2=r.cls if r is not None else "no timing",
                         v2_veto=r.v2_veto if r is not None else "", P3=r.P3 if r is not None else None,
                         A_s=r.A_s if r is not None else None, A_sig=r.A_sig if r is not None else None,
                         dBIC=r.dBIC if r is not None else None, rho=r.rho if r is not None else None,
                         verdict=vet.loc[tic, "verdict"] if tic in vet.index else ""))
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "results" / "ncvz_fa2_results.csv", index=False)
    N = len(df)
    out = {"N": N, "timed": int(df.timed.sum())}
    for v in ("v1", "v2"):
        fl = df[f"cls_{v}"].isin(["LTTE", "LTTE-long"])
        pfa = fl & (df.verdict == "pipeline-false-alarm")
        real = fl & df.verdict.str.startswith("real")
        out[v] = {"flagged": int(fl.sum()), "flag_rate_CI": cp(int(fl.sum()), N),
                  "pipeline_false_alarms": int(pfa.sum()), "false_alarm_rate_CI": cp(int(pfa.sum()), N),
                  "real_or_plausible": int(real.sum()), "flagged_tics": [int(x) for x in df.tic[fl]]}
    (HERE / "results" / "ncvz_fa2_summary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    print(df[df.cls_v1.isin(["LTTE", "LTTE-long"]) | df.cls_v2.isin(["LTTE", "LTTE-long"])].to_string())


if __name__ == "__main__":
    main()
