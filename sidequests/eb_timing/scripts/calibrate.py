"""
Calibrate the O-C uncertainty model against published third-body solutions, and score the
Northern CVZ recall test.

Reference solutions:
  - NCVZ recall set: Borkovits+2026 L1 (pure LTTE; P2, A from a sin i, e2 with errors) and D1
    (P2, e2 only), from data/b26_triples.csv; our fits in results/noise_ncvz.csv (new model)
    and results/oc_ncvz/<TIC>_fit.json (old formal errors).
  - southern/validation systems with published LTTE (validate.PUB) are reported but not used for
    the calibration fit (held out).
Pull = (ours - published) / sqrt(sigma_ours^2 + sigma_pub^2), computed for P3 and A of every
system whose P3 is within 10% of the published P2 ('P3 recovered'); systems outside that
window are recall failures, not pulls.
CAL: the factor multiplying our raw bootstrap sigmas that makes the robust standard deviation
(1.4826 MAD) of the pooled L1 P3+A pulls equal to 1.
Output: results/calibration.csv (per-system pulls), results/calibration_summary.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import validate as va  # noqa: E402


def pulls(x, sx, y, sy):
    return (x - y) / np.sqrt(sx ** 2 + np.nan_to_num(sy) ** 2)


def stats(z):
    z = np.asarray(z, float)
    z = z[np.isfinite(z)]
    if not len(z):
        return {}
    return dict(n=int(len(z)), rms=float(np.sqrt(np.mean(z ** 2))), robust_sd=float(1.4826 * np.median(np.abs(z))),
                frac_lt1=float(np.mean(np.abs(z) < 1)), frac_lt2=float(np.mean(np.abs(z) < 2)),
                median=float(np.median(z)))


def main():
    b = pd.read_csv(HERE / "data" / "b26_triples.csv")
    b = b[b.solution == 1].set_index("tic")
    nf = pd.read_csv(HERE / "results" / "noise_ncvz.csv").set_index("tic")
    rec = pd.read_csv(HERE / "data" / "ncvz_recall_list.csv")
    rows = []
    for r in rec.itertuples():
        tic, grp = int(r.TIC), r.group
        f_old = HERE / "results" / "oc_ncvz" / f"{tic}_fit.json"
        row = dict(tic=tic, group=grp, P2_pub=b.loc[tic, "P2"], P2_pub_err=b.loc[tic, "P2_err"],
                   A_pub=b.loc[tic, "A_ltte_s"], A_pub_err=b.loc[tic, "A_ltte_err_s"])
        if tic not in nf.index or not f_old.exists():
            row["status"] = "no timing (insufficient light curves)"
            rows.append(row)
            continue
        n = nf.loc[tic]
        old = json.loads(f_old.read_text())["ltte"]
        row.update(n_sectors=n.n_sectors, cls=n.cls, P3=n.P3, P3_sig=n.P3_sig, A=n.A_s, A_sig=n.A_sig, e3=n.e3,
                   P3_old_err=old["P3_err"], A_old_err=old["asini_err_s"], dBIC=n.dBIC, snr=n.snr, rho=n.rho)
        ok = abs(n.P3 / row["P2_pub"] - 1) < 0.10
        row["P3_recovered"] = bool(ok)
        row["flagged"] = n.cls in ("LTTE", "LTTE-long")
        row["status"] = ("recovered" if ok and row["flagged"] else "P3 ok, not flagged" if ok
                         else "flagged, wrong P3" if row["flagged"] else "missed")
        if ok:
            row["pullP_old"] = pulls(n.P3, row["P3_old_err"], row["P2_pub"], row["P2_pub_err"])
            row["pullP_new_raw"] = pulls(n.P3, n.P3_sig, row["P2_pub"], row["P2_pub_err"])
            if grp == "L1":
                row["pullA_old"] = pulls(n.A_s, row["A_old_err"], row["A_pub"], row["A_pub_err"])
                row["pullA_new_raw"] = pulls(n.A_s, n.A_sig, row["A_pub"], row["A_pub_err"])
        rows.append(row)
    df = pd.DataFrame(rows)
    L1 = df[(df.group == "L1") & df.get("P3_recovered", pd.Series(False, index=df.index)).fillna(False).astype(bool)]

    # calibration factor on L1 pooled pulls (raw sigmas assumed uncalibrated, CAL=1 in noise_fit)
    def robust_sd(k):
        zP = (L1.P3 - L1.P2_pub) / np.sqrt((k * L1.P3_sig) ** 2 + L1.P2_pub_err.fillna(0) ** 2)
        zA = (L1.A - L1.A_pub) / np.sqrt((k * L1.A_sig) ** 2 + L1.A_pub_err.fillna(0) ** 2)
        z = np.concatenate([zP, zA])
        return 1.4826 * np.median(np.abs(z[np.isfinite(z)]))

    lo, hi = 0.1, 50.0
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        if robust_sd(mid) > 1:
            lo = mid
        else:
            hi = mid
    cal = float(np.sqrt(lo * hi))
    df["pullP_new_cal"] = np.where(df.get("P3_recovered", False) == True,  # noqa: E712
                                   pulls(df.P3, cal * df.P3_sig, df.P2_pub, df.P2_pub_err), np.nan)
    df["pullA_new_cal"] = np.where((df.get("P3_recovered", False) == True) & (df.group == "L1"),  # noqa: E712
                                   pulls(df.A, cal * df.A_sig, df.A_pub, df.A_pub_err), np.nan)
    df.to_csv(HERE / "results" / "calibration.csv", index=False)
    l1 = df[df.group == "L1"]
    summ = dict(
        CAL=cal,
        recall_L1=dict(n=int(len(l1)), timed=int(l1.status.ne("no timing (insufficient light curves)").sum()),
                       counts=l1.status.value_counts().to_dict()),
        recall_D1=dict(n=int((df.group == "D1").sum()), counts=df[df.group == "D1"].status.value_counts().to_dict()),
        pulls_L1_old=stats(np.concatenate([l1.pullP_old.dropna(), l1.pullA_old.dropna()])),
        pulls_L1_new_raw=stats(np.concatenate([l1.pullP_new_raw.dropna(), l1.pullA_new_raw.dropna()])),
        pulls_L1_new_cal=stats(np.concatenate([l1.pullP_new_cal.dropna(), l1.pullA_new_cal.dropna()])),
        pulls_D1_P_old=stats(df[df.group == "D1"].pullP_old.dropna()) if "pullP_old" in df else {},
        pulls_D1_P_new_cal=stats(df[df.group == "D1"].pullP_new_cal.dropna()),
    )
    # held-out check: the 7 validation systems, timed on the FULL baseline (to S98) and fitted
    # with the calibrated model in the southern run (results/noise_fit.csv), against B26 where
    # B26 has them (P2 for all groups; A only for L groups). 229771234 = B26's 229771231.
    vf = HERE / "results" / "noise_fit.csv"
    if vf.exists():
        nv = pd.read_csv(vf).set_index("tic")
        alias = {229771234: 229771231}
        vrows = []
        for tic in [424461577, 219900027, 229771234, 259006185, 356896561]:
            bt = alias.get(tic, tic)
            if tic not in nv.index or bt not in b.index:
                continue
            q, pub = nv.loc[tic], b.loc[bt]
            vr = dict(tic=tic, b26_group=pub.group, P3=q.P3, P3_sig=q.P3_sig, P2_pub=pub.P2, P2_pub_err=pub.P2_err,
                      pullP=pulls(q.P3, q.P3_sig, pub.P2, pub.P2_err))
            if str(pub.group).startswith("L"):
                vr.update(A=q.A_s, A_sig=q.A_sig, A_pub=pub.A_ltte_s, A_pub_err=pub.A_ltte_err_s,
                          pullA=pulls(q.A_s, q.A_sig, pub.A_ltte_s, pub.A_ltte_err_s))
            vrows.append(vr)
        vdf = pd.DataFrame(vrows)
        vdf.to_csv(HERE / "results" / "calibration_validation.csv", index=False)
        summ["validation_heldout"] = dict(pulls=stats(np.concatenate([vdf.pullP.dropna(), vdf.get("pullA", pd.Series(dtype=float)).dropna()])),
                                          note="validation systems use the S1/S14-98 baseline, B26 used S14-86")
        print(vdf.round(2).to_string())
    (HERE / "results" / "calibration_summary.json").write_text(json.dumps(summ, indent=1, default=float))
    print(json.dumps(summ, indent=1, default=float))
    with pd.option_context("display.width", 250):
        print(df[["tic", "group", "status", "n_sectors", "P2_pub", "P3", "P3_sig", "A_pub", "A", "A_sig", "pullP_old",
                  "pullP_new_cal", "pullA_old", "pullA_new_cal", "dBIC", "snr", "rho"]].round(2).to_string())


if __name__ == "__main__":
    main()
