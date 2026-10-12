"""
Re-rank the 1,191 selected EBs (data/sample_all.csv, selected == True) by expected
light-travel-time sensitivity instead of raw eclipse count. Does NOT touch the pre-registered
pilot list or sample_top200.csv.

Per-eclipse timing error for a V-shaped eclipse of depth delta and total duration W, with
white noise sigma_1h (relative flux in 1 hour):  sigma_T = K * sigma_1h * sqrt(1h * W) / (2 delta)
(slope 2 delta / W on both branches; independent of cadence while cadence << W/2).
  sigma_1h(Tmag) = sqrt( (s10 * 10^(0.2 (Tmag - 10)))^2 + floor^2 ), s10 = 230 ppm, floor 60 ppm
  (rough TESS 1-hour CDPP; only the ranking matters). W = Wp-2g * P (catalogue width, phase
  units); delta = Dp-2g.
K is calibrated on the 27 systems we timed (median measured sigma_pri / predicted).
LTTE amplitude error from N_eff eclipses: sigma_A ~ sigma_T * sqrt(2 / N_eff),
  N_eff = N_ecl * (1 + min(1, (Ds/Dp)^2))  (secondaries count if comparably deep).
Ranking: ascending sigma_A (smallest detectable light-travel amplitude). Spots are not
modelled (Prsa has no Teff), so cool short-period systems will still be over-ranked.
Output: data/sample_ranked_precision.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
D = HERE / "data"


def sigma_1h(T):
    return np.sqrt((230e-6 * 10 ** (0.2 * (T - 10))) ** 2 + 60e-6 ** 2)


def predict(df):
    W = (df["Wp-2g"] * df.Per).clip(lower=0.01)
    delta = df["Dp-2g"].clip(lower=1e-3, upper=0.9)
    return sigma_1h(df.Tmag) * np.sqrt(W / 24.0) / (2 * delta)   # days (W in d, 1 h = 1/24 d)


def main():
    a = pd.read_csv(D / "sample_all.csv")
    a["sigT_pred_raw_s"] = predict(a) * 86400
    # calibrate on measured systems
    meas = []
    for f in sorted((HERE / "results" / "oc").glob("*_fit.json")):
        r = json.loads(f.read_text())
        meas.append((r["tic"], r["median_sig_pri_s"]))
    m = pd.DataFrame(meas, columns=["TIC", "sig_meas_s"]).merge(a[["TIC", "sigT_pred_raw_s"]], on="TIC")
    K = float(np.median(m.sig_meas_s / m.sigT_pred_raw_s))
    scatter = float(np.std(np.log10(m.sig_meas_s / (K * m.sigT_pred_raw_s))))
    a["sigT_pred_s"] = K * a.sigT_pred_raw_s
    ratio = (a["Ds-2g"] / a["Dp-2g"]).fillna(0)
    a["N_eff"] = a.N_ecl * (1 + np.minimum(1, ratio ** 2))
    a["sigA_pred_s"] = a.sigT_pred_s * np.sqrt(2 / a.N_eff)
    sel = a[a.selected].sort_values("sigA_pred_s").reset_index(drop=True)
    sel.insert(0, "rank_precision", np.arange(1, len(sel) + 1))
    old = pd.read_csv(D / "sample_top200.csv")[["TIC", "rank"]].rename(columns={"rank": "rank_eclipse_count"})
    sel = sel.merge(old, on="TIC", how="left")
    cols = ["rank_precision", "rank_eclipse_count", "TIC", "RAJ2000", "DEJ2000", "Tmag", "Per", "Morph", "Dp-2g",
            "Ds-2g", "Wp-2g", "n_sectors", "span_d", "N_ecl", "N_eff", "sigT_pred_s", "sigA_pred_s", "region",
            "known_multiple"]
    with open(D / "sample_ranked_precision.csv", "w") as fh:
        fh.write(f"# K={K:.3f} calibrated on {len(m)} timed systems; log10 scatter of measured/predicted = {scatter:.2f} dex\n")
        sel[cols].to_csv(fh, index=False, float_format="%.6g")
    print(f"K={K:.3f}, scatter {scatter:.2f} dex over {len(m)} systems")
    print(m.assign(pred=K * m.sigT_pred_raw_s).round(1).to_string())
    print(sel[cols[:3] + ["Tmag", "Per", "Dp-2g", "sigT_pred_s", "sigA_pred_s", "region"]].head(25).to_string())
    print("old top-20 pilot ranks under new ranking:",
          sorted(sel[sel.rank_eclipse_count <= 20].rank_precision.tolist()))
    print("overlap new top-200 with old top-200:", int(sel.head(200).rank_eclipse_count.notna().sum()))


if __name__ == "__main__":
    main()
