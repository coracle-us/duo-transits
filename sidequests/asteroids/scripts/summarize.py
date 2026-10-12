"""
Summary of validation + pilot: success vs V, and the noise floor.
Inputs: data/validation_results_classified.csv, data/pilot_results_classified.csv, data/pilot_manual.csv
Outputs: plots/success_vs_V.png, data/success_vs_V.csv, data/pilot_vs_lcdb.csv
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
D = HERE / "data"


def main():
    v = pd.read_csv(D / "validation_results_classified.csv")
    v["sample"] = "validation (LCDB U=3)"
    v["correct"] = v.auto_pass & (v.match == "P")
    p = pd.read_csv(D / "pilot_results_classified.csv")
    p["sample"] = "pilot (Trojans, S72)"
    s = pd.read_csv(D / "pilot_sample.csv")[["pdes", "P_lcdb_h", "U_lcdb"]].rename(columns={"pdes": "id"})
    s["id"] = s.id.astype(int)
    man = pd.read_csv(D / "pilot_manual.csv").rename(columns={"number": "id"})
    p = p.merge(s, on="id", how="left").merge(man[["id", "eye_verdict", "P_rot_h_adopted", "U_eye"]], on="id")
    r = p.P2_h / p.P_lcdb_h
    p["vs_lcdb"] = np.where(p.P_lcdb_h.isna(), "no LCDB entry",
                            np.where(np.abs(r - 1) < 0.01, "matches LCDB",
                                     np.where(np.abs(r - 2) < 0.02, "2x LCDB",
                                              np.where(np.abs(r - 0.5) < 0.005, "0.5x LCDB", "differs"))))
    p[["id", "label", "set", "V_med", "P2_h", "sP2_h", "amp2", "boot_fap", "frac_contam", "auto_pass",
       "eye_verdict", "U_eye", "P_lcdb_h", "U_lcdb", "vs_lcdb"]].to_csv(D / "pilot_vs_lcdb.csv", index=False)
    allr = pd.concat([v, p], ignore_index=True)
    bins = [10, 16, 17, 18, 19.5]
    allr["Vbin"] = pd.cut(allr.V_med, bins, right=False)
    rows = []
    for (smp, vb), g in allr.groupby(["sample", "Vbin"], observed=True):
        rows.append(dict(sample=smp, V_bin=str(vb), N=len(g), N_auto_pass=int(g.auto_pass.sum()),
                         median_resid_rms_mag=round(float(g.resid_rms_mag.median()), 3),
                         median_frac_contam=round(float(g.frac_contam.median()), 2)))
    tab = pd.DataFrame(rows)
    tab.to_csv(D / "success_vs_V.csv", index=False)
    print(tab.to_string(index=False))
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for smp, mk in [("validation (LCDB U=3)", "o"), ("pilot (Trojans, S72)", "s")]:
        g = allr[allr["sample"] == smp]
        for ok, col in [(True, "tab:blue"), (False, "tab:red")]:
            h = g[g.auto_pass == ok]
            axs[0].scatter(h.V_med, h.resid_rms_mag, marker=mk, c=col, s=28, alpha=0.8,
                           label=f"{smp}, {'pass' if ok else 'fail'}")
            axs[1].scatter(h.V_med, h.frac_contam, marker=mk, c=col, s=28, alpha=0.8)
    axs[0].set_yscale("log")
    axs[0].set_xlabel("predicted V (Horizons)")
    axs[0].set_ylabel("residual rms, 30-min bins (mag)")
    axs[0].legend(fontsize=7)
    axs[1].axhline(0.4, color="k", ls=":", lw=0.8)
    axs[1].set_xlabel("predicted V (Horizons)")
    axs[1].set_ylabel("fraction of good frames star-contaminated")
    fig.suptitle("TESS FFI moving-target photometry: automatic acceptance vs brightness", fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / "plots" / "success_vs_V.png", dpi=90)


if __name__ == "__main__":
    main()
