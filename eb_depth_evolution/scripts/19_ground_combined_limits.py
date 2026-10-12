"""Combine WASP and KELT scan results (inverse-variance) per time offset; summarise depth limits for the primary and secondary
in 2004-2008 over plausible offset windows. Expected offset of the 2006-08 primary from the TESS linear ephemeris: linear
extrapolation of the K2 O-C (-102 min at E=-124) to E~-250 gives about -0.14 d."""
import numpy as np, pandas as pd
R = pd.read_csv("data/ground_eclipse_scan.csv")
out = []
for k in ["P", "S"]:
    a = R[(R.survey == "WASP") & (R.kind == k)].set_index("delta_d"); b = R[(R.survey == "KELT") & (R.kind == k)].set_index("delta_d")
    j = a[["depth", "err"]].join(b[["depth", "err"]], lsuffix="_w", rsuffix="_k", how="outer")
    w1 = 1/j.err_w**2; w2 = 1/j.err_k**2
    wsum = w1.fillna(0) + w2.fillna(0)
    j["depth"] = (j.depth_w.fillna(0)*w1.fillna(0) + j.depth_k.fillna(0)*w2.fillna(0))/wsum
    j["err"] = 1/np.sqrt(wsum)
    j.loc[wsum == 0, ["depth", "err"]] = np.nan
    j["kind"] = k; out.append(j.reset_index())
    for lo, hi in [(-0.3, 0.1), (-0.6, 0.3), (-0.95, 0.35), (-1.5, 1.5)]:
        h = j[(j.index >= lo) & (j.index <= hi)]
        cov = h.err.notna().mean()
        lim = np.nanmax(h.depth + 3*h.err)
        print(f"{k}: offsets [{lo:+.2f},{hi:+.2f}] d: combined coverage {cov*100:.0f}%, max fitted depth {np.nanmax(h.depth):.4f}, "
              f"worst 3-sigma upper limit {lim:.4f}, median 3-sigma limit {np.nanmedian(h.depth+3*h.err):.4f}")
pd.concat(out).to_csv("data/ground_combined_scan.csv", index=False)
