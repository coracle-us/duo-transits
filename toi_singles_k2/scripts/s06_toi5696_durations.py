"""
TOI-5696 (TIC 393546540 = EPIC 248618202) check: which TESS single does the K2 C14 event pair with?

ExoFOP pairs the K2 C14 event (BTJD ~925.87) with the S45/46 event as TOI-5696.01 (P_max = 1637.27 d,
T14 = 5.98 h) and lists the S72 event as a separate single, TOI-5696.03 (T14 = 8.09 h; the PHT CTOI note
says its duration differs from the S46 single). Our K2 box search with the .03 duration finds the
K2 C14 event, so here we fit a trapezoid (eclipse_times.fit_event) to each event in each reduction and
compare T14 and depth.
K2 at 30-min cadence: the model is supersampled (x15) to avoid long-cadence smearing biasing T14.

Output: toi_singles_k2/derived/toi5696_durations.csv, toi_singles_k2/plots/toi5696_events.png
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "toi_singles_k2"
sys.path.insert(0, str(ROOT / "k2singles" / "scripts"))
sys.path.insert(0, str(OUT / "scripts"))
import eclipse_times as et  # noqa: E402
import s03_k2_search as s3  # noqa: E402

TIC, EPIC = 393546540, 248618202


def fit_super(t, f, guess, D, cad, ns=15):
    m = np.abs(t - guess) < 1.5 * D + 0.3
    tt, ff = t[m], f[m]
    oot = np.abs(tt - guess) > 0.6 * D
    ff = ff / np.median(ff[oot])
    off = (np.arange(ns) - (ns - 1) / 2) / ns * cad

    def model(p):
        t0, depth, t14, fin, a, b = p
        ts = (tt[:, None] + off[None, :]).ravel()
        mm = et.trapezoid(ts, t0, depth, t14, fin).reshape(len(tt), ns).mean(axis=1)
        return mm * (a + b * (tt - guess))

    p0 = [guess, max(0.0005, 1 - np.min(ff)), D, 0.2, 1.0, 0.0]
    lo = [guess - 0.5 * D, 0, 0.3 * D, 0.01, 0.5, -1]
    hi = [guess + 0.5 * D, 1, 2.0 * D, 0.5, 1.5, 1]
    r = least_squares(lambda p: ff - model(p), p0, bounds=(lo, hi), x_scale="jac")
    res = ff - model(r.x)
    chi2r = float(np.sum(res ** 2) / max(1, len(tt) - 6))
    cov = np.linalg.inv(r.jac.T @ r.jac) * chi2r
    beta = et.beta_of(res, max(1, int(round(0.25 * r.x[2] / cad))))
    e = np.sqrt(np.diag(cov)) * beta
    return dict(T=r.x[0], sigma_T=e[0], depth_ppm=r.x[1] * 1e6, sigma_depth_ppm=e[1] * 1e6, T14_h=r.x[2] * 24,
                sigma_T14_h=e[2] * 24, fin=r.x[3], beta=beta, n=len(tt), _tt=tt, _ff=ff, _m=model(r.x))


def tess_lcs():
    import lightkurve as lk
    ddir = OUT / "data" / "tess" / str(TIC)
    ddir.mkdir(parents=True, exist_ok=True)
    out = {}
    for author in ["SPOC", "QLP"]:
        sr = lk.search_lightcurve(f"TIC {TIC}", mission="TESS", author=author)
        for i in range(len(sr)):
            row = sr[i]
            sec = int(row.table["sequence_number"][0]) if "sequence_number" in row.table.colnames else None
            if sec not in (45, 46, 72):
                continue
            exp = float(row.table["exptime"][0])
            if author == "SPOC" and exp > 120:
                continue
            try:
                lc = row.download(download_dir=str(ddir))
            except Exception as e:  # noqa: BLE001
                print("download fail", author, sec, e)
                continue
            col = "pdcsap_flux" if author == "SPOC" else ("kspsap_flux" if "kspsap_flux" in lc.columns else "flux")
            lc = lc.remove_nans() if col == "flux" else lc.select_flux(col).remove_nans()
            out[f"{author}{int(exp)}s_S{sec}"] = (np.asarray(lc.time.value, float), np.asarray(lc.flux.value, float) /
                                                 np.nanmedian(lc.flux.value), exp / 86400)
    return out


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows, panels = [], []
    pdc = s3.k2lib.k2_lightcurves(EPIC)
    k2 = {"K2PDC": s3._segments(pdc), "K2SAP": s3.load_sap(EPIC), "EVEREST": s3.load_hlsp(EPIC, "EVEREST"),
          "K2SFF": s3.load_hlsp(EPIC, "K2SFF")}
    for name, lc in k2.items():
        r = fit_super(lc.time.to_numpy(), lc.flux.to_numpy(), 925.88, 7 / 24, 0.02043)
        rows.append(dict(event="K2 C14", product=name, **{k: v for k, v in r.items() if not k.startswith("_")}))
        panels.append((f"K2 C14 {name}", r))
    tl = tess_lcs()
    for name, (t, f, cad) in tl.items():
        sec = int(name.split("_S")[1])
        guess = 3284.714 if sec == 72 else 2563.13
        if not (np.abs(t - guess) < 0.3).any():
            # find the TOI-5696.01 epoch in S45/46 by scanning
            continue
        r = fit_super(t, f, guess, 7 / 24, cad, ns=max(1, int(round(cad * 1440 / 2))))
        lab = "TESS S72 (.03)" if sec == 72 else f"TESS S{sec} (.01)"
        rows.append(dict(event=lab, product=name, **{k: v for k, v in r.items() if not k.startswith("_")}))
        panels.append((f"{lab} {name}", r))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "derived" / "toi5696_durations.csv", index=False)
    pd.set_option("display.width", 220)
    print(d.round(4).to_string())
    fig, axs = plt.subplots(len(panels), 1, figsize=(7, 2.0 * len(panels)), sharex=True)
    for ax, (lab, r) in zip(np.atleast_1d(axs), panels):
        ax.plot((r["_tt"] - r["T"]) * 24, r["_ff"], ".", ms=2, color="k")
        ax.plot((r["_tt"] - r["T"]) * 24, r["_m"], "-", color="tab:red")
        ax.set_title(f"{lab}: T14 {r['T14_h']:.2f} h, depth {r['depth_ppm']:.0f} ppm", fontsize=8)
    np.atleast_1d(axs)[-1].set_xlabel("hours from mid-transit")
    fig.tight_layout()
    fig.savefig(OUT / "plots" / "toi5696_events.png", dpi=110)


if __name__ == "__main__":
    main()
