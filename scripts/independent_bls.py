"""
Independent re-derivation of vetted periods (CLAUDE.md rule: a claimed period must be
re-derived by a method that has not seen the first analysis).

Deliberately shares no code with alias_pipeline.py / vet_candidate.py:
  - own light-curve loader: one product per sector (SPOC 2-min > TESS-SPOC > QLP), quality==0
  - own detrending: 1-day rolling median per contiguous data segment, then 30-min binning
  - blind BLS period search (astropy BoxLeastSquares) over 5-200 d, durations bracketing the
    catalog duration; NO transit times or alias grid are given to the search.
Reports the top BLS peaks and whether the best peak matches the claimed period (or a
simple harmonic of it).

Usage:  python scripts/independent_bls.py TIC P_claim [TIC P_claim ...]
Output: results/vetting/independent_bls.csv and results/vetting/<TIC>_bls.png
"""
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares

ROOT = Path(__file__).resolve().parent.parent
LC = ROOT / "data" / "raw" / "lightcurves"
CAT = ROOT / "data" / "catalog"
OUT = ROOT / "results" / "vetting"
RANK = {"SPOC_120": 0, "TESS-SPOC": 1, "QLP": 2}


def load(tic):
    best = {}
    for f in (LC / str(tic)).glob("s*.csv.gz"):
        m = re.fullmatch(r"s(\d+)_(.+)_(\d+)s\.csv\.gz", f.name)
        if not m or int(m.group(1)) > 200:
            continue
        key = "SPOC_120" if (m.group(2) == "SPOC" and m.group(3) == "120") else m.group(2)
        if key not in RANK:
            continue
        s = int(m.group(1))
        if s not in best or RANK[key] < RANK[best[s][0]]:
            best[s] = (key, f)
    t_all, f_all = [], []
    for s, (key, f) in sorted(best.items()):
        d = pd.read_csv(f, usecols=["time", "flux", "quality"])
        d = d[(d.quality == 0) & np.isfinite(d.flux) & np.isfinite(d.time)].sort_values("time")
        if len(d) < 100:
            continue
        seg = np.concatenate([[0], np.cumsum(np.diff(d.time.values) > 0.5)])
        for g in np.unique(seg):
            x = d[seg == g]
            if key == "QLP":           # raw SAP ramps at segment starts
                x = x[x.time - x.time.min() > 0.5]
            if x.time.max() - x.time.min() < 1.5:
                continue
            ser = pd.Series(x.flux.values, index=pd.to_datetime(x.time.values, unit="D"))
            trend = ser.rolling("1D", center=True, min_periods=10).median().values
            fl = x.flux.values / trend
            ok = np.isfinite(fl)
            # 30-min bins
            tb = np.floor(x.time.values[ok] * 48)
            df = pd.DataFrame({"tb": tb, "t": x.time.values[ok], "f": fl[ok]}).groupby("tb").mean()
            t_all.append(df.t.values)
            f_all.append(df.f.values)
    t, f = np.concatenate(t_all), np.concatenate(f_all)
    # clip upward outliers only
    med, mad = np.median(f), 1.4826 * np.median(np.abs(f - np.median(f)))
    ok = f < med + 5 * mad
    return t[ok], f[ok]


def run(tic, P_claim, cand):
    c = cand.loc[tic]
    t, f = load(tic)
    D = c.Dur_h / 24
    durations = np.array([0.6, 0.8, 1.0, 1.25]) * D
    baseline = t.max() - t.min()
    pmin, pmax = 5.0, 200.0
    # uniform frequency grid: a frequency error df shifts the transit phase by df*baseline
    # cycles = df*baseline*P days; require < D_min/4 at the shortest period searched
    df = durations.min() / (4 * baseline * pmin)
    freqs = np.arange(1 / pmax, 1 / pmin, df)
    periods = np.sort(1 / freqs)
    dy = np.full_like(f, 1.4826 * np.median(np.abs(f - np.median(f))))  # robust per-bin scatter
    bls = BoxLeastSquares(t, f, dy)
    res = bls.power(periods, durations, objective="snr")
    order = np.argsort(res.power)[::-1]
    # top distinct peaks (merge within 0.5%)
    peaks = []
    for i in order:
        p = res.period[i]
        if all(abs(p - q) / q > 0.005 for q, _ in peaks):
            peaks.append((p, res.power[i]))
        if len(peaks) == 5:
            break
    p_best = peaks[0][0]
    ratio = p_best / P_claim
    harm = min([1, 2, 3, 0.5, 1 / 3], key=lambda h: abs(ratio - h))
    match = abs(ratio - 1) < 0.001
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.5))
    ax[0].plot(res.period, res.power, lw=0.5)
    ax[0].axvline(P_claim, color="r", ls=":", lw=1)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("period (d)")
    ax[0].set_title(f"TIC {tic}: BLS (red: claimed {P_claim:.4f} d)", fontsize=9)
    i = np.argmax(res.power)
    ph = ((t - res.transit_time[i]) / p_best + 0.5) % 1 - 0.5
    ax[1].plot(ph * p_best * 24, (f - 1) * 1e3, ".", ms=1.5)
    ax[1].set_xlim(-3 * c.Dur_h, 3 * c.Dur_h)
    ax[1].set_xlabel("hours from mid-transit")
    ax[1].set_title(f"folded at BLS best {p_best:.5f} d", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / f"{tic}_bls.png", dpi=85)
    plt.close(fig)
    return dict(TIC=tic, P_claim=P_claim, P_bls=round(p_best, 5), bls_snr=round(peaks[0][1], 1),
                second_peak=round(peaks[1][0], 4), second_snr=round(peaks[1][1], 1),
                ratio=round(ratio, 5), nearest_harmonic=harm, exact_match=match,
                n_points=len(t), baseline_d=round(baseline, 1),
                top5=" ".join(f"{p:.4f}({s:.0f})" for p, s in peaks))


def main():
    cand = pd.read_csv(CAT / "candidates.csv").set_index("TIC")
    a = sys.argv[1:]
    rows = [run(int(a[i]), float(a[i + 1]), cand) for i in range(0, len(a), 2)]
    out = pd.DataFrame(rows)
    path = OUT / "independent_bls.csv"
    if path.exists():
        old = pd.read_csv(path)
        out = pd.concat([old[~old.TIC.isin(out.TIC)], out], ignore_index=True)
    out.to_csv(path, index=False)
    pd.set_option("display.width", 250)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()


def alias_scan(tic, cand, n_claim):
    """Second, less blind check: BLS SNR (independent detrending) evaluated only at the
    duotransit aliases P_n = dT/n (+/- a small window for ephemeris error)."""
    c = cand.loc[tic]
    t, f = load(tic)
    dy = np.full_like(f, 1.4826 * np.median(np.abs(f - np.median(f))))
    bls = BoxLeastSquares(t, f, dy)
    D = c.Dur_h / 24
    dT = c.Tc2 - c.Tc1
    rows = []
    for n in range(1, int(dT / 5) + 1):
        P = dT / n
        grid = P * (1 + np.linspace(-2e-4, 2e-4, 21))
        r = bls.power(grid, np.array([0.8, 1.0]) * D, objective="snr")
        rows.append((n, P, float(r.power.max())))
    s = pd.DataFrame(rows, columns=["n", "P", "snr"]).sort_values("snr", ascending=False)
    rank = int(np.where(s.n.values == n_claim)[0][0]) + 1
    return s, rank
