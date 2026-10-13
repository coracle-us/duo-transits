"""
1. Inventory of TESS light-curve products (all authors) for every target TIC (lightkurve search on MAST).
2. Empirical 1-hour noise calibration: download one light curve per target where a TESS-SPOC or QLP FFI
   product exists, flatten (1-day window), bin to 1 h, and take 1.4826*MAD of the binned flux.
   Fit log10(sigma_1h) vs Tmag and compare with the project's old model 220 ppm * 10^(0.2 (T - 9.1)).

Outputs: data/tess_lc_inventory.csv, data/noise_calibration.csv, data/noise_model.txt
"""
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parents[1]
CACHE = HERE / "data" / "lk_cache"


def inventory(tics):
    import lightkurve as lk
    out = []
    for tic in tics:
        for attempt in range(3):
            try:
                r = lk.search_lightcurve(f"TIC {tic}", mission="TESS")
                break
            except Exception:  # noqa: BLE001
                time.sleep(5)
                r = None
        if r is None or len(r) == 0:
            out.append(dict(TIC=tic, author="", exptime=np.nan, sector=np.nan))
            continue
        t = r.table
        for i in range(len(t)):
            m = str(t["mission"][i])
            out.append(dict(TIC=tic, author=str(t["author"][i]), exptime=float(t["exptime"][i]),
                            sector=int(m.split()[-1]) if "Sector" in m else np.nan))
    return pd.DataFrame(out)


def sigma_1h(lc):
    lc = lc.remove_nans()
    if len(lc) < 200:
        return np.nan
    lc = lc.normalize()
    dt = np.nanmedian(np.diff(lc.time.value))
    win = int(round(1.0 / dt)) | 1
    flat = lc.flatten(window_length=max(win, 5))
    b = flat.bin(time_bin_size=1 / 24)
    f = b.flux.value
    f = f[np.isfinite(f)]
    if len(f) < 50:
        return np.nan
    return 1.4826 * np.nanmedian(np.abs(f - np.nanmedian(f))) * 1e6


def main():
    tg = pd.read_csv(HERE / "data" / "kepler_singles_targets.csv")
    stars = tg.drop_duplicates("TIC")[["TIC", "KIC", "Tmag"]]
    invf = HERE / "data" / "tess_lc_inventory.csv"
    if invf.exists():
        inv = pd.read_csv(invf)
    else:
        inv = inventory(stars.TIC.tolist())
        inv.to_csv(invf, index=False)
    print(inv.groupby("author").agg(n=("TIC", "size"), stars=("TIC", "nunique")))

    import lightkurve as lk
    rows = []
    CACHE.mkdir(exist_ok=True)
    ncf = HERE / "data" / "noise_calibration.csv"
    for s in ([] if ncf.exists() else stars.itertuples()):
        for author in ["TESS-SPOC", "QLP", "SPOC"]:
            sub = inv[(inv.TIC == s.TIC) & (inv.author == author)]
            if not len(sub):
                continue
            # one recent sector (200-s era if possible) and one early sector
            secs = sorted(sub.sector.dropna().astype(int).unique())
            pick = sorted({secs[-1], secs[0]})
            for sec in pick:
                try:
                    r = lk.search_lightcurve(f"TIC {s.TIC}", mission="TESS", author=author, sector=sec)
                    lc = r[0].download(download_dir=str(CACHE))
                    col = "pdcsap_flux" if "pdcsap_flux" in lc.colnames else (
                        "kspsap_flux" if "kspsap_flux" in lc.colnames else "flux")
                    lc = lc.select_flux(col)
                    q = lc.quality.value == 0 if hasattr(lc, "quality") else slice(None)
                    sig = sigma_1h(lc[q])
                except Exception as e:  # noqa: BLE001
                    print("fail", s.TIC, author, sec, e)
                    sig = np.nan
                rows.append(dict(TIC=s.TIC, KIC=s.KIC, Tmag=s.Tmag, author=author, sector=sec, sigma_1h_ppm=sig))
            break  # one author per star
    if rows:
        pd.DataFrame(rows).to_csv(ncf, index=False)
    nc = pd.read_csv(ncf)
    nc["sigma_1h_ppm"] = pd.to_numeric(nc.sigma_1h_ppm, errors="coerce")
    g = nc[np.isfinite(nc.sigma_1h_ppm) & (nc.sigma_1h_ppm > 0)]
    # robust fit log10 sigma = a + b (T - 13)
    x, y = g.Tmag.values - 13, np.log10(g.sigma_1h_ppm.values)
    A = np.vstack([np.ones_like(x), x]).T
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    for _ in range(3):
        res = y - A @ coef
        keep = np.abs(res) < 3 * 1.4826 * np.median(np.abs(res))
        coef = np.linalg.lstsq(A[keep], y[keep], rcond=None)[0]
    old = lambda T: 220 * 10 ** (0.2 * (T - 9.1))  # noqa: E731
    lines = [f"empirical fit (N={keep.sum()} of {len(g)}): log10 sigma_1h[ppm] = {coef[0]:.3f} + {coef[1]:.3f} (Tmag - 13)",
             "Tmag  empirical  old_model  (ppm, 1 h)"]
    for T in [10, 11, 12, 13, 14, 15, 16]:
        lines.append(f"{T:4d}  {10 ** (coef[0] + coef[1] * (T - 13)):9.0f}  {old(T):9.0f}")
    by = g.assign(Tbin=np.floor(g.Tmag)).groupby(["Tbin"]).sigma_1h_ppm.agg(["size", "median"])
    lines.append("binned medians:\n" + by.round(0).to_string())
    by2 = g.assign(era=np.where(g.sector >= 56, "S56+ (200 s)", "S14-55")).groupby("era").apply(
        lambda d: np.median(np.log10(d.sigma_1h_ppm) - (coef[0] + coef[1] * (d.Tmag - 13))))
    lines.append("median log residual by era:\n" + by2.round(3).to_string())
    (HERE / "data" / "noise_model.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
