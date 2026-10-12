"""
Point-by-point comparison of our eclipse times with an external list of times of minima
(e.g. the Borkovits+2026 NCVZ minima, "available on request"), plus an O-C fit of the external
times with the SAME noise model (ocfit.analyse) so P3/A can be compared like for like.

External file format (CSV): TIC, T, type[, err]
  T     BJD_TDB (values > 2.4e6 are converted to BTJD = BJD - 2457000) or BTJD
  type  'p'/'s' or 'I'/'II' or 1/2 (primary/secondary)
  err   1-sigma uncertainty in days (optional; if missing, the robust scatter of the external
        O-C about a running median is used)
Usage:
  python sidequests/eb_timing/scripts/compare_minima.py --ext=FILE [--oc=results/oc_ncvz] [TIC ...]
  python sidequests/eb_timing/scripts/compare_minima.py --selftest=TIC [--oc=...]
     (self-test: our own times + 20 s Gaussian noise + a 30 s offset play the external list)
Output: results/minima_compare.csv (one row per TIC) and figs/<TIC>_minima_cmp.png
Matching: our linear ephemeris (P, T0, secondary phase from <TIC>_fit.json) assigns each
external minimum an epoch; it is matched to our eclipse of the same type and epoch.
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ocfit  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
SECTORS = pd.read_csv(HERE / "data" / "sector_table.csv")


def to_btjd(t):
    t = np.asarray(t, float)
    return np.where(t > 2.4e6, t - 2457000.0, t)


def norm_type(x):
    s = str(x).strip().lower()
    return "s" if s in ("s", "ii", "2", "sec", "secondary") else "p"


def sector_of(t):
    idx = np.searchsorted(SECTORS.t_start.values, t, side="right") - 1
    idx = np.clip(idx, 0, len(SECTORS) - 1)
    inside = t <= SECTORS.t_end.values[idx] + 1.0
    return np.where(inside, SECTORS.sector.values[idx], -1)


def compare(tic, ext, ocdir):
    fit = json.loads((ocdir / f"{tic}_fit.json").read_text())
    P, T0, ph2 = fit["P"], fit["T0"], fit.get("ph_sec") or 0.5
    ours = pd.read_csv(ocdir / f"{tic}_oc.csv")
    ours = ours[ours.flag.isna() | (ours.flag == "")]
    e = ext[ext.TIC == tic].copy()
    e["T"] = to_btjd(e["T"])
    e["kind"] = e["type"].map(norm_type)
    ph = np.where(e.kind == "s", ph2, 0.0)
    e["epoch"] = np.round((e["T"] - T0 - ph * P) / P).astype(int)
    e["oc"] = e["T"] - (T0 + (e.epoch + ph) * P)
    m = e.merge(ours[["kind", "epoch", "T", "oc", "sig_s"]], on=["kind", "epoch"], suffixes=("_ext", "_ours"))
    diff = (m.T_ext - m.T_ours) * 86400
    row = dict(tic=tic, n_ext=len(e), n_matched=len(m), mean_diff_s=diff.mean(), median_diff_s=diff.median(),
               rms_diff_s=diff.std(), robust_diff_s=ocfit.robust_sigma(diff),
               median_sig_ours_s=(m.sig_s * 86400).median())
    if "err" in e and e["err"].notna().any():
        row["median_sig_ext_s"] = e["err"].median() * 86400
        z = diff / np.hypot(m.sig_s * 86400, m["err"] * 86400) if "err" in m else np.nan
        row["rms_pull"] = float(np.sqrt(np.mean(z ** 2))) if np.ndim(z) else np.nan
    # same-noise-model O-C fit of the external times
    if "err" not in e or e["err"].isna().all():
        r = e.oc - e.oc.rolling(15, center=True, min_periods=3).median()
        e["err"] = max(ocfit.robust_sigma(r), 1e-5)
    e["sector"] = sector_of(e["T"].values)
    e["sig"] = e["err"]
    if len(e) >= 30:
        fe = ocfit.analyse(e[["T", "kind", "sector", "sig", "oc"]], boot=60)
        row.update({f"ext_{k}": fe[k] for k in ("P3", "P3_sig", "A_s", "A_sig", "e3", "cls")})
    fo = ocfit.analyse(ours.drop(columns=["sig"]).rename(columns={"sig_s": "sig"})[["T", "kind", "sector", "sig", "oc"]], boot=60)
    row.update({f"ours_{k}": fo[k] for k in ("P3", "P3_sig", "A_s", "A_sig", "e3", "cls")})
    fig, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
    for k, c in (("p", "C0"), ("s", "C3")):
        mk = m.kind == k
        ax[0].plot(m.T_ours[mk], m.oc_ours[mk] * 1440, ".", color=c, ms=2, label=f"ours {k}")
        ax[0].plot(m.T_ours[mk], m.oc_ext[mk] * 1440, "x", color="k", ms=2)
        ax[1].plot(m.T_ours[mk], diff[mk], ".", color=c, ms=2)
    ax[0].set_ylabel("O-C (min)")
    ax[1].set_ylabel("ext - ours (s)")
    ax[1].set_xlabel("BTJD")
    ax[0].set_title(f"TIC {tic}: external (x) vs ours (.)", fontsize=9)
    fig.savefig(HERE / "figs" / f"{tic}_minima_cmp.png", dpi=90, bbox_inches="tight")
    plt.close(fig)
    return row


def main():
    args = sys.argv[1:]
    ocdir, ext, tics = HERE / "results" / "oc_ncvz", None, []
    for a in args:
        if a.startswith("--oc="):
            ocdir = HERE / a.split("=", 1)[1]
        elif a.startswith("--ext="):
            ext = pd.read_csv(a.split("=", 1)[1], comment="#")
        elif a.startswith("--selftest="):
            tic = int(a.split("=", 1)[1])
            tics = [tic]
            rng = np.random.default_rng(1)
            o = pd.read_csv(ocdir / f"{tic}_oc.csv")
            o = o[o.flag.isna() | (o.flag == "")]
            ext = pd.DataFrame({"TIC": tic, "T": o["T"] + 2457000.0 + 30 / 86400 + rng.normal(0, 20 / 86400, len(o)),
                                "type": o.kind.map({"p": "I", "s": "II"}), "err": 20 / 86400})
        else:
            tics.append(int(a))
    if ext is None:
        raise SystemExit("need --ext=FILE or --selftest=TIC")
    if not tics:
        tics = sorted(set(ext.TIC.astype(int)))
    rows = []
    for tic in tics:
        if not (ocdir / f"{tic}_oc.csv").exists():
            print(f"TIC {tic}: no O-C of ours in {ocdir}")
            continue
        r = compare(tic, ext, ocdir)
        print(r, flush=True)
        rows.append(r)
    name = "minima_compare_selftest.csv" if any(a.startswith("--selftest") for a in args) else "minima_compare.csv"
    pd.DataFrame(rows).to_csv(HERE / "results" / name, index=False)


if __name__ == "__main__":
    main()
