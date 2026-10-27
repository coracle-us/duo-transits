#!/usr/bin/env python3
"""Ground-photometry check of predicted transit epochs of the KIC 8462852 (Boyajian's star)
TESS single-transit companion. Run on your own machine (the cloud session can't reach these archives).

Sources (each optional; missing ones are skipped):
  --aavso FILE   CSV from AAVSO Data Download for "KIC 8462852" (you download it yourself while
                 signed in to AAVSO; nothing here signs in). Uses JD, Magnitude, Uncertainty, Band,
                 Observer Code columns. Request JD 2459400 to today.
  --asassn       pull the ASAS-SN Sky Patrol v2 light curve with pyasassn (no account needed).
  --ztf          pull ZTF public light curves from IRSA (note: the star is near or past ZTF
                 saturation, so expect few clean points; catflags==0 only).
  --tess         ask TESSCut which TESS sectors (incl. FFI-only ones after S82) cover the target.

For each epoch window below it slides a trial mid-transit time tc over the window in 0.01 d steps
and asks, per observer/band group, whether points in the transit core (|t-tc| < 0.20 d) are
>= 0.012 mag fainter than that group's own median over tc +- 15 d (excluding |t-tc| < 0.6 d).
A tc is 'excluded' if the summed chi^2 penalty for a 1.1% transit is > 25 across groups,
'consistent' if the data prefer a dip by > 25, else 'untested'. Results: epoch_coverage.csv,
one PNG per window, and a printed summary.
"""
import argparse, io, sys
import numpy as np, pandas as pd

RA, DEC = 301.564386, 44.456886          # TIC 185336364 (from the TESS FITS header)
T0 = 2458730.164                          # BJD, Madurga-Favieres et al. 2026 (MNRAS 550, stag1273)
DEPTH_MAG = 0.012                         # 1.1% in magnitudes
CORE = 0.20                               # d, half-width of the full-depth part used for the test
WINDOWS = [  # name, n, Pmin, Pmax  (P ranges: Madurga-Favieres Table 1; 776-794 d from kepler_alias_scan.py)
    ("W3 n=+2 (Oct 2026)", 2, 1301.0, 1303.6),
    ("W3 n=+1 (Mar 2023)", 1, 1301.0, 1303.6),
    ("A776 n=+1 (Oct 2021)", 1, 776.15, 776.61),
    ("A776 n=+2 (Dec 2023)", 2, 776.15, 776.61),
    ("A790 n=+1 (Nov 2021)", 1, 790.48, 790.68),
    ("A790 n=+2 (Jan 2024)", 2, 790.48, 790.68),
    ("A792 n=+1 (Nov 2021)", 1, 792.06, 793.42),
    ("A792 n=+2 (Jan 2024)", 2, 792.06, 793.42),
    ("W4 n=+1 (May 2023)", 1, 1345.640, 1345.656),
    ("W5 n=+1 (Jun 2023)", 1, 1367.3, 1368.3),
]

def load_aavso(path):
    d = pd.read_csv(path, low_memory=False)
    d = d.rename(columns={c: c.strip() for c in d.columns})
    d = d[pd.to_numeric(d["Magnitude"], errors="coerce").notna()]
    if "Fainter Than" in d: d = d[d["Fainter Than"].astype(str).str.strip() != "1"]
    out = pd.DataFrame({"jd": d["JD"].astype(float), "mag": d["Magnitude"].astype(float),
                        "err": pd.to_numeric(d.get("Uncertainty"), errors="coerce").fillna(0.01),
                        "group": "AAVSO:" + d["Observer Code"].astype(str) + ":" + d["Band"].astype(str)})
    return out[out.err < 0.05]

def load_asassn():
    from pyasassn.client import SkyPatrolClient      # PyPI package "skypatrol" (the old "pyasassn" 0.6.4 no longer works)
    # 1", not 10": a 10" cone holds ~110 ASAS-SN ids (mostly Swift-catalogue positions) whose curves would be merged in
    lcs = SkyPatrolClient().cone_search(RA, DEC, radius=1, units="arcsec", download=True)
    df = lcs.data
    ids = df["asas_sn_id"].unique()
    if len(ids) != 1: sys.exit(f"expected one ASAS-SN source within 1 arcsec, got {list(ids)}")
    print(f"ASAS-SN id {ids[0]}: {len(df)} points")
    df = df[(df["mag_err"] < 0.05) & (df["quality"] == "G")]
    return pd.DataFrame({"jd": df["jd"], "mag": df["mag"], "err": df["mag_err"],
                         "group": "ASASSN:" + df["phot_filter"].astype(str) + ":" + df["camera"].astype(str)})

def load_ztf():
    import requests
    url = ("https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
           f"?POS=CIRCLE%20{RA}%20{DEC}%200.0008&BANDNAME=g,r&FORMAT=csv&BAD_CATFLAGS_MASK=32768")
    d = pd.read_csv(io.StringIO(requests.get(url, timeout=120).text))
    d = d[d["catflags"] == 0]
    return pd.DataFrame({"jd": d["mjd"] + 2400000.5, "mag": d["mag"], "err": d["magerr"],
                         "group": "ZTF:" + d["filtercode"].astype(str)})

def tess_sectors():
    from astroquery.mast import Tesscut
    from astropy.coordinates import SkyCoord
    print(Tesscut.get_sectors(coordinates=SkyCoord(RA, DEC, unit="deg")))

def mad(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))

def test_window(data, name, n, pmin, pmax):
    """Per group, the core points are averaged into one mean and compared with the group's
    baseline median. Its error is the scatter of the group's nightly means in the baseline
    (nights = integer JD, which starts at noon UTC), divided by sqrt(nights in the core).
    Points from one night share that night's zero-point error, so treating them as
    independent (dchi2_per_point, the original statistic) overstates exclusions for
    time-series observers. With fewer than 5 baseline nights the single-point scatter is used.
    """
    tcs = np.arange(T0 + n * pmin, T0 + n * pmax + 1e-9, 0.01)
    res = []
    for tc in tcs:
        pen, pen_pp, used = 0.0, 0.0, []
        for g, d in data.groupby("group"):
            near = d[np.abs(d.jd - tc) < 15]
            core = near[np.abs(near.jd - tc) < CORE]
            base = near[np.abs(near.jd - tc) > 0.6]
            if len(core) == 0 or len(base) < 10: continue
            zp = np.median(base.mag)
            sig_pt = np.sqrt(core.err.values ** 2 + mad(base.mag.values) ** 2)
            r = core.mag.values - zp
            pen_pp += np.sum(((r - DEPTH_MAG) / sig_pt) ** 2 - (r / sig_pt) ** 2)
            nightly = base.groupby(np.floor(base.jd)).mag.mean()
            n_core = np.floor(core.jd).nunique()
            sig = (mad(nightly.values) if len(nightly) >= 5 else mad(base.mag.values)) / np.sqrt(n_core)
            sig = max(sig, 1e-3)
            rm = r.mean()
            pen += ((rm - DEPTH_MAG) / sig) ** 2 - (rm / sig) ** 2
            used.append(f"{g}({len(core)} pts, {n_core} night{'s' if n_core > 1 else ''})")
        st = "excluded" if pen > 25 else ("consistent" if pen < -25 else "untested")
        st_pp = "excluded" if pen_pp > 25 else ("consistent" if pen_pp < -25 else "untested")
        res.append((name, tc, st, pen, st_pp, pen_pp, ";".join(used)))
    return pd.DataFrame(res, columns=["window", "tc_bjd", "status", "dchi2", "status_per_point",
                                      "dchi2_per_point", "groups"])

def plot_window(data, name, n, pmin, pmax, fn):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    a, b = T0 + n * pmin - 3, T0 + n * pmax + 3
    fig, ax = plt.subplots(figsize=(9, 4))
    for g, d in data[(data.jd > a - 12) & (data.jd < b + 12)].groupby("group"):
        base = d[(d.jd < a) | (d.jd > b)]
        zp = np.median(base.mag) if len(base) else np.median(d.mag)
        w = d[(d.jd > a) & (d.jd < b)]
        if len(w): ax.errorbar(w.jd - 2460000, w.mag - zp, w.err, fmt=".", ms=3, label=g)
    ax.axvspan(T0 + n * pmin - 2460000, T0 + n * pmax - 2460000, color="r", alpha=0.1)
    ax.axhline(DEPTH_MAG, color="r", ls="--", lw=0.8)
    ax.invert_yaxis(); ax.set_xlabel("BJD - 2460000"); ax.set_ylabel("mag - group median")
    ax.set_title(name); ax.legend(fontsize=6, ncol=3)
    fig.tight_layout(); fig.savefig(fn, dpi=120); plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aavso"); ap.add_argument("--asassn", action="store_true")
    ap.add_argument("--ztf", action="store_true"); ap.add_argument("--tess", action="store_true")
    a = ap.parse_args()
    parts = []
    if a.aavso: parts.append(load_aavso(a.aavso))
    if a.asassn: parts.append(load_asassn())
    if a.ztf: parts.append(load_ztf())
    if a.tess: tess_sectors()
    if not parts: sys.exit("no photometry loaded")
    data = pd.concat(parts, ignore_index=True)
    data.to_csv("ground_photometry_merged.csv", index=False)
    out = []
    for w in WINDOWS:
        r = test_window(data, *w); out.append(r)
        plot_window(data, *w, fn="win_" + w[0].split()[0] + "_n" + str(w[1]) + ".png")
        c = r.status.value_counts(normalize=True).to_dict()
        p = r.status_per_point.value_counts(normalize=True).to_dict()
        print(f"{w[0]:24s} excluded {c.get('excluded',0):5.0%}  consistent {c.get('consistent',0):5.0%}  untested {c.get('untested',0):5.0%}"
              f"   | per-point statistic: excluded {p.get('excluded',0):4.0%} consistent {p.get('consistent',0):4.0%}")
    pd.concat(out).to_csv("epoch_coverage.csv", index=False)

if __name__ == "__main__":
    main()
