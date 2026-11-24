"""
AAS82225 revision: the TESS-SPOC (FFI) TCEs of two Table 1 stars called new, from their Sector 70 DV XML reports
(found by dv_listing.py in review 2, downloaded 2026-11-24 with Jon's OK into data/literature/dv_tic<TIC>/; not edited):
  HD 220967 = EPIC 246194159 = TIC 301328232
  EPIC 246089124 = TIC 49603732
For each TCE: its period and epoch, and the flux at each of its predicted transit times inside S70, measured in our
cached S70 light curves (TESS-SPOC PDCSAP and SAP, quality 0; QLP), as the median in a window of the TCE pulse
duration minus the median of the two flanking windows of the same width, 1-2 pulse durations away. Also whether the
TCE period fits our eclipse times (results/ephemerides/<EPIC>.json).
Output: results/literature/dv_tce_s70.txt
Usage:  python -I k2singles/revision2/dv_tce_s70.py
"""
import glob
import json
import os
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
K2 = os.path.dirname(HERE)
LC = os.path.expanduser("~/Documents/GitHub/coracle/duo-transits-cache/raw/lightcurves")
NS = "{http://www.nasa.gov/2018/TESS/DV}"
# secondary-eclipse phase: from the joint fit where the ephemeris file has one, else the table note (build_table.py)
STARS = {301328232: (246194159, "HD 220967 = EPIC 246194159", None), 49603732: (246089124, "EPIC 246089124", 0.494)}
out = []


def say(s=""):
    print(s)
    out.append(s)


def load(tic, prod):
    f = os.path.join(LC, str(tic), f"s070_{prod}_200s.csv.gz")
    d = pd.read_csv(f)
    return d[(d.quality == 0) & np.isfinite(d.flux)] if prod == "TESS-SPOC" else d[np.isfinite(d.flux)]


def dip(t, f, tc, w, num=False):
    f = f / np.nanmedian(f)
    inn = np.abs(t - tc) < w / 2
    flank = (np.abs(t - tc) > w) & (np.abs(t - tc) < 2 * w)
    if inn.sum() < 5 or flank.sum() < 5:
        return np.nan if num else f"no data (n_in {inn.sum()}, n_flank {flank.sum()})"
    if num:
        return 1e3 * (np.median(f[flank]) - np.median(f[inn]))
    return f"{1e3 * (np.median(f[flank]) - np.median(f[inn])):+8.2f} ppt  (n_in {inn.sum()})"


for tic, (epic, name, phase) in STARS.items():
    eph = json.load(open(os.path.join(K2, "results", "ephemerides", f"{epic}.json")))
    j = eph.get("joint", eph)
    P0, T0, ps = j["P"], j["T0"], j.get("phase_sec", phase)
    say(f"######## {name} (TIC {tic}); our P = {P0:.6f} d, primary T0 = {T0:.4f}, secondary phase {ps:.3f}")
    ours = {}
    for k in range(int((3200 - T0) / P0), int((3240 - T0) / P0) + 2):
        ours[T0 + k * P0] = "our primary"
        ours[T0 + (k + ps) * P0] = "our secondary"
    ev = sorted(e["T"] for e in eph["events"])
    sp = {p: load(tic, p) for p in ("TESS-SPOC", "QLP")}
    lo, hi = sp["TESS-SPOC"].time.min(), sp["TESS-SPOC"].time.max()
    say(f"S70 TESS-SPOC quality-0 data span {lo:.2f}-{hi:.2f}")
    for f in sorted(glob.glob(os.path.join(HERE, "data", "literature", f"dv_tic{tic}", "*_dvr.xml"))):
        for c in ET.parse(f).getroot().iter(NS + "planetCandidate"):
            P, E = float(c.get("orbitalPeriodInDays")), float(c.get("epochTjd"))
            w = float(c.get("trialTransitPulseDurationInHours")) / 24
            say(f"\nTCE {c.get('planetNumber')}: P = {P:.6f} d, epoch {E:.4f}, pulse {24 * w:.1f} h, "
                f"MES {float(c.get('maxMultipleEventSigma')):.1f}, max SES {float(c.get('maxSingleEventSigma')):.1f}, "
                f"transits observed {c.get('observedTransitCount')} of {c.get('expectedTransitCount')}")
            for k in range(int(np.floor((lo - 1 - E) / P)), int(np.ceil((hi + 1 - E) / P)) + 1):
                tc = E + k * P
                lab = next((f"{v} ({tc - t:+.2f} d)" for t, v in ours.items() if abs(tc - t) < 0.3), "")
                say(f"  {tc:9.3f}  PDCSAP {dip(sp['TESS-SPOC'].time.values, sp['TESS-SPOC'].pdcsap_flux.values, tc, w)}"
                    f" | SAP {dip(sp['TESS-SPOC'].time.values, sp['TESS-SPOC'].sap_flux.values, tc, w)}"
                    f" | QLP {dip(sp['QLP'].time.values, sp['QLP'].flux.values, tc, w)}  {lab}")
            # control: the same statistic every 0.1 d in S70, at least 1 d from our eclipses
            grid = [g for g in np.arange(lo, hi, 0.1) if min(abs(g - t) for t in ours) > 1]
            ctl = np.array([dip(sp["TESS-SPOC"].time.values, sp["TESS-SPOC"].pdcsap_flux.values, g, w, True)
                            for g in grid])
            ctl = ctl[np.isfinite(ctl)]
            say(f"  control (PDCSAP, {len(ctl)} times 0.1 d apart, >1 d from our eclipses): median {np.median(ctl):+.2f}"
                f" ppt, 95th percentile {np.percentile(ctl, 95):+.2f} ppt, max {ctl.max():+.2f} ppt")
            res = [((t - E) / P - round((t - E) / P)) * P for t in ev]
            say("  our eclipse times, offset from the TCE ephemeris (d): " + ", ".join(f"{x:+.2f}" for x in res))
            say(f"  our period / TCE period = {P0 / P:.4f}")
    say()
open(os.path.join(HERE, "results", "literature", "dv_tce_s70.txt"), "w").write("\n".join(out) + "\n")
