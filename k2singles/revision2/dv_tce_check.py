"""
AAS82225 revision: the SPOC TCEs of TOI-5560 (EPIC 211924561 = TIC 238606535), from the DV XML reports on MAST
(downloaded 2026-11-23 with Jon's OK into data/literature/dv_tic238606535/; the files are not edited).
For each TCE: its epoch and period, which of its transit times fall on our transits or on events we rejected
(independent_211924561/README.md and s06_vet_tess_events.py), and whether the TCE period fits our four transits.
Usage:  python -I k2singles/revision2/dv_tce_check.py
"""
import glob
import os
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
DV = os.path.join(HERE, "data", "literature", "dv_tic238606535")
OURS = {179.89774: "C5 transit", 1163.3554: "C16 transit", 2568.29833: "S46 transit", 3270.76459: "S72 transit"}
REJECTED = {2509.81: "S44 rejected (scattered-light step)", 2539.52: "S45 rejected (QLP only, gap edge)",
            2541.05: "S45 vetted, not a transit (< 0.5 ppt, on a steep background ramp)",
            3275.46: "S72 rejected (TESS-SPOC PDC edge)"}
TESS = (2491.0, 3290.0)  # S44 start to S72 end
SPANS = {("s0042", "s0046"): (2419.0, 2583.0), ("s0072", "s0072"): (3260.18, 3285.59)}  # each run's data, approximately
S72_START = 3260.18      # first S72 cadence (independent_211924561/inventory.txt)


def name(t, tol=0.3):
    for d in (OURS, REJECTED):
        for k, v in d.items():
            if abs(t - k) < tol:
                return f"{v} ({t - k:+.2f} d)"
    if abs(t - S72_START) < 0.5:
        return f"{t - S72_START:+.2f} d from the first S72 cadence; not an event in our search"
    return "not an event in our search"


for f in sorted(glob.glob(os.path.join(DV, "*_dvr.xml"))):
    run = os.path.basename(f).split("-")[1:3]
    for c in ET.parse(f).getroot().iter("{http://www.nasa.gov/2018/TESS/DV}planetCandidate"):
        P, T0 = float(c.get("orbitalPeriodInDays")), float(c.get("epochTjd"))
        print(f"\nSPOC run {run[0]}-{run[1]}, TCE {c.get('planetNumber')}: P = {P:.6f} d, epoch {T0:.4f}, "
              f"MES {float(c.get('maxMultipleEventSigma')):.1f}, transits observed {c.get('observedTransitCount')} "
              f"of {c.get('expectedTransitCount')} expected")
        k0, k1 = int((TESS[0] - T0) // P), int((TESS[1] - T0) // P) + 1
        hits = [T0 + k * P for k in range(k0, k1 + 1) if TESS[0] <= T0 + k * P <= TESS[1]]
        lo, hi = SPANS.get(tuple(run), TESS)  # within the run's own data: every predicted time; else only named ones
        for t in hits:
            n = name(t)
            if (lo <= t <= hi and tuple(run) in SPANS) or n != "not an event in our search":
                print(f"   {t:9.3f}  {n}")
        res = {v: ((t - T0) / P - round((t - T0) / P)) * P for t, v in OURS.items()}
        worst = max(abs(x) for x in res.values())
        print("   our transits, offset from the TCE ephemeris (d): "
              + ", ".join(f"{v.split()[0]} {x:+.2f}" for v, x in res.items())
              + f"  -> {'fits' if worst < 0.1 else 'does not fit'} all four")
