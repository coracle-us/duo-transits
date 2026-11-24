"""
AAS82225 revision: the SPOC TCEs of the two Table 1 stars whose only prior periods are SPOC TCEs, from the DV XML
reports on MAST (downloaded 2026-11-23 with Jon's OK into data/literature/dv_tic<TIC>/; the files are not edited):
  TOI-5560 = EPIC 211924561 = TIC 238606535 (runs S42-46, S72, S2-S72)
  EPIC 211489484 = TIC 6951741 (run S1-S46)
For each TCE: its epoch and period, which of its transit times fall on our events or on events we rejected, and
whether the TCE period fits all our primary events. Our events: independent_211924561/README.md and
s06_vet_tess_events.py (TOI-5560); independent_table/scripts/s7_final.py (211489484).
Usage:  python -I k2singles/revision2/dv_tce_check.py
"""
import glob
import os
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
NS = "{http://www.nasa.gov/2018/TESS/DV}"
STARS = {
    238606535: dict(
        name="TOI-5560 = EPIC 211924561", period=140.49401,
        ours={179.89774: "C5 transit", 1163.3554: "C16 transit", 2568.29833: "S46 transit", 3270.76459: "S72 transit"},
        other={2509.81: "S44 rejected (scattered-light step)", 2539.52: "S45 rejected (QLP only, gap edge)",
               2541.05: "S45 vetted, not a transit (< 0.5 ppt, on a steep background ramp)",
               3275.46: "S72 rejected (TESS-SPOC PDC edge)"},
        span=(2491.0, 3290.0),
        runs={("s0042", "s0046"): (2500.19, 2578.70), ("s0072", "s0072"): (3260.18, 3285.59)},  # S44-46, S72 data
        starts={3260.18: "the first S72 cadence"}),
    6951741: dict(
        name="EPIC 211489484", period=82.25843,
        ours={192.621: "C5 primary", 1261.983: "C18 primary", 1508.757: "S7 primary", 2249.087: "S34 primary",
              2578.118: "S46 primary", 3236.180: "S71 primary"},
        other={208.605: "C5 secondary", 1277.968: "C18 secondary", 3252.17: "S71 secondary"},
        span=(1325.0, 3290.0), runs={}, starts={}),
}


def name(t, s, tol=0.3):
    for d in (s["ours"], s["other"]):
        for k, v in d.items():
            if abs(t - k) < tol:
                return f"{v} ({t - k:+.2f} d)"
    for k, v in s["starts"].items():
        if abs(t - k) < 0.5:
            return f"{t - k:+.2f} d from {v}; not an event in our search"
    return "not an event in our search"


for tic, s in STARS.items():
    print(f"\n######## {s['name']} (TIC {tic}); our period {s['period']} d")
    for f in sorted(glob.glob(os.path.join(HERE, "data", "literature", f"dv_tic{tic}", "*_dvr.xml"))):
        run = tuple(os.path.basename(f).split("-")[1:3])
        for c in ET.parse(f).getroot().iter(NS + "planetCandidate"):
            P, T0 = float(c.get("orbitalPeriodInDays")), float(c.get("epochTjd"))
            print(f"\nSPOC run {run[0]}-{run[1]}, TCE {c.get('planetNumber')}: P = {P:.6f} d, epoch {T0:.4f}, "
                  f"MES {float(c.get('maxMultipleEventSigma')):.1f}, transits observed {c.get('observedTransitCount')} "
                  f"of {c.get('expectedTransitCount')} expected, pulse {c.get('trialTransitPulseDurationInHours')} h")
            if P > 1:
                lo, hi = s["runs"].get(run, s["span"])
                k0, k1 = int((s["span"][0] - T0) // P), int((s["span"][1] - T0) // P) + 1
                for t in (T0 + k * P for k in range(k0, k1 + 1)):
                    n = name(t, s)
                    if (run in s["runs"] and lo <= t <= hi) or n != "not an event in our search":
                        print(f"   {t:9.3f}  {n}")
            if P < 1:
                # a sub-day period cannot be the eclipse period: every ~P of data would hold an event, and our
                # searches find isolated events only. Just record how it relates to our period.
                print(f"   sub-day period; our period / TCE period = {s['period'] / P:.4f}")
                continue
            res = {v: ((t - T0) / P - round((t - T0) / P)) * P for t, v in s["ours"].items()}
            worst = max(abs(x) for x in res.values())
            print("   our events, offset from the TCE ephemeris (d): "
                  + ", ".join(f"{v.split()[0]} {x:+.3f}" for v, x in res.items())
                  + f"  -> {'fits' if worst < 0.1 else 'does not fit'} all of them")
            r = s["period"] / P
            print(f"   our period / TCE period = {r:.4f}")
