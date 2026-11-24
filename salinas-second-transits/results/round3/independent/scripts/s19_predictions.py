"""Predicted transits from today to the end of the published TESS pointing schedule, and which of them
TESS will observe (tess-point: is the star on silicon in that sector?). For TIC 107113048 it also lists
which future TESS transits would tell P from the periods that E0/E15/E36 alone allow (3P, 3P/2, 3P/5).

Times: T(E) = T0_piv + (E - Epiv) P, sigma^2 = sT0_piv^2 + ((E - Epiv) sP)^2, times the ephemeris scale.
Sector windows: halfway between tess-point's sector mid-times (actual sectors have a ~1 d gap in the
middle and start/end a little differently, so "covered" means "inside the sector window").
Run with: uv run --with tess-point python scripts/s19_predictions.py
Output: out/predictions.txt"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from astropy.time import Time
import tess_stars2px as tp
from common import TARGETS, OUT

NOW = Time.now().jd - 2457000
spd = tp.TESS_Spacecraft_Pointing_Data()
mid = spd.midtimes - 2457000
secs = np.arange(1, len(mid) + 1)
# drop special pointings appended out of time order (e.g. S135, the 3I/ATLAS pointing in Jan 2026)
keep = mid >= np.maximum.accumulate(mid)
mid, secs = mid[keep], secs[keep]
edges = np.concatenate([[mid[0] - 13.7], (mid[1:] + mid[:-1]) / 2, [mid[-1] + 13.7]])
LAST = edges[-1]
L = []
pr = L.append
pr(f"today BTJD {NOW:.2f} ({Time(NOW + 2457000, format='jd').iso[:10]}); tess-point schedule ends with "
   f"S{secs[-1]} at BTJD {LAST:.1f} ({Time(LAST + 2457000, format='jd').iso[:10]})")


def iso(b):
    return Time(b + 2457000, format="jd").iso[:16]


def sectors_on(ra, dec):
    r = tp.tess_stars2px_function_entry(1, ra, dec)
    return {int(s): (int(c), int(cc)) for s, c, cc in zip(r[3], r[4], r[5])}


def sector_at(t):
    i = np.searchsorted(edges, t) - 1
    if i >= len(secs):
        return None
    return int(secs[i]) if 0 <= i < len(secs) else None


for tic, tg in TARGETS.items():
    e = json.load(open(os.path.join(OUT, f"ephem_{tic}.json")))
    P, sP, Ep, Tp, sTp, sc = e["P"], e["sP"], e["Epiv"], e["T0_piv"], e["sT0_piv"], max(1.0, e["scale"])
    on = sectors_on(tg["ra"], tg["dec"])
    WIN = {int(s_): (edges[i], edges[i + 1]) for i, s_ in enumerate(secs)}
    on = {s_: v for s_, v in on.items() if s_ in WIN}
    fut = sorted(s for s in on if WIN[s][1] > NOW)
    pr(f"\nTIC {tic}: sectors >= S95 with the star on silicon: {sorted(s for s in on if s >= 95)}")
    pr(f"TIC {tic}: P = {P:.6f} d; future sectors with the star on silicon: "
       + (", ".join(f"S{s} (cam {on[s][0]}, {iso(WIN[s][0])[:10]} to {iso(WIN[s][1])[:10]})" for s in fut) or "none"))
    E0 = int(np.ceil((NOW - Tp) / P)) + Ep
    E1 = int(np.floor((LAST - Tp) / P)) + Ep
    rows = []
    for E in range(E0, E1 + 1):
        t = Tp + (E - Ep) * P
        s = sc * np.hypot(sTp, (E - Ep) * sP)
        sct = sector_at(t)
        rows.append((E, t, s, sct, sct in on))
    obs = [r for r in rows if r[4]]
    pr(f"  next transits ({len(rows)} before the schedule ends; {len(obs)} inside a TESS sector on silicon):")
    for E, t, s, sct, ok in rows[:4] + [r for r in rows[4:] if r[4]]:
        pr(f"    E{E:4d}: BTJD {t:.4f} +- {s*1440:5.1f} min  ({iso(t)} UTC)  "
           + (f"TESS S{sct}" if ok else ""))
    if tic == 107113048:
        pr("  alternative periods allowed by E0/E15/E36 if E29 were not a transit:")
        T0 = e["T0_E0"]
        for name, Pa in [("3P", 3 * P), ("3P/2", 1.5 * P), ("3P/5", 0.6 * P)]:
            ts = T0 + Pa * np.arange(np.ceil((NOW - T0) / Pa), np.floor((LAST - T0) / Pa) + 1)
            in_tess = [(t, sector_at(t)) for t in ts if sector_at(t) in on]
            # the true-P transit at the same time?
            txt = []
            for t, s_ in in_tess:
                n = (t - T0) / P
                txt.append(f"S{s_} BTJD {t:.2f} ({'also a P epoch' if abs(n - round(n)) * P < 0.1 else 'not a P epoch'})")
            pr(f"    {name} = {Pa:.3f} d: TESS transits {'; '.join(txt) if txt else 'none'}")
        Pt = [(t, sector_at(t)) for E, t, s, sct, ok in rows if ok]
        for t, s_ in Pt:
            n3 = [((t - T0) / (k * P)) for k in (3, 1.5, 0.6)]
            alt = [nm for nm, x in zip(("3P", "3P/2", "3P/5"), n3) if abs(x - round(x)) * P < 0.1]
            pr(f"    P transit in S{s_} at BTJD {t:.2f}: also predicted by {', '.join(alt) if alt else 'none of them'}")

txt = "\n".join(L)
open(os.path.join(OUT, "predictions.txt"), "w").write(txt + "\n")
print(txt)
