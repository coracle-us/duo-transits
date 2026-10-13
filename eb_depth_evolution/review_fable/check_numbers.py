"""
Adversarial fact-check of eb_depth_evolution/rnaas/rnaas.tex: recompute every cheap number.
Run from the repo root:  .venv/Scripts/python eb_depth_evolution/review_fable/check_numbers.py
Output: eb_depth_evolution/review_fable/check_numbers_out.txt
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.constants import G, c, M_sun, R_sun, au
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from ebmodel import flux, ecosw_for_phase  # noqa: E402

OUT = []


def p(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s)
    print(s.encode("ascii", "replace").decode())


# ---------------------------------------------------------------- 1. inclinations, rates, onsets
p("=" * 100)
p("1. INCLINATION HISTORY: ours (data/inclination_trend.csv) vs blind (independent/results.csv)")
R = pd.read_csv(ROOT / "data/inclination_trend.csv")
B = pd.read_csv(ROOT / "independent/results.csv")
base = R[(R.variant == "base") & (R.ph14 == 0.32)].iloc[0]
t = {"2014": 2014.356, "2021": 2021.804, "2023": 2023.866, "2026": 2026.041}
p("our epochs used by 18_inclination_trend.py (decimal yr):", {k: round(Time(v + 2457000, format="jd").decimalyear, 3)
                                                              for k, v in {"2014": -211.557, "2021": 2510.37, "2023": 3261.6, "2026": 4056.38}.items()})
for lab, i14, i21, i23, i26 in [("ours base", base.i2014, base.i2021, base.i2023, base.i2026),
                                 ("ours min", R.i2014.min(), R.i2021.min(), R.i2023.min(), R.i2026.min()),
                                 ("ours max", R.i2014.max(), R.i2021.max(), R.i2023.max(), R.i2026.max()),
                                 ("text", 85.7, 87.3, 88.0, 88.9),
                                 ("blind", *B.i_deg.values)]:
    r1 = (i21 - i14) / (t["2021"] - t["2014"]); r2 = (i23 - i21) / (t["2023"] - t["2021"]); r3 = (i26 - i23) / (t["2026"] - t["2023"])
    r4 = (i26 - i21) / (t["2026"] - t["2021"])
    p(f"  {lab:10s} i = {i14:.2f} {i21:.2f} {i23:.2f} {i26:.2f} | di/dt 14-21 {r1:.3f}  21-23 {r2:.3f}  23-26 {r3:.3f}  | 21-26 {r4:.3f} deg/yr")
p("  TEXT says: 0.20 (2014-2021) and 0.31 (2021-2026).  Our own base model gives 0.36 for 2021-2026; the text's own")
p("  inclinations (87.3 -> 88.9 over 4.2 yr) imply 0.38. Only the blind analysis gives 0.31.")
p(f"  ours: i2026 - i2023 = {base.i2026 - base.i2023:.3f};  blind: {B.i_deg[3] - B.i_deg[2]:.3f}  (difference {base.i2026 - base.i2023 - (B.i_deg[3] - B.i_deg[2]):.2f} deg)")
p(f"  ours variant ranges: 2023->2026 rate {R.rate_23_26.min():.3f}-{R.rate_23_26.max():.3f}; blind 0.316 +- 0.005 (variants 0.303-0.370): NO overlap of ranges")

# onset / totality with simple extrapolations
iP_ours, iS_ours = base.iP, base.iS
iP_blind, iS_blind, itot_blind = 84.415, 86.523, 89.59
p("\n  Primary onset (i = i_graze_P) by back-extrapolation from the 2014 point:")
for lab, i14, rate, ith in [("ours base, 2014-21 rate", base.i2014, base.rate_14_21, iP_ours),
                            ("ours max i2014, min rate", R.i2014.max(), R.rate_14_21.min(), R.iP.min()),
                            ("ours min i2014, max rate", R.i2014.min(), R.rate_14_21.max(), R.iP.max()),
                            ("blind nominal", B.i_deg[0], 0.197, iP_blind),
                            ("blind, TESS-era rate", B.i_deg[0], 0.308, iP_blind)]:
    p(f"    {lab:28s}: 2014.36 - ({i14:.2f}-{ith:.2f})/{rate:.3f} = {2014.356 - (i14 - ith) / rate:.1f}")
p(f"    ours global-linear fit (18_inclination_trend.py lin_tP_onset): {R.lin_tP_onset.min():.1f}-{R.lin_tP_onset.max():.1f}")
p("    TEXT: 2007-2011. OK as a union of methods, but the method is not stated.")
p("\n  Totality dates. Threshold ours 89.25 (REPORT), blind 89.59 (89.20-89.95). From each i2026 at each 2021-26 rate:")
for lab, i26, rate in [("ours base", base.i2026, 0.357), ("ours min", R.i2026.min(), 0.34), ("ours max", R.i2026.max(), 0.38),
                       ("blind", B.i_deg[3], 0.310), ("text i2026=88.9", 88.9, 0.31), ("text i2026=88.5", 88.5, 0.31)]:
    for ith in (89.2, 89.25, 89.59, 89.9):
        p(f"    {lab:16s} i2026={i26:.2f} rate={rate:.3f} -> i_tot={ith:5.2f}: {t['2026'] + (ith - i26) / rate:.1f}")
p("    TEXT: 'around 2028-2031'. With the text's own i2026 = 88.9 and i_tot = 89.2 the date is 2027.0; our model gives 2027-2028.")
p("    REPORT (ours) said 2027-29; blind 2027.5-2031. An honest union is 2027-2031.")

# ---------------------------------------------------------------- 2. depths and S/P
p("\n" + "=" * 100)
p("2. DEPTHS (results_depth_events.csv)")
T = pd.read_csv(ROOT / "results_depth_events.csv")
T["season"] = np.where(T.year < 2015, "2014", np.where(T.year < 2022, "2021", np.where(T.year < 2024, "2023", "2026")))
g = T[T.coverage > 0.8].groupby(["season", "kind"]).depth_modelfree.agg(["mean", "std", "size"])
p(g.round(4).to_string())
for s in ["2021", "2023"]:
    p(f"  S/P {s}: {g.loc[(s, 'S'), 'mean'] / g.loc[(s, 'P'), 'mean']:.3f}")
p("  TEXT: 20.6+-0.1, 27.4+-0.2, 35.6; 4.4+-0.1, 12.5+-0.1; S/P 0.21 -> 0.46. OK. The +- are eclipse-to-eclipse std only.")
p("  Reduction systematics (QLP vs own aperture) are 1-2% absolute (blind: 'ap/red' 0.016-0.022), not quoted in the text.")
p("  Blind depths (model, r=2 px): 20.1/27.0/34.5 and 4.3/12.0 -> S/P 0.216 / 0.446.")

# ---------------------------------------------------------------- 3. ground limits
p("\n" + "=" * 100)
p("3. GROUND LIMITS (data/ground_combined_scan.csv), offsets -0.6..+0.3 d")
Gc = pd.read_csv(ROOT / "data/ground_combined_scan.csv")
for k in ["P", "S"]:
    h = Gc[(Gc.kind == k) & (Gc.delta_d >= -0.6) & (Gc.delta_d <= 0.3)]
    lim = h.depth + 3 * h.err
    p(f"  {k}: coverage {h.err.notna().mean() * 100:.0f}%, worst 3-sigma limit {np.nanmax(lim):.4f}, median {np.nanmedian(lim):.4f}, max fitted depth {np.nanmax(h.depth):.4f}")
p("  TEXT: < 4.5% (3 sigma), typically < 1.6%. OK.")

# ---------------------------------------------------------------- 4. the 'turnaround' argument
p("\n" + "=" * 100)
p("4. TURNAROUND ARGUMENT: V-band depth at the quadratic i(2006.9); inclination allowed by the ground limits")
C = pd.read_csv(ROOT / "data/depth_vs_inclination.csv")
p(f"  quad_i2006.9 over variants: {R['quad_i2006.9'].min():.2f}-{R['quad_i2006.9'].max():.2f} deg")
p(f"  -> V-band primary depth: {np.interp(R['quad_i2006.9'].min(), C.inc, C.dP_V):.3f}-{np.interp(R['quad_i2006.9'].max(), C.inc, C.dP_V):.3f}")
p(f"  lin_i2006.9: {R['lin_i2006.9'].min():.2f}-{R['lin_i2006.9'].max():.2f} -> depth {np.interp(R['lin_i2006.9'].min(), C.inc, C.dP_V):.3f}")
for d in (0.045, 0.016):
    p(f"  inclination at which the V-band primary depth = {d * 100:.1f}%: {np.interp(d, C.dP_V, C.inc):.2f} deg")
p(f"  V-band depth at i2014 = 85.69-85.86: {np.interp(85.69, C.inc, C.dP_V):.3f}-{np.interp(85.86, C.inc, C.dP_V):.3f}")
p("  So the ground limits say i(2004-08) < ~85.5 (worst case) or < ~85.0 (typical). They exclude the specific quadratic through")
p("  the four points (minimum i ~ 85.45 near 2010). They do NOT exclude every turnaround: e.g. a minimum below ~85.0 in 2008-2012")
p("  followed by a faster rise to 85.7 by 2014 is not excluded by WASP/KELT. The text should say 'the quadratic trend', not 'a trend'.")
# quadratic through base points: minimum
ti = np.array([2014.356, 2021.808, 2023.866, 2026.041]); ii = np.array([base.i2014, base.i2021, base.i2023, base.i2026])
q = np.polyfit(ti - 2020, ii, 2); tmin = -q[1] / (2 * q[0]) + 2020
p(f"  quadratic through base points: minimum i = {np.polyval(q, tmin - 2020):.2f} at {tmin:.1f}; d2i/dt2 = {2 * q[0]:.4f} deg/yr^2; rate at 2026: {np.polyval(np.polyder(q), 6.041):.3f}")

# ---------------------------------------------------------------- 5. timing
p("\n" + "=" * 100)
p("5. TIMING: S-P intervals (data/sec_minus_pri_interval.csv, data/eclipse_times.csv)")
I = pd.read_csv(ROOT / "data/sec_minus_pri_interval.csv")
p(I.round(5).to_string(index=False))
d21 = (I.interval_d[1] - I.interval_d[0]) * 1440; d23 = (I.interval_d[3] - I.interval_d[2]) * 1440
p(f"  2021 one-orbit change: {d21:+.1f} +- {np.hypot(I.err_min[0], I.err_min[1]):.1f} min;  2023 E34->E35: {d23:+.1f} +- {np.hypot(I.err_min[2], I.err_min[3]):.1f} min")
ev = pd.read_csv(ROOT / "data/eclipse_times.csv")
s36 = ev[(ev.kind == "S") & (ev.cycle == 36)].tc.values[0]; p36 = ev[(ev.kind == "P") & (ev.cycle == 36)].tc.values[0]
p35 = ev[(ev.kind == "P") & (ev.cycle == 35)].tc.values[0]
p(f"  2023 E35->E36 using the PARTIAL S36 (excluded from timing in REPORT.txt): {((s36 - p36) - I.interval_d[3]) * 1440:+.1f} min")
p("  TEXT: 'by about 50 min within one orbit in 2023'. That number only exists if the partial S36 (egress only, two time")
p("  estimates differing by 15 min) is used; our own REPORT.txt excludes it. Without it, the 2023 one-orbit change is -30 +- 15 min.")
mean_chg = (I.interval_d[2:].mean() - I.interval_d[:2].mean()) * 1440 / (I.cycle[2:].mean() - I.cycle[:2].mean())
p(f"  mean change 2021->2023: {mean_chg:+.2f} min/orbit (text 3-4; blind -4.2). OK.")
# (P+S)/2 check
P = ev[(ev.kind == "P") & (ev.src != "K2 EVEREST")].set_index("cycle").tc; S = ev[ev.kind == "S"].set_index("cycle").tc
Pm = 21.774799
p(f"  2021: P1-P0 - Pm = {(P[1] - P[0] - Pm) * 1440:+.1f} min; S1-S0 - Pm = {(S[1] - S[0] - Pm) * 1440:+.1f} min; mean {((P[1] - P[0] + S[1] - S[0]) / 2 - Pm) * 1440:+.1f} min")
p(f"  2023: P36-P35 - Pm = {(P[36] - P[35] - Pm) * 1440:+.1f} min; S35-S34 - Pm = {(S[35] - S[34] - Pm) * 1440:+.1f} min; mean {((P[36] - P[35] + S[35] - S[34]) / 2 - Pm) * 1440:+.1f} min")
p("  TEXT: 'their mean stays on the period' - true to 0.2 min in 2021, to ~5 min in 2023.")
# linear ephemeris check
E = P.index.values.astype(float); A = np.vstack([np.ones_like(E), E]).T
coef, *_ = np.linalg.lstsq(A, P.values, rcond=None); res = P.values - A @ coef
p(f"  unweighted TESS-primary linear fit: P={coef[1]:.6f} T0={coef[0]:.5f}; rms {np.std(res) * 1440:.1f} min  (text P=21.774799, T0=2488.58926, rms 10 min)")
p(f"  secular Delta(e cos w) 2021->2023 ours: -0.25953 -> -0.26595 = {-0.26595 + 0.25953:+.4f}; blind -0.2593 -> -0.2658. Same.")
from ebmodel import sec_phase  # noqa: E402
esw = 0.2312; ec0 = -0.2595
dph = (sec_phase(np.hypot(ec0 + 1e-4, esw), np.arctan2(esw, ec0 + 1e-4)) - sec_phase(np.hypot(ec0, esw), np.arctan2(esw, ec0))) / 1e-4
p(f"  Per-orbit Delta(e cos w) for -26 min: d(phase)/d(ecosw) = {dph:.3f}; 26 min = {26 / 1440 / Pm:.5f} in phase -> Delta(ecosw) = {26 / 1440 / Pm / dph:.5f} per orbit (REPORT: 0.0013)")

# ---------------------------------------------------------------- 6. apsidal motion
p("\n" + "=" * 100)
p("6. CLASSICAL + GR APSIDAL MOTION (R/a 0.03-0.04, e 0.35, M_tot 4.3 Msun, P 21.775 d)")
Pd = 21.774799 * u.d; M = 4.3 * M_sun; e = 0.35
n = (2 * np.pi / Pd).to(1 / u.s).value  # rad/s
GM = (G * M).to(u.m ** 3 / u.s ** 2).value; cc = c.to(u.m / u.s).value
wdot_gr = 3 * GM ** (2 / 3) * n ** (5 / 3) / (cc ** 2 * (1 - e ** 2))  # rad/s
p(f"  GR: {np.degrees(wdot_gr) * 365.25 * 86400:.4f} deg/yr  (REPORT: 0.003 deg/yr)")
k2 = 0.005  # apsidal constant, A stars (log k2 ~ -2.3)
fe = (1 + 1.5 * e ** 2 + e ** 4 / 8) / (1 - e ** 2) ** 5
for ra in (0.03, 0.04):
    tidal = 2 * 15 * k2 * ra ** 5 * 1.0 * fe            # fraction of 2pi per orbit, two equal stars, q=1
    Rs = ra * (0.248 * au).to(u.m).value        # a ~ 53 Rsun
    om_n = 1e5 / Rs / n                           # v_rot = 100 km/s
    rot = 2 * k2 * ra ** 5 * om_n ** 2 * 2 / (1 - e ** 2) ** 2
    per_yr = 365.25 / 21.774799
    p(f"  R/a={ra}: tidal {tidal * 360 * per_yr:.2e} deg/yr; rotational (v=100 km/s, Omega/n={om_n:.0f}) {rot * 360 * per_yr:.2e} deg/yr")
p("  Observed: secular dw/dt 0.35 (ours) to 0.9 (blind) deg/yr between seasons; within-season 0.3-0.6 deg per ORBIT (5-10 deg/yr).")
p("  'Orders of magnitude too slow' is correct (factor >100). The text never states the observed secular apsidal rate.")

# ---------------------------------------------------------------- 7. nodal precession -> outer period
p("\n" + "=" * 100)
p("7. NODAL PRECESSION ESTIMATE: P_nodal ~ (4/3)(M/m3)(P2^2/P1);  di/dt_max ~ 2 pi i_mut / P_nodal")
P1 = 21.774799 / 365.25
for didt in (0.2, 0.31, 0.4):
    for im in (2, 5, 10):
        Pn = 2 * np.pi * np.radians(im) / np.radians(didt)
        for q3 in (0.05, 0.1, 0.3):
            P2 = np.sqrt(Pn * P1 * 0.75 * q3) * 365.25
            p(f"  di/dt={didt} deg/yr, i_mut={im:2d} deg -> P_nodal={Pn:6.0f} yr; m3/M={q3}: P2={P2:5.0f} d")
p("  TEXT '~100-600 d' is reproduced for i_mut 2-10 deg, m3/M 0.05-0.3. Note i_mut >= ~2.5 deg is required by the >= 4.8 deg swing")
p("  (from < 84.4 in ~2008 to > 89.2 by ~2030): the amplitude of i(t) is 2 i_mut (for a nearly edge-on outer orbit).")
# stability
p(f"  Mardling-Aarseth (e_out=0, q_out~0.1, i_mut small): a2/a1 > 2.8*(1.1)^0.4 = {2.8 * 1.1 ** 0.4:.2f} -> P2/P1 > {(2.8 * 1.1 ** 0.4) ** 1.5 / np.sqrt(1.1):.1f} -> P2 > {(2.8 * 1.1 ** 0.4) ** 1.5 / np.sqrt(1.1) * 21.775:.0f} d")
# dynamical ETV amplitude
for q3, P2 in [(0.1, 100), (0.1, 200), (0.3, 300), (0.05, 100)]:
    A = (1 / (2 * np.pi)) * q3 * 21.775 ** 2 / P2 * 1440
    p(f"  dynamical ETV amplitude ~ (1/2pi)(m3/M)(P1^2/P2): m3/M={q3}, P2={P2} d -> {A:.0f} min (x f(e) ~ 1-3)")
p("  A 26-min change of T_S-T_P in ONE orbit, vs 3-4 min/orbit averaged over 2 yr, needs the perturbation to vary on a few-orbit")
p("  timescale, i.e. P2 of order 100-200 d, the LOW end of the quoted range. Worth saying.")

# ---------------------------------------------------------------- 8. Gaia RUWE argument
p("\n" + "=" * 100)
p("8. GAIA: photocentre wobble of the inner binary around the triple's barycentre, d = 661 pc")
d_pc = 661.0
for m3 in (0.1, 0.2, 0.3, 0.5):
    for P2 in (100, 200, 400, 600):
        Mt = 4.3 + m3; a_au = (Mt * (P2 / 365.25) ** 2) ** (1 / 3)
        a_cm = a_au * m3 / Mt  # tertiary light negligible
        p(f"  m3={m3} Msun, P2={P2} d: a_out={a_au:.2f} AU, binary wobble semi-amplitude {a_cm / d_pc * 1000:.3f} mas")
p("  Gaia DR3 per-transit AL precision at G=11 is ~0.1-0.2 mas; wobbles >~0.2 mas (0.5 Msun at P2 >~ 200 d) would raise RUWE.")
p("  RUWE 1.007 therefore disfavours m3 >~ 0.4 Msun for P2 >~ 200 d; it says nothing about m3 <~ 0.2 Msun at any P2 or any mass at P2 ~ 100 d.")
p("  'Mildly favours a low-mass tertiary' is fair; 'or a short outer period' should be added.")

# ---------------------------------------------------------------- 9. mass / radius sanity
p("\n" + "=" * 100)
p("9. MASS/RADIUS: R1/a=0.0297, R2/a=0.0391, G=10.95, A_G=0.37, plx -> 661 pc, Teff 7400-7700 K")
MG = 10.95 - 5 * np.log10(d_pc / 10) - 0.37
p(f"  M_G(system) = {MG:.2f} (REPORT 1.49)")
for Teff, BC in [(7400, 0.08), (7700, 0.06)]:
    Mbol = MG + BC; L = 10 ** ((4.74 - Mbol) / 2.5)
    R2sum = L / (Teff / 5772) ** 4  # R1^2 + R2^2 in Rsun^2 (same Teff)
    R1 = np.sqrt(R2sum / (1 + 1.32 ** 2)); R2 = 1.32 * R1
    a = R1 / 0.0297; Mtot = (a / 215.03) ** 3 / (21.774799 / 365.25) ** 2
    p(f"  Teff={Teff}: L={L:.1f} Lsun, R1={R1:.2f}, R2={R2:.2f} Rsun, a={a:.1f} Rsun, M_tot={Mtot:.2f} Msun")
p("  TEXT: 1.6 and 2.1 Rsun, 4.3 Msun. Reproduced within the Teff/BC/extinction slop (R1 1.5-1.6, M_tot 3.5-4.4).")
p("  Note M_tot is derived from the radii and R/a, so 'consistent with two A stars' is a self-consistency check, not independent.")
p("  Blind radii (R1/a 0.032, R2/a 0.037) would give R1~R2~1.8 Rsun and a similar M_tot.")

# ---------------------------------------------------------------- 10. predictions
p("\n" + "=" * 100)
p("10. PREDICTED PRIMARIES (P=21.774799, T0=BTJD 2488.58926), geocentric UTC with light-travel time removed")
COORD = SkyCoord(97.98012 * u.deg, 25.48207 * u.deg); GEO = EarthLocation.from_geocentric(0, 0, 0, unit=u.m)
for ncyc in (85, 86, 87, 88):
    tb = Time(2457000 + 2488.58926 + ncyc * 21.774799, format="jd", scale="tdb", location=GEO)
    lt = tb.light_travel_time(COORD, kind="barycentric")
    utc = (tb - lt).utc
    p(f"  E={ncyc}: BJD_TDB {tb.jd:.4f}  ltt {lt.to(u.min).value:+.1f} min  UTC {utc.iso[:19]}  (text: see draft)")
p("  Text: Oct 25 22:39, Nov 16 17:12, Dec 8 11:46, Dec 30 06:21. iso[:16] truncates seconds (fine at +-30 min).")
# ephemeris extrapolation error
chi2, dof = 483, 4
sigP = 0.000013 * np.sqrt(chi2 / dof)
p(f"  ephemeris error at E=88 (16 cycles past E72): ~{sigP * 16 * 1440:.1f} min (inflated), ETV +-20 min -> +-30 min OK")

# ---------------------------------------------------------------- 11. eclipse duration
p("\n" + "=" * 100)
p("11. PRIMARY ECLIPSE DURATION from the base model (first to last contact)")
d = json.load(open(ROOT / "data/physical_fit_base.json")); pp = d["params"]
ec = ecosw_for_phase(0.32817, pp["esinw"]); e_ = np.hypot(ec, pp["esinw"]); w_ = np.arctan2(pp["esinw"], ec)
tt = np.linspace(-0.4, 0.4, 8001)
for inc in (87.26, 87.91, 88.77, 89.3, 89.6, 90.0):
    f = flux(tt, 0, "P", 21.774799, e_, w_, np.radians(inc), pp["r1"], pp["r2"], pp["J"])
    inside = tt[f < 1 - 1e-6]
    p(f"  i={inc:5.2f}: duration {(inside.max() - inside.min()) * 24:.2f} h, depth {(1 - f.min()) * 100:.1f}%")
p("  TEXT: 'Each lasts about 7 h'. The model gives ~8 h now and ~8.5 h near totality; blind 'duration' 7.6-8.4 h. Say 'about 8 h'.")

# ---------------------------------------------------------------- 12. K2 secondary limit
p("\n" + "=" * 100)
p("12. K2 2014 SECONDARY LIMIT: ours < 0.03% (3 sigma); blind < 0.08% (99% false-alarm of a box scan); K2 note (AAS82225) says < 0.1%.")
p("    Three different numbers for the same null. The blind's trials-corrected limit is the defensible one; < 0.1% is consistent with AAS82225.")

# ---------------------------------------------------------------- 13. DASCH
p("\n" + "=" * 100)
p("13. DASCH: 1206 points, eclipse duty cycle ~8 h / 21.77 d = {:.1%} -> {:.0f} points expected in a primary eclipse;".format(8 / 24 / 21.775, 1206 * 8 / 24 / 21.775))
p("    8/1206 > 0.3 mag vs 0.64% in comparisons (7.7 expected). A 35% eclipse (0.47 mag) persisting would add ~10-15 faint points. OK.")

Path(__file__).with_name("check_numbers_out.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
