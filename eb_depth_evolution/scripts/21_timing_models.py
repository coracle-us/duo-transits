"""Timing analysis. Eclipse times from data/eclipse_times.csv (own TESScut r=1.5 + K2 EVEREST). The partial S72 secondary
(egress only, KvW and model times differ by 15 min) is excluded from timing.
 (1) TESS primaries, linear: chi2, rms.  (2) primaries incl. K2: linear vs quadratic with a jitter term fitted so chi2/dof = 1.
 (3) Secondary-minus-primary interval per cycle vs time (dilution- and LTTE-independent; tracks e*cos(w)).
 (4) Predicted primary times for the 2026-27 season."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.time import Time
ev = pd.read_csv("data/eclipse_times.csv")
ev = ev[~((ev.kind == "S") & (ev.cycle == 36))].copy()
P = ev[ev.kind == "P"].sort_values("cycle"); S = ev[ev.kind == "S"].sort_values("cycle")


def wls(E, T, s, deg):
    A = np.vander(E, deg + 1, increasing=True); W = 1 / s ** 2
    C = np.linalg.inv(A.T @ (A * W[:, None])); p = C @ (A.T @ (W * T)); r = T - A @ p
    return p, C, r, float(np.sum(r ** 2 * W))


def with_jitter(E, T, s, deg):
    lo, hi = 0.0, 0.1
    for _ in range(60):
        j = 0.5 * (lo + hi); p, C, r, chi2 = wls(E, T, np.hypot(s, j), deg)
        if chi2 / (len(E) - deg - 1) > 1: lo = j
        else: hi = j
    return j, p, C, r


Et = P[P.src != "K2 EVEREST"]
p, C, r, chi2 = wls(Et.cycle.values * 1.0, Et.tc.values, Et.e_tc_floor.values, 1)
print(f"(1) TESS primaries linear: P={p[1]:.6f}+-{np.sqrt(C[1,1]):.6f} T0={p[0]:.5f}; chi2={chi2:.0f} dof={len(Et)-2}; rms O-C={np.std(r)*1440:.1f} min")
E, T, s = P.cycle.values * 1.0, P.tc.values, P.e_tc_floor.values
for deg in (1, 2):
    j, pp, CC, rr = with_jitter(E, T, s, deg)
    print(f"(2) primaries incl. K2, degree {deg}: jitter needed {j*1440:.1f} min; coefficients {pp}; errors {np.sqrt(np.diag(CC))}")
    if deg == 2:
        Pm = pp[1]; dPdE = 2 * pp[2]
        print(f"    quadratic term c2={pp[2]:.3e}+-{np.sqrt(CC[2,2]):.3e} d ({pp[2]/np.sqrt(CC[2,2]):.1f} sigma); dP/dt = {dPdE/Pm*365.25*86400:.2f} s/yr")
# K2 offset w.r.t. TESS-only linear ephemeris, with ephemeris error inflated by sqrt(chi2/dof)
k2 = P[P.src == "K2 EVEREST"].iloc[0]; Ek = k2.cycle
pred = p[0] + p[1] * Ek; sp = np.sqrt(C[0, 0] + Ek ** 2 * C[1, 1] + 2 * Ek * C[0, 1]); infl = np.sqrt(chi2 / (len(Et) - 2))
print(f"    K2 primary O-C vs TESS linear: {(k2.tc-pred)*1440:+.1f} min; ephemeris error {sp*1440:.1f} min (x{infl:.1f} for scatter -> {sp*infl*1440:.1f} min): "
      f"{(k2.tc-pred)/(sp*infl):.1f} sigma")
# (3) S - P interval per cycle
Pfit = np.poly1d([p[1], p[0]])
rows = []
for _, s_ in S.iterrows():
    same = P[P.cycle == s_.cycle]
    if len(same):
        tp = same.tc.values[0]; how = "observed P same cycle"; etp = same.e_tc_floor.values[0]
    else:
        nxt = P[P.cycle == s_.cycle + 1]
        tp = nxt.tc.values[0] - p[1]; how = "P(next cycle) - P_mean"; etp = np.hypot(nxt.e_tc_floor.values[0], 15 / 1440)
    rows.append(dict(cycle=s_.cycle, year=Time(s_.tc + 2457000, format="jd").decimalyear, interval_d=s_.tc - tp,
                     err_min=np.hypot(s_.e_tc_floor, etp) * 1440, ref=how))
I = pd.DataFrame(rows); print("(3) secondary-minus-primary interval:"); print(I.round(5).to_string(index=False))
d1 = (I.interval_d.values[1] - I.interval_d.values[0]) * 1440; d2 = (I.interval_d.values[3] - I.interval_d.values[2]) * 1440
dl = (I.interval_d.values[2:].mean() - I.interval_d.values[:2].mean()) * 1440 / (I.cycle.values[2:].mean() - I.cycle.values[:2].mean())
print(f"    change within 2021 (1 cycle): {d1:+.1f} min; within 2023 (1 cycle): {d2:+.1f} min; mean 2021->2023: {dl:+.2f} min/cycle "
      f"(total {(I.interval_d.values[2:].mean()-I.interval_d.values[:2].mean())*1440:+.0f} min)")
I.to_csv("data/sec_minus_pri_interval.csv", index=False)
# figure
fig, ax = plt.subplots(2, 1, figsize=(10, 8))
for kind, g, col in [("primary", P, "b"), ("secondary", S, "g")]:
    ref = p[0] + p[1] * g.cycle + (0.0 if kind == "primary" else 0.3302 * p[1])
    ax[0].errorbar(g.cycle, (g.tc - ref) * 1440, g.e_tc_floor * 1440, fmt="o", color=col, label=kind + (" (secondary offset by const)" if kind == "secondary" else ""))
ax[0].set_xlabel("cycle (TESS linear ephemeris P=%.6f d)" % p[1]); ax[0].set_ylabel("O-C (min)"); ax[0].axhline(0, color="k", lw=0.5); ax[0].legend()
ax[0].set_title("EPIC 202126877 eclipse timing; K2 2014 at cycle -124")
ax[1].errorbar(I.year, (I.interval_d - I.interval_d.mean()) * 1440, I.err_min, fmt="o")
ax[1].set_xlabel("year"); ax[1].set_ylabel("T_S - T_P - mean (min)"); ax[1].set_title("secondary-minus-primary interval (tracks e cos w)")
plt.tight_layout(); plt.savefig("figures/timing_oc.png", dpi=90)
# (4) predictions for the 2026-27 observing season (RA 6.5h: Oct-Mar). Uncertainty: jitter (~12 min) + ephemeris
print("(4) predicted primary minima, 2026-10 .. 2027-04 (TESS-only linear ephemeris; expect +-30 min from ETVs):")
for n in range(int((Time("2026-10-01").jd - 2457000 - p[0]) / p[1]) + 1, int((Time("2027-04-30").jd - 2457000 - p[0]) / p[1]) + 1):
    t = p[0] + p[1] * n; print(f"    E={n}  BTJD {t:.4f}  {Time(t + 2457000, format='jd').iso[:16]} UTC")
