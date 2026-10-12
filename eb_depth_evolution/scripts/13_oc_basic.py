"""Eclipse times -> O-C.  Primary and secondary minima from own TESScut r=1.5 photometry (cross-checked with SPOC/QLP), K2 primary
from EVEREST.  Fits: (1) linear ephemeris to TESS primaries; (2) same for secondaries; (3) secondary-minus-primary phase vs time."""
import numpy as np, pandas as pd
ev = pd.read_csv("data/tess_events_all_products.csv")
own = ev[(ev["product"] == "own r=1.5") & (ev.covers_min)].copy()
own = own[~((own.kind == "S") & (own.cycle == -1))]          # S43 cycle -1 secondary: only the egress tail at sector start
k2 = pd.read_csv("data/k2_primary_fits.csv").set_index("product")
rows = own[["kind", "cycle", "tc", "e_tc", "depth", "e_depth", "w", "sector", "coverage"]].copy()
rows["src"] = "TESS own r=1.5"
rows = pd.concat([rows, pd.DataFrame([dict(kind="P", cycle=-124, tc=k2.loc["EVEREST", "tc"], e_tc=k2.loc["EVEREST", "e_tc"],
                    depth=k2.loc["EVEREST", "depth"], e_depth=k2.loc["EVEREST", "e_depth"], w=k2.loc["EVEREST", "w"], sector=0, coverage=1.0, src="K2 EVEREST")])])
# the S72 secondary is partial (egress only, coverage 0.62): keep but flag
rows = rows.sort_values("tc").reset_index(drop=True)
rows["e_tc_floor"] = np.sqrt(rows.e_tc**2 + (1.0/1440)**2)
rows.to_csv("data/eclipse_times.csv", index=False)
def linfit(E, T, s):
    A = np.vstack([np.ones_like(E), E]).T; W = 1/s**2
    C = np.linalg.inv(A.T @ (A*W[:, None])); p = C @ (A.T @ (W*T)); r = T - A @ p
    return p, C, r, np.sum(r**2*W)
P = rows[(rows.kind == "P") & (rows.src != "K2 EVEREST")]
p, C, r, chi2 = linfit(P.cycle.values.astype(float), P.tc.values, P.e_tc_floor.values)
print(f"TESS primaries linear: T0={p[0]:.5f}+-{np.sqrt(C[0,0]):.5f}  P={p[1]:.6f}+-{np.sqrt(C[1,1]):.6f}  chi2={chi2:.1f} dof={len(P)-2}")
for c, t, rr, s in zip(P.cycle, P.tc, r, P.e_tc_floor): print(f"   E={c:4d} T={t:.5f} O-C={rr*1440:+7.1f} +- {s*1440:.1f} min")
pred = p[0] + p[1]*(-124); sp = np.sqrt(C[0,0] + 124**2*C[1,1] - 2*124*C[0,1])
k = rows[rows.src == "K2 EVEREST"].iloc[0]
print(f"K2 primary O-C vs TESS-primary ephemeris: {(k.tc - pred)*1440:+.1f} min, ephemeris error {sp*1440:.1f} min -> {(k.tc-pred)/np.hypot(sp, k.e_tc):.1f} sigma")
S = rows[rows.kind == "S"]
ps, Cs, rs, chi2s = linfit(S.cycle.values.astype(float), S.tc.values, S.e_tc_floor.values)
print(f"TESS secondaries linear: T0={ps[0]:.5f} P={ps[1]:.6f}+-{np.sqrt(Cs[1,1]):.6f} chi2={chi2s:.1f} dof={len(S)-2}")
for c, t, rr, s in zip(S.cycle, S.tc, rs, S.e_tc_floor): print(f"   E={c:4d} T={t:.5f} O-C={rr*1440:+7.1f} +- {s*1440:.1f} min")
print(f"Period difference P_sec - P_pri = {(ps[1]-p[1])*1440:.2f} +- {np.sqrt(Cs[1,1]+C[1,1])*1440:.2f} min/cycle")
# secondary phase per epoch, using nearest primaries in same season
for yr, sel in [("2021", (S.cycle <= 2)), ("2023", (S.cycle >= 34))]:
    ss = S[sel]; ph = []
    for _, s_ in ss.iterrows():
        tp = p[0] + p[1]*s_.cycle + (r[np.argmin(np.abs(P.cycle.values - s_.cycle))])   # observed primary nearest in cycle
        ph.append((s_.tc - tp)/p[1])
    print(f"{yr}: secondary phase (T_S - T_P(observed nearest primary))/P = {np.mean(ph):.5f} (individual {np.round(ph,5)})")
