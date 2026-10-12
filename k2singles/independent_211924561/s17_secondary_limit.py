"""Secondary-eclipse limit at phase 0.5 for P=70.247 d (circular): weighted mean template amplitude
over the observed phase-0.5 epochs (= odd epochs of the P/2 grid), using s14's fit_amp and empirical sigmas."""
import json, os, numpy as np
import s14_alias_analysis as A
P = 70.24701; T0 = 179.89774
rows = []
for name, c in A.cache.items():
    n0 = int(np.ceil((c["tmin"] - T0) / P - 0.5)) - 1; n1 = int(np.floor((c["tmax"] - T0) / P)) + 1
    for n in range(n0, n1 + 1):
        tc = T0 + (n + 0.5) * P
        a, cv = A.fit_amp(c["t"], c["f"], tc, c["exp"], c["u"])
        if np.isfinite(a) and cv >= 0.8:
            depth_full = 1e3 * (1 - A.tfit.model(np.array([0.]), 0, A.k_, A.b_, A.aR_, c["u"])[0])
            rows.append(dict(ds=name, tc=round(tc, 3), A=float(a), sA=float(c["sigA"]), ppt=float(a * depth_full), sppt=float(c["sigA"] * depth_full)))
# one reduction per epoch-family for independence: pick smallest sigma per tc
best = {}
for r in rows:
    if r["tc"] not in best or r["sppt"] < best[r["tc"]]["sppt"]:
        best[r["tc"]] = r
v = list(best.values())
w = np.array([1 / r["sppt"] ** 2 for r in v]); d = np.array([r["ppt"] for r in v])
m = np.sum(w * d) / np.sum(w); s = 1 / np.sqrt(np.sum(w))
out = dict(epochs=v, mean_ppt=m, err_ppt=s, upper3sig_ppt=m + 3 * s if m > 0 else 3 * s)
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(A.HERE, "secondary_limit_P70.json"), "w"), indent=1)
