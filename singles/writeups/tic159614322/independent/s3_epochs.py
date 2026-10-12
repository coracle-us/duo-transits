"""Predict every transit epoch for n = 1..12 that falls within (or near) TESS data; measure coverage per product."""
import json, os
import numpy as np
from lcio import load, list_files, OUT

fk = json.load(open(os.path.join(OUT, "fit_known.json")))["primary"]
TC1, S1 = fk["tc1"][0], max(fk["tc1"][1:])
TC2, S2 = fk["tc2"][0], max(fk["tc2"][1:])
DT = TC2 - TC1
T14 = fk["T14_hr"][0] / 24.0
HALF = T14 / 2

sectors = sorted(set(s for s, p, c, f in list_files()))
spans = {}
for s in sectors:
    lo, hi = 1e9, -1e9
    for p in ["SPOC", "TESS-SPOC", "QLP"]:
        d = load(s, p, "flux", strict_quality=False)
        if d is not None:
            lo, hi = min(lo, d.time.min()), max(hi, d.time.max())
    spans[s] = (lo, hi)


def coverage(s, p, tc, strict=True):
    d = load(s, p, "flux", strict_quality=strict)
    if d is None:
        return None
    t = d.time.values
    cad = d.cad.iloc[0] / 86400
    phase = t - tc
    intr = np.abs(phase) < HALF
    # fraction of in-transit window covered, plus pre/post baseline
    grid = np.arange(-HALF, HALF, cad)
    cov_grid = np.array([np.any(np.abs(phase - g) < cad * 0.75) for g in grid])
    ing = np.abs(phase + HALF) < 0.03  # within ~45 min of 1st contact
    egr = np.abs(phase - HALF) < 0.03
    pre = (phase < -HALF) & (phase > -HALF - 0.4)
    post = (phase > HALF) & (phase < HALF + 0.4)
    first_half = (phase > -HALF) & (phase < 0)
    second_half = (phase >= 0) & (phase < HALF)
    return dict(n_in=int(intr.sum()), frac_in=float(cov_grid.mean()), ingress=bool(ing.any()), egress=bool(egr.any()),
                n_pre=int(pre.sum()), n_post=int(post.sum()), n_first_half=int(first_half.sum()),
                n_second_half=int(second_half.sum()))


rows = []
for n in range(1, 13):
    P = DT / n
    kmin = int(np.floor((min(v[0] for v in spans.values()) - TC1) / P)) - 1
    kmax = int(np.ceil((max(v[1] for v in spans.values()) - TC1) / P)) + 1
    for k in range(kmin, kmax + 1):
        tc = TC1 + k * P
        f = k / n
        sig = np.hypot(S1 * (1 - f), S2 * f)
        for s, (lo, hi) in spans.items():
            if lo - 1.0 < tc < hi + 1.0:
                r = dict(n=n, P=P, k=k, tc=tc, sig_tc=sig, sector=s, known=(k == 0 or k == n))
                for p in ["SPOC", "TESS-SPOC", "QLP"]:
                    c = coverage(s, p, tc)
                    if c is not None:
                        r[p] = c
                    c2 = coverage(s, p, tc, strict=False)
                    if c2 is not None:
                        r[p + "_allq"] = c2
                rows.append(r)

json.dump(dict(TC1=TC1, TC2=TC2, DT=DT, T14_d=T14, spans=spans, epochs=rows),
          open(os.path.join(OUT, "epochs.json"), "w"), indent=1)
print("dT", DT, "T14 d", T14)
for s, v in spans.items():
    print("sector", s, v)
for r in rows:
    line = f"n={r['n']:2d} P={r['P']:8.3f} k={r['k']:3d} tc={r['tc']:9.3f}±{r['sig_tc']:.3f} S{r['sector']} known={r['known']}"
    for p in ["SPOC", "TESS-SPOC", "QLP"]:
        if p in r:
            c = r[p]
            line += f" | {p}: in={c['n_in']} frac={c['frac_in']:.2f} I={int(c['ingress'])} E={int(c['egress'])} pre={c['n_pre']} post={c['n_post']}"
            ca = r[p + "_allq"]
            if ca["n_in"] != c["n_in"]:
                line += f" (allq in={ca['n_in']} frac={ca['frac_in']:.2f})"
    print(line)
