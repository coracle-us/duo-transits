"""
Fable review 3: recompute every number in k2singles/rnaas/rnaas.tex that can be recomputed cheaply.
Reads only; writes check_numbers_out.txt in this folder.
Usage: .venv/Scripts/python k2singles/review_fable3/check_numbers.py
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

ROOT = Path(__file__).resolve().parents[2]
K2 = ROOT / "k2singles"
EPH = K2 / "results" / "ephemerides"
OUT = []


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)


# ---------------------------------------------------------------- 1. blind vs fit agreement
blind1 = {  # independent_table/REPORT.md (P, sigma)
    248811085: (49.25012, 0.00002), 246696483: (182.2452, 0.0003), 247967714: (185.6267, 0.00003),
    246089124: (73.96241, 0.00002), 211892898: (72.5100, 0.0014), 211503363: (122.9049, 0.0007),
    211390677: (184.3103, 0.0010), 211411112: (62.59479, 0.00005), 211489484: (82.25843, 0.00007),
    202071902: (22.904130, 0.000002), 202126877: (21.77521, None), 202072917: (34.745606, 0.000007),
}
blind2 = {246194159: (73.245174, 0.000005), 212012030: (220.32178, 0.00003)}
blind3 = {211924561: (140.49401, 0.00011)}
blind = {**blind1, **blind2, **blind3}
tab = pd.read_csv(K2 / "rnaas" / "table.csv")
p("== 1. blind vs our fit, in sigma (quadrature)")
for _, r in tab.iterrows():
    e = int(r.EPIC)
    j = json.loads((EPH / f"{e}.json").read_text())
    P, sP = r.P_fit, r.sigma_P
    bP, bs = blind[e]
    if bs is None:
        p(f"  {e}: ours {P:.6f} (chi2/dof {j['chi2']:.1f}/{j['dof']}), blind {bP} (no error; ETVs) diff {P-bP:+.5f} d")
        continue
    nsig = (P - bP) / np.hypot(sP, bs)
    flag = "  <-- >2 sigma" if abs(nsig) > 2 else ""
    p(f"  {e}: ours {P:.7f}+-{sP:.7f} blind {bP}+-{bs} -> {nsig:+.2f} sigma; fit chi2/dof {j['chi2']:.2f}/{j['dof']}{flag}")
    if j["dof"] > 0 and j["chi2"] / j["dof"] > 3:
        p(f"      NOTE chi2/dof = {j['chi2']/j['dof']:.1f}: if sigma_P is unscaled the true error is x{np.sqrt(j['chi2']/j['dof']):.1f}")

# ---------------------------------------------------------------- 2. counts
p("\n== 2. counts")
cat = pd.read_csv(K2 / "data" / "k2_singles.csv")
p(f"  catalogue rows {len(cat)}, unique EPIC {cat.EPIC.nunique()}, campaigns: {sorted(cat.campaign.astype(str).unique())}")
s = cat[(cat.n_tess_observed > 0)]
p(f"  TESS-observed {len(s)}; SNR>=7 {(s.tess_snr_per_transit >= 7).sum()}; +depth<1e6 {((s.tess_snr_per_transit >= 7) & (s.depth_ppm < 1e6)).sum()}")
sm = pd.read_csv(K2 / "results" / "full3" / "summary.csv")
p(f"  full3 rows {len(sm)}; n_repeats>0 {(sm.n_repeats > 0).sum()}; max_support>=2 {(sm.max_support >= 2).sum()}")
p(f"  246194159 in full3: n_repeats {sm[sm.EPIC == 246194159].n_repeats.tolist()} max_support {sm[sm.EPIC == 246194159].max_support.tolist()}")

# ---------------------------------------------------------------- 3. TOI-5560.01 numbers
p("\n== 3. TOI-5560.01")
j = json.loads((EPH / "211924561.json").read_text())
P, sP, T0, sT0 = j["P"], j["sigma_P"], j["T0"], j["sigma_T0"]
p(f"  P {P:.5f}+-{sP:.5f}; P/2 {P/2:.6f}+-{sP/2:.6f}; T0 BJD {T0+2457000:.4f}+-{sT0:.4f}; cycles {j['cycles']} gcd {np.gcd.reduce(j['cycles'][1:])}")
depths = [e["depth"] * 1e3 for e in j["events"]]
durs = [e["T14_h"] for e in j["events"]]
p(f"  trapezoid depths ppt {np.round(depths, 2)}; T14 h {np.round(durs, 2)} (mean {np.mean(durs):.2f})")
star = SkyCoord(132.652546, 19.091365, unit="deg")
for label, n, per in (("shared, cycle 30 of 140.49", 30, P), ("70-d only, cycle 61 of 70.247", 61, P / 2)):
    tc = T0 + n * per
    stc = np.hypot(sT0, n * (sP if per == P else sP / 2))
    t = Time(tc + 2457000, format="jd", scale="tdb")
    # BJD_TDB -> geocentric UTC: subtract the light-travel (Roemer) delay toward the star
    ltt = t.light_travel_time(star, kind="barycentric", location=EarthLocation.of_site("greenwich"))
    t_geo = (t - ltt).utc
    p(f"  {label}: BTJD {tc:.4f}+-{stc:.4f} ({stc*1440:.1f} min); BJD_TDB as if UTC: {t.utc.iso[:16]}; "
      f"geocentric UTC (minus {ltt.to(u.min).value:+.1f} min light time): {t_geo.iso[:16]}")

# M2 sin i limit, reproduce gaia_rv_limit.py
G = 6.674e-11
MSUN, DAY = 1.989e30, 86400.0
from scipy import optimize, stats
RV_ERR, N, AMP = 0.601236, 22, 7.698715
M1 = 1.424


def m2sini(K, Pd, m1=M1):
    f = (K * 1e3) ** 3 * Pd * DAY / (2 * np.pi * G) / MSUN
    return optimize.brentq(lambda m2: m2 ** 3 / (m1 + m2) ** 2 - f, 1e-6, 10)


sig_ep = (RV_ERR - 0.11) * np.sqrt(2 * N / np.pi)
s2max = sig_ep ** 2 * stats.chi2.ppf(0.997, N - 1) / (N - 1)
Kb = np.sqrt(2 * max(s2max - sig_ep ** 2, 0))
Ka = AMP / 2
for Pd in (P, P / 2):
    p(f"  P={Pd:.3f}: K_a={Ka:.2f} -> M2sini<{m2sini(Ka, Pd):.3f}; K_b={Kb:.2f} -> M2sini<{m2sini(Kb, Pd):.3f} Msun")
p(f"  per-epoch scatter {sig_ep:.2f} km/s; expected max-min of 22 Gaussian draws ~ {sig_ep*3.8:.1f} km/s (robust amp is {AMP})")

# ---------------------------------------------------------------- 4. EPIC 211503363
p("\n== 4. EPIC 211503363")
j = json.loads((EPH / "211503363.json").read_text())
RJ_per_RSUN = 695700 / 71492
for d in (0.00337, 0.0037, 0.00378):
    p(f"  depth {d*1e3:.2f} ppt: R2 = {np.sqrt(d)*3.60*RJ_per_RSUN:.2f} RJ (R*=3.60), {np.sqrt(d)*4.098*RJ_per_RSUN:.2f} RJ (R*=4.098); "
      f"Kovacs-style R2 {np.sqrt(0.004019)*4.098:.3f} Rsun")
p(f"  P {j['P']:.5f}+-{j['sigma_P']:.5f}; depths {[round(e['depth']*1e3,2) for e in j['events']]} ppt; T14 {[round(e['T14_h'],1) for e in j['events']]} h")

# ---------------------------------------------------------------- 5. HD 220967 eccentricity
p("\n== 5. HD 220967 e cos w")
phi = 0.86723
dphi = phi - 0.5
p(f"  first-order: e cos w = (pi/2)*(phi-0.5) = {np.pi/2*dphi:.4f}")


def sec_phase(e, w, i=np.pi / 2):
    """Exact phase of secondary conjunction (nu = 3pi/2 - w) relative to primary (nu = pi/2 - w)."""
    def M_of(nu):
        E = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu / 2))
        return E - e * np.sin(E)
    dM = (M_of(3 * np.pi / 2 - w) - M_of(np.pi / 2 - w)) % (2 * np.pi)
    return dM / (2 * np.pi)


# duration ratio gives e sin w: T_sec/T_pri ~ (1 - e sin w)/(1 + e sin w) for near-i=90 (sign depends on convention)
ratio = 13.85 / 7.41
esinw_abs = (ratio - 1) / (ratio + 1)
p(f"  duration ratio {ratio:.3f} -> |e sin w| ~ {esinw_abs:.3f} (first-order, equal-chord approx)")
# solve exactly for e cos w, given e sin w = +-esinw_abs
for sgn in (+1, -1):
    def f(ecosw):
        e = np.hypot(ecosw, sgn * esinw_abs)
        w = np.arctan2(sgn * esinw_abs, ecosw)
        return sec_phase(e, w) - phi
    try:
        ec = optimize.brentq(f, 0.01, 0.95)
        e = np.hypot(ec, esinw_abs)
        p(f"  exact, e sin w = {sgn*esinw_abs:+.3f}: e cos w = {ec:.4f}, e = {e:.3f}")
    except ValueError as ex:
        p(f"  exact, e sin w = {sgn*esinw_abs:+.3f}: no root ({ex})")
# pure e cos w (e sin w = 0) exact
ec0 = optimize.brentq(lambda x: sec_phase(x, 0.0) - phi, 0.01, 0.95)
p(f"  exact with e sin w = 0: e cos w = e = {ec0:.4f}  (lower bound on e)")

# HD 220967 joint fit from the 5 events in the ephemeris JSON (common P, separate T0 for primary/secondary)
j = json.loads((EPH / "246194159.json").read_text())
ev = j["events"]
prim = [(e["T"], e["sigma_T"]) for e in ev if e["kind"] == "primary"]
sec = [(e["T"], e["sigma_T"]) for e in ev if e["kind"] == "secondary"]
P0 = 73.2452
A, b, w = [], [], []
for T, sT in prim:
    n = round((T - prim[0][0]) / P0)
    A.append([1, 0, n]); b.append(T); w.append(1 / sT)
for T, sT in sec:
    n = round((T - sec[0][0]) / P0)
    A.append([0, 1, n]); b.append(T); w.append(1 / sT)
A, b, w = np.array(A, float), np.array(b), np.array(w)
Aw, bw = A * w[:, None], b * w
x, *_ = np.linalg.lstsq(Aw, bw, rcond=None)
cov = np.linalg.inv(Aw.T @ Aw)
res = (A @ x - b) * w
p(f"  joint fit of the 5 JSON eclipse times: P = {x[2]:.7f} +- {np.sqrt(cov[2,2]):.7f}; chi2 {np.sum(res**2):.2f} for {len(b)-3} dof; "
  f"secondary phase {((x[1]-x[0])/x[2]) % 1:.5f}")
p(f"  table.csv / build_table.py hard-code P = 73.245178 +- 0.000007 (claimed source: full2_triage.md, which actually says 73.242)")
p(f"  durations: primary {[round(e['T14_h'],2) for e in ev if e['kind']=='primary']} h, secondary {[round(e['T14_h'],2) for e in ev if e['kind']=='secondary']} h; "
  f"depths {[round(e['depth']*100,1) for e in ev]} %")
for lab, n0, per in (("next primary", None, None),):
    pass
# next eclipses from the JSON ephemeris
for kind, T0k in (("primary", x[0]), ("secondary", x[1])):
    n = int(np.ceil((Time("2026-10-12").jd - 2457000 - T0k) / x[2]))
    for k in range(n, n + 2):
        tc = T0k + k * x[2]
        p(f"    next {kind}: BTJD {tc:.3f} = {Time(tc+2457000, format='jd', scale='tdb').utc.iso[:16]} (BJD_TDB, uncorrected)")

# ---------------------------------------------------------------- 6. EPIC 202126877
p("\n== 6. EPIC 202126877: K2 event vs TESS ephemeris")
j = json.loads((EPH / "202126877.json").read_text())
P, T0 = j["P"], j["T0"]
k2T = -211.5569  # blind fit (independent_table results.csv, EVEREST); catalogue epoch -211.5659
for n in (-124, -123, -125):
    pred = T0 + n * P
    p(f"  cycle {n}: predicted {pred:.4f}, K2 event {k2T}: O-C = {(k2T-pred)*24:+.2f} h  (sigma_pred ~ {np.hypot(j['sigma_T0'], n*j['sigma_P'])*24:.2f} h formal)")
p(f"  JSON ephemeris is TESS-only (cycles {j['cycles']}); chi2/dof {j['chi2']:.0f}/{j['dof']}; O-C min {np.round(j['oc_min'],1)}")
p(f"  blind ephemeris (P 21.77521, T0 2488.5793): cycle -124 -> {2488.5793-124*21.77521:.4f}, O-C {(k2T-(2488.5793-124*21.77521))*24:+.2f} h")
pr = [(e["sector"], e["product"], round(e["depth"] * 100, 1), round(e["T14_h"], 1)) for e in j["events"] if e["kind"] == "primary"]
se = [(e["sector"], e["product"], round(e["depth"] * 100, 1), round(e["T14_h"], 1)) for e in j["events"] if e["kind"] == "secondary"]
p(f"  primaries (sector, product, depth %, T14 h): {pr}")
p(f"  secondaries: {se}")
p(f"  ratio S71-72/S43-45: primary {28.5/20.5:.2f}, secondary {12.7/4.3:.2f}  (pure dilution would give the same ratio)")

# ---------------------------------------------------------------- 7. VSX periods vs sector length
p("\n== 7. VSX periods: eclipses per 27-d sector: 3.842 d -> {:.1f}; 8.697 d -> {:.1f}".format(27 / 3.842, 27 / 8.697))

# ---------------------------------------------------------------- 8. 202071902 vs Kovalev
p("\n== 8. 202071902: ours 22.9041282+-0.0000030 vs Kovalev 22.859+-0.015 -> {:.2f} sigma; vs Guo 22.89 (no error given here)".format((22.9041282 - 22.859) / 0.015))
p("   211892898: ours 72.5089+-0.0014 vs Gaia 72.5336+-0.0288 -> {:.2f} sigma".format((72.5336 - 72.5089) / np.hypot(0.0014, 0.0288)))
p("   211411112: ours 62.59466+-0.00018 vs Sandquist 62.59484+-0.00003 -> {:.2f} sigma".format((62.59484 - 62.59466) / np.hypot(0.00018, 0.00003)))

# ---------------------------------------------------------------- 9. word count (own implementation; does NOT touch the tex)
p("\n== 9. word count")
tex = (K2 / "rnaas" / "rnaas.tex").read_text()


def words(s):
    s = re.sub(r"%.*", "", s)
    s = re.sub(r"\\cite[tp]\{[^}]*\}", "CITE", s)
    s = re.sub(r"\$[^$]*\$", "X", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", "", s)
    s = re.sub(r"[{}~]", " ", s)
    return len([w for w in s.split() if re.search(r"[A-Za-z0-9]", w)])


ab = tex.split("\\begin{abstract}")[1].split("\\end{abstract}")[0]
body = tex.split("\\end{abstract}")[1].split("\\begin{acknowledgments}")[0]
tbl = re.search(r"\\begin\{deluxetable\*\}.*?\\end\{deluxetable\*\}", body, flags=re.S).group(0)
body_nt = body.replace(tbl, "")
body_nt = re.sub(r"\\keywords\{.*?\n", "", body_nt)
ai = tex.split("\\section*{AI disclosure statement}")[1].split("\\begin{thebibliography}")[0]
p(f"  abstract {words(ab)} (limit 150); body excl. table/ack {words(body_nt)} (limit 1350); table text {words(tbl)}; AI statement {words(ai)}")
p(f"  body + AI statement + table caption/notes would be {words(body_nt)+words(ai)+words(tbl)}")

(K2 / "review_fable3" / "check_numbers_out.txt").write_text("\n".join(OUT) + "\n")
