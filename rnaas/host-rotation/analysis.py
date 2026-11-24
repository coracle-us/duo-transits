"""
Host-rotation retest of Binnenfeld et al. (arXiv 2611.08267, "B26"), following RULES.md.

  .venv/bin/python analysis.py mma      step 2: reproduce B26 on MMA  -> results/mma.json, results/mma.txt
  .venv/bin/python analysis.py santos   step 3: Santos retest (refuses to run unless results/mma.json passed)
                                         -> results/santos.json, results/santos.txt, results/detection_bias.csv
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RES = HERE / "results"
B26 = np.array([2.846, -0.769, -1.211, -0.559])
NBOOT = 2000
RNG = np.random.default_rng(20261124)


def latest(pattern):
    fs = sorted(DATA.glob(pattern))
    if not fs:
        raise SystemExit(f"missing data/{pattern}; run fetch.py")
    return fs[-1]


# ---------------------------------------------------------------- inputs
def kic_positions():
    pos = []
    for f in DATA.glob("*_table*_*.csv"):
        pos.append(pd.read_csv(f, usecols=["KIC", "_RAJ2000", "_DEJ2000"]))
    k = pd.read_csv(latest("koi_cumulative_*.csv"))[["kepid", "ra", "dec"]]
    k.columns = ["KIC", "_RAJ2000", "_DEJ2000"]
    return pd.concat(pos + [k]).dropna().drop_duplicates("KIC").set_index("KIC")


def unitvec(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.vstack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]).T


def gaia_match(radius=1.0):
    """KIC -> nearest Gaia DR3 source (G < 17.5) within `radius` arcsec, Gaia moved back to J2000 with its proper
    motion. Also records how many sources lie within 2 arcsec."""
    cache = DATA / "gaia_kic_match.csv"
    if cache.exists():
        return pd.read_csv(cache).set_index("KIC")
    from scipy.spatial import cKDTree
    g = pd.read_csv(latest("gaia_field_*.csv.gz"))
    dt_yr = 2016.0 - 2000.0
    dec0 = g.dec - g.pmdec.fillna(0.0) * dt_yr / 3.6e6
    ra0 = g.ra - g.pmra.fillna(0.0) * dt_yr / 3.6e6 / np.cos(np.radians(g.dec))
    tree = cKDTree(unitvec(ra0.values, dec0.values))
    kp = kic_positions()
    xyz = unitvec(kp._RAJ2000.values, kp._DEJ2000.values)
    chord = lambda arcsec: 2 * np.sin(np.radians(arcsec / 3600) / 2)  # noqa: E731
    dist, idx = tree.query(xyz, k=1, distance_upper_bound=chord(radius))
    n2 = np.array([len(x) for x in tree.query_ball_point(xyz, chord(2.0))])
    ok = np.isfinite(dist)
    m = g.iloc[idx[ok]].reset_index(drop=True)
    m.insert(0, "KIC", kp.index.values[ok])
    m["sep_arcsec"] = np.degrees(2 * np.arcsin(dist[ok] / 2)) * 3600
    m["n_within_2arcsec"] = n2[ok]
    m.to_csv(cache, index=False)
    print(f"Gaia match: {ok.sum()} of {len(kp)} KIC stars within {radius} arcsec")
    return m.set_index("KIC")


def hosts():
    k = pd.read_csv(latest("koi_cumulative_*.csv"))
    conf = set(k.loc[k.koi_disposition == "CONFIRMED", "kepid"])
    cand = set(k.loc[k.koi_disposition.isin(["CONFIRMED", "CANDIDATE"]), "kepid"])
    return {"C": conf, "C+K": cand}


def mma():
    """All MMA-analysed stars: KIC, Prot (NaN if none). 2013 KOI table + 2014 table1 + table2."""
    m13 = pd.read_csv(latest("mma2013_table1_*.csv"))[["KIC", "Prot"]]
    m13.loc[~(m13.Prot > 0), "Prot"] = np.nan
    t1 = pd.read_csv(latest("mma2014_table1_*.csv"))[["KIC", "Prot"]]
    t2 = pd.read_csv(latest("mma2014_table2_*.csv"))[["KIC"]].assign(Prot=np.nan)
    d = pd.concat([m13, t1, t2]).sort_values("Prot").drop_duplicates("KIC")  # a KIC with a period keeps it
    return d.set_index("KIC")


def santos(drop_flagged=False):
    a = pd.read_csv(latest("santos2019_table3_*.csv")).rename(columns={"Fl1": "cpcb"})[["KIC", "Prot", "cpcb"]]
    b = pd.read_csv(latest("santos2021_table1_*.csv")).rename(columns={"flag1": "cpcb"})[["KIC", "Prot", "cpcb"]]
    c = pd.read_csv(latest("santos2019_table4_*.csv"))[["KIC"]].assign(Prot=np.nan, cpcb=np.nan)
    e = pd.read_csv(latest("santos2021_table2_*.csv"))[["KIC"]].assign(Prot=np.nan, cpcb=np.nan)
    d = pd.concat([a, b, c, e]).sort_values("Prot").drop_duplicates("KIC")
    if drop_flagged:
        d = d[~(d.cpcb == 1)]
    return d.set_index("KIC")[["Prot"]]


# ---------------------------------------------------------------- cuts
def with_gaia(cat, g):
    d = cat.join(g, how="inner")
    d = d[(d.teff_gspphot >= 4500) & (d.teff_gspphot <= 6000)]
    ok = (d.parallax > 0) & d.ag_gspphot.notna() & d.ebpminrp_gspphot.notna() & d.bp_rp.notna()
    d = d[ok].copy()
    d["MG0"] = d.phot_g_mean_mag - d.ag_gspphot + 5 * np.log10(d.parallax) - 10
    d["c0"] = d.bp_rp - d.ebpminrp_gspphot
    d["dist"] = 1000.0 / d.parallax
    return d


def ms_locus(d):
    """Cubic M_G0((BP-RP)0) through 0.05-mag binned medians, three clipping passes at +-1 mag."""
    keep = np.ones(len(d), bool)
    for _ in range(3):
        x, y = d.c0.values[keep], d.MG0.values[keep]
        edges = np.arange(np.floor(x.min() * 20) / 20, x.max() + 0.05, 0.05)
        idx = np.digitize(x, edges)
        bx, by = [], []
        for i in np.unique(idx):
            s = idx == i
            if s.sum() >= 20:
                bx.append(np.median(x[s]))
                by.append(np.median(y[s]))
        coef = np.polyfit(bx, by, 3)
        keep = np.abs(d.MG0.values - np.polyval(coef, d.c0.values)) <= 1.0
    return coef


def apply_cuts(d, ms_coef, dist_col="dist"):
    ms = np.abs(d.MG0 - np.polyval(ms_coef, d.c0)) <= 1.0
    return d[ms & (d[dist_col] <= 1300) & (d.mh_gspphot > -0.4)].copy()


# ---------------------------------------------------------------- statistic
def xvar(teff):
    return (np.asarray(teff) - 5500.0) / 1000.0


def median_cubic(teff, P, iters=60):
    """L1 (median) regression of ln P on (x, x^2, x^3) by iteratively reweighted least squares."""
    x = xvar(teff)
    X = np.vstack([np.ones_like(x), x, x ** 2, x ** 3]).T
    y = np.log(np.asarray(P))
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    for _ in range(iters):
        w = 1.0 / np.maximum(np.abs(y - X @ b), 1e-4)
        sw = np.sqrt(w)
        bn = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
        if np.max(np.abs(bn - b)) < 1e-7:
            b = bn
            break
        b = bn
    return b


def pexp(coef, teff):
    x = xvar(teff)
    return np.exp(coef[0] + coef[1] * x + coef[2] * x ** 2 + coef[3] * x ** 3)


def compare(s, hostset, lo, hi, coef_fixed=None, nboot=NBOOT):
    """s: cut sample with Prot. Fit on non-hosts 1<P<=30 (unless coef_fixed); compare lo<P<=hi."""
    h = s.index.isin(hostset)
    fitset = s[~h & (s.Prot > 1) & (s.Prot <= 30)]
    coef = median_cubic(fitset.teff_gspphot, fitset.Prot) if coef_fixed is None else np.asarray(coef_fixed)
    win = (s.Prot > lo) & (s.Prot <= hi)
    H, N = s[win & h], s[win & ~h]
    rH = H.Prot.values / pexp(coef, H.teff_gspphot.values)
    rN = N.Prot.values / pexp(coef, N.teff_gspphot.values)
    out = dict(n_hosts=int(len(H)), n_nonhosts=int(len(N)), coef=[round(float(c), 4) for c in coef],
               med_r_hosts=float(np.median(rH)) if len(H) else None, med_r_nonhosts=float(np.median(rN)),
               n_fit=int(len(fitset)))
    if len(H) < 5:
        return out
    out["offset_pct"] = 100 * (out["med_r_hosts"] / out["med_r_nonhosts"] - 1)
    ks = stats.ks_2samp(rH, rN)
    out.update(ks_D=float(ks.statistic), ks_p=float(ks.pvalue))
    # bootstrap: resample hosts and non-hosts; refit the cubic on resampled fit-set non-hosts
    ft, fp = fitset.teff_gspphot.values, fitset.Prot.values
    boots, coefs = [], []
    for _ in range(nboot):
        if coef_fixed is None:
            i = RNG.integers(0, len(ft), len(ft))
            c = median_cubic(ft[i], fp[i])
            coefs.append(c)
        else:
            c = coef
        bh = RNG.integers(0, len(H), len(H))
        bn = RNG.integers(0, len(N), len(N))
        rh = H.Prot.values[bh] / pexp(c, H.teff_gspphot.values[bh])
        rn = N.Prot.values[bn] / pexp(c, N.teff_gspphot.values[bn])
        boots.append(100 * (np.median(rh) / np.median(rn) - 1))
    boots = np.array(boots)
    out.update(ci68=[float(np.percentile(boots, 16)), float(np.percentile(boots, 84))],
               ci95=[float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))])
    if coefs:
        out["coef_boot_sd"] = [round(float(v), 4) for v in np.std(coefs, axis=0)]
    return out


# ---------------------------------------------------------------- steps
def build(catalog, g, ms_coef=None):
    d = with_gaia(catalog, g)
    if ms_coef is None:
        ms_coef = ms_locus(d)
    return d, apply_cuts(d, ms_coef), ms_coef


def fmt(name, o):
    s = f"{name}: hosts {o['n_hosts']}, non-hosts {o['n_nonhosts']}"
    if "offset_pct" in o:
        s += (f"; median r {o['med_r_hosts']:.3f} vs {o['med_r_nonhosts']:.3f}; offset {o['offset_pct']:+.1f}% "
              f"(68% {o['ci68'][0]:+.1f} to {o['ci68'][1]:+.1f}; 95% {o['ci95'][0]:+.1f} to {o['ci95'][1]:+.1f}); "
              f"KS D {o['ks_D']:.3f}, p {o['ks_p']:.2g}; coef {o['coef']}")
    return s


def step_mma():
    RES.mkdir(exist_ok=True)
    g = gaia_match()
    H = hosts()
    allm, cut, ms_coef = build(mma(), g)
    per = cut[cut.Prot.notna()]
    lines, out = [], dict(ms_coef=list(ms_coef), n_gaia_matched=int(mma().index.isin(g.index).sum()))
    lines.append(f"MS locus (cubic in (BP-RP)0): {np.round(ms_coef, 4).tolist()}")
    for hname, hs in H.items():
        region = per  # retained region before any period cut
        o = dict(region_stars=int(len(region)), region_hosts=int(region.index.isin(hs).sum()))
        o["p1_20"] = compare(per, hs, 1, 20)
        o["p20_30"] = compare(per, hs, 20, 30)
        o["p1_20_B26coef"] = compare(per, hs, 1, 20, coef_fixed=B26)
        out[hname] = o
        lines.append(f"[{hname}] region before period cut: {o['region_hosts']} hosts / {o['region_stars']} stars "
                     f"(B26: 400 / 12525)")
        lines.append("  " + fmt("1<P<=20 own fit ", o["p1_20"]) + f"; coef sd {o['p1_20'].get('coef_boot_sd')}")
        lines.append("  " + fmt("1<P<=20 B26 coef", o["p1_20_B26coef"]))
        lines.append("  " + fmt("20<P<=30        ", o["p20_30"]))
    # pass criterion (RULES.md): counts within 10% of 188 / 7274 and offset in 8.7-15.9%
    verdict = {}
    for hname in H:
        p = out[hname]["p1_20"]
        ok = (abs(p["n_hosts"] / 188 - 1) <= 0.10 and abs(p["n_nonhosts"] / 7274 - 1) <= 0.10
              and 8.7 <= p.get("offset_pct", -99) <= 15.9)
        verdict[hname] = bool(ok)
    out["pass"] = verdict
    lines.append(f"PASS (counts within 10%, offset in 8.7-15.9%): {verdict}")
    # sensitivity: GSP-Phot distance
    cut2 = apply_cuts(allm, ms_coef, dist_col="distance_gspphot")
    per2 = cut2[cut2.Prot.notna()]
    for hname, hs in H.items():
        o = compare(per2, hs, 1, 20, nboot=500)
        out[hname]["p1_20_gspphot_distance"] = o
        lines.append(f"  [{hname}] sensitivity, GSP-Phot distance: " + fmt("1<P<=20", o))
    (RES / "mma.json").write_text(json.dumps(out, indent=1))
    (RES / "mma.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def detection_bias(analysed, hostset, label):
    """Per 250-K Teff bin: fraction of hosts and of non-hosts with a period, among all analysed stars passing the
    Gaia cuts (cuts 1-4)."""
    rows = []
    edges = np.arange(4500, 6001, 250)
    h = analysed.index.isin(hostset)
    for lo, hi in zip(edges[:-1], edges[1:]):
        b = (analysed.teff_gspphot >= lo) & (analysed.teff_gspphot < hi if hi < 6000 else analysed.teff_gspphot <= hi)
        for grp, m in (("hosts", h), ("non-hosts", ~h)):
            s = analysed[b & m]
            rows.append(dict(catalog=label, teff_lo=int(lo), teff_hi=int(hi), group=grp, n_analysed=int(len(s)),
                             n_with_period=int(s.Prot.notna().sum()),
                             frac=float(s.Prot.notna().mean()) if len(s) else np.nan))
    return rows


def step_santos():
    mj = RES / "mma.json"
    if not mj.exists():
        raise SystemExit("run the MMA step first")
    mres = json.loads(mj.read_text())
    passed = [k for k, v in mres["pass"].items() if v]
    if not passed:
        raise SystemExit("MMA reproduction did not pass (RULES.md); Santos step refused")
    # primary host definition: the passing one whose counts are closest to B26's
    primary = min(passed, key=lambda k: abs(mres[k]["p1_20"]["n_hosts"] - 188) / 188
                  + abs(mres[k]["p1_20"]["n_nonhosts"] - 7274) / 7274)
    ms_coef = np.array(mres["ms_coef"])
    g = gaia_match()
    H = hosts()
    m_all = apply_cuts(with_gaia(mma(), g), ms_coef)
    s_all = apply_cuts(with_gaia(santos(), g), ms_coef)
    s_unfl = apply_cuts(with_gaia(santos(drop_flagged=True), g), ms_coef)
    mper, sper = m_all[m_all.Prot.notna()], s_all[s_all.Prot.notna()]
    out, lines = dict(primary_hosts=primary), [f"primary host definition: {primary} (both run)"]
    bias = []
    for hname in (primary,) + tuple(k for k in H if k != primary):
        hs = H[hname]
        o = {}
        o["all_1_20"] = compare(sper, hs, 1, 20)
        coefS = o["all_1_20"]["coef"]
        o["all_1_20_B26coef"] = compare(sper, hs, 1, 20, coef_fixed=B26)
        o["all_20_30"] = compare(sper, hs, 20, 30)
        o["all_1_20_noCPCB"] = compare(s_unfl[s_unfl.Prot.notna()], hs, 1, 20, nboot=500)
        both = sper.index.intersection(mper.index)
        o["both_1_20_santosP"] = compare(sper.loc[both], hs, 1, 20, coef_fixed=coefS)
        o["both_1_20_mmaP"] = compare(mper.loc[both], hs, 1, 20, coef_fixed=mres[hname]["p1_20"]["coef"])
        pr = (sper.loc[both, "Prot"] / mper.loc[both, "Prot"])
        win = (mper.loc[both, "Prot"] > 1) & (mper.loc[both, "Prot"] <= 20)
        ph, pn = pr[win & pr.index.isin(hs)], pr[win & ~pr.index.isin(hs)]
        o["both_period_ratio"] = dict(n_hosts=int(len(ph)), n_nonhosts=int(len(pn)),
                                      median_hosts=float(ph.median()), median_nonhosts=float(pn.median()),
                                      frac_within_10pct_hosts=float((abs(ph - 1) < 0.1).mean()),
                                      frac_within_10pct_nonhosts=float((abs(pn - 1) < 0.1).mean()),
                                      ks_p=float(stats.ks_2samp(ph, pn).pvalue))
        only = sper.index.difference(mper.index)
        o["santos_only_1_20"] = compare(sper.loc[only], hs, 1, 20, coef_fixed=coefS)
        out[hname] = o
        lines.append(f"[{hname}]")
        for k in ("all_1_20", "all_1_20_B26coef", "all_20_30", "all_1_20_noCPCB", "both_1_20_santosP",
                  "both_1_20_mmaP", "santos_only_1_20"):
            lines.append("  " + fmt(f"{k:18s}", o[k]))
        lines.append(f"  both, P_Santos/P_MMA (MMA 1-20 d): {o['both_period_ratio']}")
        bias += detection_bias(apply_cuts(with_gaia(mma(), g), ms_coef), hs, f"MMA[{hname}]")
        bias += detection_bias(s_all, hs, f"Santos[{hname}]")
    b = pd.DataFrame(bias)
    b.to_csv(RES / "detection_bias.csv", index=False)
    tot = b.groupby(["catalog", "group"])[["n_analysed", "n_with_period"]].sum()
    tot["frac"] = tot.n_with_period / tot.n_analysed
    lines.append("detection fraction (all Teff bins):\n" + tot.to_string())
    (RES / "santos.json").write_text(json.dumps(out, indent=1))
    (RES / "santos.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    {"mma": step_mma, "santos": step_santos}[sys.argv[1]]()
