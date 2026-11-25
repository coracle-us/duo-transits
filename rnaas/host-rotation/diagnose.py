"""
Diagnose the non-host excess in the MMA reproduction of B26 (arXiv 2611.08267). See DIAGNOSIS.md.

Each variant changes ONE thing relative to the RULES.md baseline and is compared with every number B26 gives
for the Kepler sample (checkpoints below). Nothing here touches Santos.

  .venv/bin/python diagnose.py            -> results/diagnose.txt, results/diagnose.json
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analysis as A  # noqa: E402

# B26 checkpoints (section 2, Fig. 4 caption, eq. 2)
B26 = dict(region=(12525, 400), far_metal_rich=(902, 9), near_metal_poor=(1398, 18),
           p1_20=(188, 7274), med_r=(0.909, 0.800), offset=13.7, ks_D=0.156, ks_p=2.4e-4,
           coef=[2.846, -0.769, -1.211, -0.559])


def kepler_targets_gaia(g_field):
    """All Kepler targets (keplerstellar, positions from the 2MASS designation) -> nearest Gaia G<17.5 within 1"."""
    from scipy.spatial import cKDTree
    k = pd.read_csv(A.latest("kic_2mass_*.csv")).dropna()
    pat = re.compile(r"J(\d\d)(\d\d)(\d\d\d\d)([+-])(\d\d)(\d\d)(\d\d\d)")
    rows = []
    for kic, des in zip(k.kepid, k.tm_designation):
        m = pat.search(des)
        if not m:
            continue
        h, mi, s, sg, d, dm, ds = m.groups()
        ra = 15 * (int(h) + int(mi) / 60 + int(s) / 100 / 3600)
        dec = (int(d) + int(dm) / 60 + int(ds) / 10 / 3600) * (1 if sg == "+" else -1)
        rows.append((kic, ra, dec))
    kp = pd.DataFrame(rows, columns=["KIC", "ra", "dec"]).drop_duplicates("KIC")
    g = g_field
    dt_yr = 2016.0 - 2000.0
    dec0 = g.dec - g.pmdec.fillna(0.0) * dt_yr / 3.6e6
    ra0 = g.ra - g.pmra.fillna(0.0) * dt_yr / 3.6e6 / np.cos(np.radians(g.dec))
    tree = cKDTree(A.unitvec(ra0.values, dec0.values))
    # 2MASS designations are truncated (RA to 0.01 s, Dec to 0.1"), so allow 1.5"
    dist, idx = tree.query(A.unitvec(kp.ra.values, kp.dec.values), k=1,
                           distance_upper_bound=2 * np.sin(np.radians(1.5 / 3600) / 2))
    ok = np.isfinite(dist)
    m = g.iloc[idx[ok]].reset_index(drop=True)
    m.insert(0, "KIC", kp.KIC.values[ok])
    return m.set_index("KIC")


def prep(d):
    d = d[(d.teff_gspphot >= 4500) & (d.teff_gspphot <= 6000)]
    ok = (d.parallax > 0) & d.ag_gspphot.notna() & d.ebpminrp_gspphot.notna() & d.bp_rp.notna()
    d = d[ok].copy()
    d["MG0"] = d.phot_g_mean_mag - d.ag_gspphot + 5 * np.log10(d.parallax) - 10
    d["c0"] = d.bp_rp - d.ebpminrp_gspphot
    d["dist"] = 1000.0 / d.parallax
    return d


def evaluate(d, ms_coef, hostset, exclude=frozenset(), width=1.0):
    """d: Teff-cut MMA stars with Gaia (prep). Returns every checkpoint."""
    d = d[~d.index.isin(exclude)]
    d = d[d.Prot.notna()]
    ms = np.abs(d.MG0 - np.polyval(ms_coef, d.c0)) <= width
    e = d[ms]
    h = e.index.isin(hostset)
    near, rich = e.dist <= 1300, e.mh_gspphot > -0.4
    out = dict(region=(int((near & rich).sum()), int((near & rich & h).sum())),
               far_metal_rich=(int((~near & rich).sum()), int((~near & rich & h).sum())),
               near_metal_poor=(int((near & ~rich).sum()), int((near & ~rich & h).sum())))
    s = e[near & rich]
    o = A.compare(s, hostset, 1, 20, nboot=1)
    out["p1_20"] = (o["n_hosts"], o["n_nonhosts"])
    out["med_r"] = (round(o["med_r_hosts"], 3), round(o["med_r_nonhosts"], 3))
    out["offset"] = round(o["offset_pct"], 1)
    out["ks_D"], out["ks_p"] = round(o["ks_D"], 3), float(f"{o['ks_p']:.2g}")
    out["coef"] = o["coef"]
    ob = A.compare(s, hostset, 1, 20, coef_fixed=A.B26, nboot=1)
    out["offset_B26coef"] = round(ob["offset_pct"], 1)
    pb = s.Prot
    out["P_bins_region"] = {k: int(v) for k, v in
                            dict(le1=(pb <= 1).sum(), p1_20=((pb > 1) & (pb <= 20)).sum(),
                                 p20_30=((pb > 20) & (pb <= 30)).sum(), gt30=(pb > 30).sum()).items()}
    out["gate"] = bool(abs(o["n_hosts"] / 188 - 1) <= 0.10 and abs(o["n_nonhosts"] / 7274 - 1) <= 0.10
                       and 8.7 <= o["offset_pct"] <= 15.9)
    return out


def line(name, o):
    def r(k):
        a, b = o[k], B26[k]
        return f"{a[0]}/{a[1]} ({a[0] / b[0] - 1:+.0%}/{a[1] / b[1] - 1:+.0%})"
    return (f"{name:34s} region {r('region')}; far rich {r('far_metal_rich')}; near poor {r('near_metal_poor')}; "
            f"1-20 d {r('p1_20')}; med r {o['med_r'][0]}/{o['med_r'][1]}; offset {o['offset']:+.1f}% "
            f"(B26 coef {o['offset_B26coef']:+.1f}%); KS D {o['ks_D']} p {o['ks_p']}; gate {'PASS' if o['gate'] else 'fail'}")


def main():
    g = A.gaia_match()
    H = A.hosts()
    conf = H["C"]
    koi = pd.read_csv(A.latest("koi_cumulative_*.csv"))
    cand_only = set(koi.loc[koi.koi_disposition == "CANDIDATE", "kepid"]) - conf
    any_koi = set(koi.kepid) - conf
    base = prep(A.mma().join(g, how="inner"))
    ms0 = A.ms_locus(base)
    res, lines = {}, []

    def run(name, d=base, ms=ms0, exclude=frozenset(), width=1.0, hostset=conf):
        o = evaluate(d, ms, hostset, exclude, width)
        res[name] = o
        lines.append(line(name, o))
        print(lines[-1], flush=True)
        return o

    lines.append("B26 checkpoints: " + json.dumps(B26))
    run("0 baseline (RULES.md)")
    # 1-2: host bookkeeping
    run("1 candidates out of non-hosts", exclude=cand_only)
    run("2 all non-confirmed KOIs out", exclude=any_koi)
    # 3-5: Gaia quality
    run("3 RUWE < 1.4", d=base[base.ruwe < 1.4])
    run("4a parallax/error > 5", d=base[base.parallax / base.parallax_error > 5])
    run("4b parallax/error > 10", d=base[base.parallax / base.parallax_error > 10])
    run("5 single Gaia source within 2\"", d=base[base.n_within_2arcsec == 1])
    # 6: MS locus fitted to all Kepler targets in the Teff range (Fig. 1 grey points), not MMA-analysed stars
    gf = pd.read_csv(A.latest("gaia_field_*.csv.gz"))
    kt = prep(kepler_targets_gaia(gf).assign(Prot=np.nan))
    ms_kt = A.ms_locus(kt)
    lines.append(f"  MS locus, MMA-analysed {np.round(ms0, 4).tolist()} ({len(base)} stars); "
                 f"all Kepler targets {np.round(ms_kt, 4).tolist()} ({len(kt)} stars)")
    print(lines[-1])
    run("6 MS locus from all Kepler targets", ms=ms_kt)
    # 7: MMA 2014 table1 only for non-KOI stars (drop 2013-KOI-table periods for non-hosts)
    m13 = set(pd.read_csv(A.latest("mma2013_table1_*.csv")).KIC)
    run("7 2013 KOI-table non-hosts out", exclude=m13 - conf)
    # 8: GSP-Phot logg >= 4.0 (not in B26; evolved-star cut by a different route)
    run("8 logg_gspphot >= 4.0", d=base[base.logg_gspphot >= 4.0])

    # ---- round 2 (added after round 1 showed the excess is non-host-only and weighted to P > 20 d)
    lines.append("round 2")
    t14 = pd.read_csv(A.latest("mma2014_table1_*.csv")).set_index("KIC")[["w", "Ref", "n_Prot"]]
    b2 = base.join(t14)
    for th in (0.35, 0.40, 0.45, 0.50):  # MMA 2014 weight; hosts come from the 2013 KOI table and have no w
        run(f"9 MMA2014 w > {th:.2f}", d=b2[b2.w.isna() | (b2.w > th)])
    run("10 MMA2014 Ref=1 out", d=b2[b2.Ref != 1])
    run("11 MMA2014 n_Prot-flagged out", d=b2[b2.n_Prot.isna()])
    for G in (15.0, 15.5, 16.0):  # brightness limit (not in B26)
        run(f"12 G < {G}", d=base[base.phot_g_mean_mag < G])
    lines.append(f"  duplicate Gaia source_id among Teff-cut MMA stars: {int(base.source_id.duplicated(keep=False).sum())}")

    # Fig. 2 of B26: detected-host fraction per bin, read off the rendered figure (+-0.1 percentage point).
    # Our host count / B26 fraction = the star count B26 would have if its hosts match ours; ratio = our excess.
    e = base[base.Prot.notna()]
    e = e[np.abs(e.MG0 - np.polyval(ms0, e.c0)) <= 1.0]
    h = e.index.isin(conf)
    fig2 = dict(
        distance=[((0, 300), 6.4), ((300, 500), 4.6), ((500, 700), 2.85), ((700, 900), 2.8), ((900, 1100), 2.35),
                  ((1100, 1300), 2.35), ((1300, 1500), 1.6), ((1500, 1800), 0.35)],
        mh=[((-1.0, -0.8), 0.95), ((-0.8, -0.6), 1.4), ((-0.6, -0.4), 1.3), ((-0.4, -0.2), 1.8), ((-0.2, 0.0), 2.9),
            ((0.0, 0.2), 5.0), ((0.2, 0.4), 5.55), ((0.4, 0.8), 2.6)])
    res["fig2"] = {}
    for panel, bins in fig2.items():
        col, sel = ("dist", e.mh_gspphot > -0.4) if panel == "distance" else ("mh_gspphot", e.dist <= 1300)
        lines.append(f"  Fig. 2 {panel} panel: bin, our hosts/stars, our %, B26 %, our stars / B26-implied stars")
        for (lo, hi), f in bins:
            b = sel & (e[col] > lo) & (e[col] <= hi)
            n, nh = int(b.sum()), int((b & h).sum())
            ratio = n / (nh / (f / 100)) if nh else float("nan")
            res["fig2"][f"{panel} {lo}..{hi}"] = dict(hosts=nh, stars=n, ours_pct=100 * nh / n, b26_pct=f, excess=ratio)
            lines.append(f"    {lo}..{hi}: {nh}/{n}, {100 * nh / n:.2f}%, {f:.2f}%, x{ratio:.2f}")
    out = dict(checkpoints=B26, ms_locus_mma=list(ms0), ms_locus_kepler=list(ms_kt), variants=res)
    (A.RES / "diagnose.json").write_text(json.dumps(out, indent=1))
    (A.RES / "diagnose.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
