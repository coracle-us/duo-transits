"""
Step 3 of the blend-ceiling test: tag blends among the FP/FA negatives (rules in RULES.md, section 3),
and apply the Gaia rule to the positives as a diagnostic.

Input:  data/sample.csv (from build_sample.py), data/toi_2026-11-24.csv (for Comments)
Cache:  data/tic_gaia.csv (TIC v8 GAIA DR2 ids), data/gaia_dr3_target.csv, data/gaia_dr3_neighbours.csv
Output: data/blend_comment_matches.csv, data/blends.csv, results/blend_fraction.txt
Usage:  .venv/bin/python -I tag_blends.py
"""
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
R = os.path.join(HERE, "results")
os.makedirs(R, exist_ok=True)

PATTERNS = [r"\bNEBs?\b",
            r"nearby\s+(eclipsing\s+binar|EB\b|star|source|target|neighbou?r)",
            r"\bcentroid\s+(offset|shift|motion)",
            r"\boff[- ]target",
            r"\bon\s+(a\s+|the\s+)?(nearby\s+)?(neighbou?r|companion|nearby)",
            r"\bblend(ed|ing)?\b",
            r"\bcontaminat(ing|ion|ed)\s+(star|source|EB)",
            r"\bsignal\s+(is\s+)?(from|on)\s+(TIC|a\s+nearby|the\s+nearby|neighbou?r)",
            r"\bassociated\s+with\s+TIC"]
NEGATIONS = [r"clear", r"\bno\b", r"\bnot\b", r"ruled out", r"excluded", r"unlikely"]
RADIUS = 42.0  # arcsec, 2 TESS pixels


def wilson(k, n, z=norm.ppf(0.975)):
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return p, c - h, c + h


def comment_rule(s):
    toi = pd.read_csv(os.path.join(D, "toi_2026-11-24.csv"), usecols=["TOI", "Comments"])
    com = dict(zip(toi.TOI, toi.Comments.fillna("")))
    rows, hit = [], {}
    for r in s.itertuples():
        segs = re.split(r"[;.]\s+", com.get(r.toi, ""))
        hit[r.toi] = False
        for seg in segs:
            for p in PATTERNS:
                m = re.search(p, seg, re.I)
                if m:
                    neg = next((n for n in NEGATIONS if re.search(n, seg, re.I)), "")
                    rows.append(dict(toi=r.toi, tic=r.tic, disposition=r.disposition, pattern=p,
                                     match=m.group(0), negated_by=neg, segment=seg))
                    hit[r.toi] |= not neg
    pd.DataFrame(rows).to_csv(os.path.join(D, "blend_comment_matches.csv"), index=False)
    return hit


def tic_gaia(tics):
    f = os.path.join(D, "tic_gaia.csv")
    if os.path.exists(f):
        return pd.read_csv(f, dtype={"gaia_dr2": "Int64"})  # 19-digit ids: never via float
    from astroquery.mast import Catalogs
    out = []
    tics = list(tics)
    for i in range(0, len(tics), 200):
        t = Catalogs.query_criteria(catalog="Tic", ID=[int(x) for x in tics[i:i + 200]]).to_pandas()
        out.append(t[["ID", "GAIA", "ra", "dec", "Tmag"]])
    t = pd.concat(out).rename(columns={"ID": "tic", "GAIA": "gaia_dr2"})
    t["tic"] = t.tic.astype(int)
    t["gaia_dr2"] = pd.array([int(x) if str(x) not in ("", "nan", "None", "--") else None for x in t.gaia_dr2],
                             dtype="Int64")
    t.to_csv(f, index=False)
    return t


def gaia_query(q, tries=6):
    """Gaia archive query with retries, cached on disk; alternates sync and async jobs (the async service failed on 2026-11-24)."""
    import hashlib
    import time
    from astroquery.gaia import Gaia
    cdir = os.path.expanduser("~/Documents/GitHub/coracle/duo-transits-cache/blend-ceiling/gaia")
    os.makedirs(cdir, exist_ok=True)
    cf = os.path.join(cdir, hashlib.sha1(q.encode()).hexdigest() + ".csv")
    if os.path.exists(cf):
        return pd.read_csv(cf, dtype={"source_id": "Int64", "dr2_source_id": "Int64", "dr3_source_id": "Int64"})
    import io
    import requests
    for k in range(tries):
        try:
            if k % 2 == 0:  # ARI Heidelberg TAP mirror of the ESA gaiadr3 tables (ESA degraded on 2026-11-24)
                resp = requests.post("https://gaia.ari.uni-heidelberg.de/tap/sync", timeout=180,
                                     data=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q))
                resp.raise_for_status()
                r = pd.read_csv(io.StringIO(resp.text))
            else:
                r = Gaia.launch_job(q).get_results().to_pandas()
            for c in ("source_id", "dr2_source_id", "dr3_source_id"):
                if c in r:
                    r[c] = r[c].astype("Int64")
            r.to_csv(cf, index=False)
            return r
        except Exception as e:
            print(f"  Gaia query failed ({str(e)[:80]}); retry {k + 1}", flush=True)
            time.sleep(20 * (k + 1))
    raise RuntimeError("Gaia query failed after retries: " + q[:120])


def gaia_targets(tg):
    f = os.path.join(D, "gaia_dr3_target.csv")
    if os.path.exists(f):
        return pd.read_csv(f, dtype={"gaia_dr2": "Int64", "dr3_source_id": "Int64"})
    ids = [int(x) for x in tg.gaia_dr2.dropna()]
    rows = []
    for i in range(0, len(ids), 25):
        q = ("SELECT n.dr2_source_id, n.dr3_source_id, n.angular_distance, g.ra, g.dec, g.phot_g_mean_mag "
             "FROM gaiadr3.dr2_neighbourhood AS n JOIN gaiadr3.gaia_source AS g ON g.source_id = n.dr3_source_id "
             f"WHERE n.dr2_source_id IN ({','.join(map(str, ids[i:i + 25]))})")
        rows.append(gaia_query(q))
    n = pd.concat(rows)
    n = n.sort_values("angular_distance").drop_duplicates("dr2_source_id")
    t = tg.merge(n, left_on="gaia_dr2", right_on="dr2_source_id", how="left")
    t["dr3_source_id"] = t.dr3_source_id.astype("Int64")
    t["match"] = np.where(t.dr3_source_id.notna(), "dr2_neighbourhood", "")
    # fallback: brightest DR3 source within 3" of the TIC position (RULES 3)
    for i in t.index[t.dr3_source_id.isna()]:
        q = ("SELECT source_id, ra, dec, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE 1=CONTAINS("
             f"POINT('ICRS', ra, dec), CIRCLE('ICRS', {t.at[i, 'ra_x'] if 'ra_x' in t else t.at[i, 'ra']}, "
             f"{t.at[i, 'dec_x'] if 'dec_x' in t else t.at[i, 'dec']}, {3 / 3600})) ORDER BY phot_g_mean_mag")
        g = gaia_query(q).dropna(subset=["phot_g_mean_mag"])
        if len(g):
            t.loc[i, "dr3_source_id"] = int(g.source_id.iloc[0])
            t.loc[i, "match"] = "brightest_within_3arcsec"
            t.loc[i, ["ra_y", "dec_y", "phot_g_mean_mag"]] = [g.ra.iloc[0], g.dec.iloc[0], g.phot_g_mean_mag.iloc[0]]
    t = t.rename(columns={"ra_y": "g_ra", "dec_y": "g_dec", "ra_x": "tic_ra", "dec_x": "tic_dec",
                          "phot_g_mean_mag": "g_target"})
    t = t[["tic", "gaia_dr2", "dr3_source_id", "match", "angular_distance", "g_ra", "g_dec", "g_target"]]
    t.to_csv(f, index=False)
    return t


def gaia_neighbours(t):
    f = os.path.join(D, "gaia_dr3_neighbours.csv")
    if os.path.exists(f):
        return pd.read_csv(f, dtype={"source_id": "Int64"})
    from concurrent.futures import ThreadPoolExecutor
    m = t.dropna(subset=["dr3_source_id", "g_ra"])

    def cone(r):
        q = ("SELECT source_id, ra, dec, phot_g_mean_mag, DISTANCE(POINT('ICRS', ra, dec), "
             f"POINT('ICRS', {r.g_ra}, {r.g_dec})) * 3600 AS sep FROM gaiadr3.gaia_source WHERE 1=CONTAINS("
             f"POINT('ICRS', ra, dec), CIRCLE('ICRS', {r.g_ra}, {r.g_dec}, {RADIUS / 3600}))")
        g = gaia_query(q)
        g = g[g.source_id.astype("Int64") != int(r.dr3_source_id)]
        g.insert(0, "tic", r.tic)
        return g

    rows = []
    with ThreadPoolExecutor(4) as ex:
        for k, g in enumerate(ex.map(cone, list(m.itertuples())), 1):
            rows.append(g)
            if k % 50 == 0:
                print(f"  Gaia cones {k}/{len(m)}", flush=True)
    n = pd.concat(rows, ignore_index=True)
    n.to_csv(f, index=False)
    return n


def main():
    s = pd.read_csv(os.path.join(D, "sample.csv"))
    s = s[s.drop_reason.isna()].copy()
    hit = comment_rule(s)
    s["blend_C"] = s.toi.map(hit)

    tg = tic_gaia(s.tic.unique())
    t = gaia_targets(tg)
    n = gaia_neighbours(t)
    n = n.dropna(subset=["phot_g_mean_mag"]).merge(t[["tic", "g_target"]], on="tic")
    n["frac"] = 1 / (1 + 10 ** (-0.4 * (n.g_target - n.phot_g_mean_mag)))  # flux_n / (flux_t + flux_n)
    best = n.groupby("tic").frac.max()
    s["gaia_unmatched"] = ~s.tic.isin(t.dropna(subset=["dr3_source_id", "g_target"]).tic)
    s["max_neighbour_frac"] = s.tic.map(best).fillna(0.0)
    s["blend_G"] = (~s.gaia_unmatched) & (s.max_neighbour_frac > s.depth_ppm * 1e-6)
    s["blend_CG"] = s.blend_C | s.blend_G
    s.to_csv(os.path.join(D, "blends.csv"), index=False)

    neg, pos = s[s.label == 0], s[s.label == 1]
    out = [f"Sample (kept TOIs): {len(pos)} positives (CP/KP), {len(neg)} negatives (FP/FA: "
           f"{(neg.disposition == 'FP').sum()} FP, {(neg.disposition == 'FA').sum()} FA)",
           f"Gaia target unmatched: {int(neg.gaia_unmatched.sum())} negatives, {int(pos.gaia_unmatched.sum())} positives", ""]
    for col, name in (("blend_C", "comment rule (C)"), ("blend_CG", "comment or Gaia (C or G)"),
                      ("blend_G", "Gaia rule alone (G), for reference")):
        k = int(neg[col].sum())
        f, lo, hi = wilson(k, len(neg))
        out.append(f"{name}: {k} of {len(neg)} negatives, f = {f:.3f} (95% Wilson {lo:.3f}-{hi:.3f}); "
                   f"AUC_max = 1 - f/2 = {1 - f / 2:.3f} ({1 - hi / 2:.3f}-{1 - lo / 2:.3f})")
    kp = int(pos.blend_G.sum())
    f, lo, hi = wilson(kp, len(pos))
    out.append(f"\nDiagnostic: Gaia rule (G) flags {kp} of {len(pos)} positives ({f:.3f}, 95% Wilson {lo:.3f}-{hi:.3f})")
    cm = pd.read_csv(os.path.join(D, "blend_comment_matches.csv"))
    cm = cm[cm.toi.isin(neg.toi)]
    out.append(f"Comment matches on negatives: {len(cm)} matches on {cm.toi.nunique()} TOIs; "
               f"{(cm.negated_by.notna() & (cm.negated_by != '')).sum()} removed by a negation")
    txt = "\n".join(out)
    print(txt)
    open(os.path.join(R, "blend_fraction.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
