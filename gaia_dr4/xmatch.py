"""
Cross-match the frozen host list with a Gaia release and compare Gaia orbits with the candidates.

Stages (each cached as a CSV in --out; delete a file to redo that stage):
  ids.csv     TIC -> Gaia DR2 ID (MAST TIC) -> DR3 source_id (gaiadr3.dr2_neighbourhood, nearest);
              if the TIC has no DR2 ID: brightest gaiadr3 source within 5.4 arcsec, id_method=positional
  source.csv  one gaia_source row per host
  nss.csv     every nss_two_body_orbit row per host (a star can have several)
  report.csv  one row per candidate x NSS solution (or one row with no solution): flags and verdicts

The ID stage always resolves to DR3. For a later release, --schema changes the tables queried in the
source and nss stages; whether DR3 source_ids carry over, or need a neighbourhood table, must be
checked against that release's documentation first (see prereg_dr4.yaml, 'open_items').

Usage: python gaia_dr4/xmatch.py --schema gaiadr3 --out gaia_dr4/results_dr3
"""
import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from core import binarity_flags, exact_ids, period_match, phase_match  # noqa: E402

HERE = Path(__file__).resolve().parent
SOURCE_COLS = ("source_id, ra, dec, parallax, phot_g_mean_mag, grvs_mag, ruwe, non_single_star, "
               "radial_velocity, radial_velocity_error, rv_nb_transits, rv_chisq_pvalue, "
               "rv_renormalised_gof, rv_amplitude_robust")
NSS_COLS = ("source_id, nss_solution_type, period, period_error, eccentricity, t_periastron, "
            "semi_amplitude_primary")


def gaia_query(sql, tries=4):
    from astroquery.gaia import Gaia
    for attempt in range(tries):
        try:
            return Gaia.launch_job(sql).get_results().to_pandas()
        except Exception as e:  # noqa: BLE001
            print(f"  retry {attempt + 1}: {str(e)[:100]}", flush=True)
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"Gaia query failed after {tries} tries: {sql[:120]}")


def chunks(seq, n):
    seq = list(seq)
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def tic_lookup(tics, chunk=50):
    """TIC -> ra, dec, Tmag, DR2 ID from the MAST TIC."""
    from astroquery.mast import Catalogs
    rows = []
    for c in chunks(tics, chunk):
        t = Catalogs.query_criteria(catalog="Tic", ID=[int(x) for x in c])
        rows.append(t["ID", "ra", "dec", "Tmag", "GAIA"].to_pandas())
        print(f"  TIC {sum(len(r) for r in rows)}/{len(tics)}", flush=True)
    d = pd.concat(rows, ignore_index=True)
    d = d.rename(columns={"ID": "tic", "GAIA": "dr2_source_id"})
    d["tic"] = d.tic.astype("int64")
    d["dr2_source_id"] = exact_ids(d.dr2_source_id)
    return d


def resolve_ids(hosts):
    t = tic_lookup(hosts.tic.tolist())
    t = hosts[["tic"]].merge(t, on="tic", how="left")
    nb = []
    for c in chunks(t.dr2_source_id.dropna().astype("int64").unique(), 200):
        nb.append(gaia_query(
            "SELECT dr2_source_id, dr3_source_id, angular_distance, magnitude_difference "
            f"FROM gaiadr3.dr2_neighbourhood WHERE dr2_source_id IN ({','.join(map(str, c))})"))
    nb = pd.concat(nb, ignore_index=True)
    n_match = nb.groupby("dr2_source_id").size().rename("n_dr3_matches")
    best = nb.sort_values("angular_distance").drop_duplicates("dr2_source_id")
    best = best.merge(n_match, on="dr2_source_id")
    # Int64 before the left merge: a missing value would otherwise turn the column into float64
    best["dr2_source_id"] = best.dr2_source_id.astype("Int64")
    best["dr3_source_id"] = best.dr3_source_id.astype("Int64")
    out = t.merge(best, on="dr2_source_id", how="left")
    out["id_method"] = out.dr3_source_id.notna().map({True: "dr2_neighbourhood", False: ""})
    for i in out.index[out.dr3_source_id.isna()]:
        r = out.loc[i]
        if pd.isna(r.ra):
            out.loc[i, "id_method"] = "not in TIC"
            continue
        g = gaia_query(
            "SELECT source_id, DISTANCE(POINT({ra},{dec}),POINT(ra,dec))*3600000 AS sep_mas "
            "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec),CIRCLE({ra},{dec},0.0015)) "
            "ORDER BY phot_g_mean_mag".format(ra=r.ra, dec=r.dec))
        if len(g):
            out.loc[i, "dr3_source_id"] = int(g.source_id.iloc[0])
            out.loc[i, "angular_distance"] = float(g.sep_mas.iloc[0])
            out.loc[i, "id_method"] = "positional"
        else:
            out.loc[i, "id_method"] = "no Gaia source within 5.4 arcsec"
    out["source_id"] = out.dr3_source_id.astype("Int64")
    return out.drop(columns="dr3_source_id")


def fetch(schema, table, cols, ids):
    parts = [gaia_query(f"SELECT {cols} FROM {schema}.{table} "
                        f"WHERE source_id IN ({','.join(map(str, c))})") for c in chunks(ids, 200)]
    d = pd.concat(parts, ignore_index=True)
    d["source_id"] = d.source_id.astype("int64")
    return d


def build_report(cands, ids, source, nss, cfg):
    th = cfg["thresholds"]
    base = cands.merge(ids[["tic", "source_id", "id_method", "n_dr3_matches"]], on="tic", how="left")
    base = base.merge(source, on="source_id", how="left")
    rows = []
    for r in base.to_dict("records"):
        flags = binarity_flags(r, ruwe_max=th["ruwe_max"])
        keep = {k: r[k] for k in ("tic", "catalogue", "cand_id", "t1_btjd", "t2_btjd", "period_d",
                                  "source_id", "id_method", "n_dr3_matches", "phot_g_mean_mag",
                                  "grvs_mag", "ruwe", "non_single_star", "rv_nb_transits",
                                  "rv_amplitude_robust")}
        keep.update(flags)
        sols = nss[nss.source_id == r["source_id"]] if pd.notna(r["source_id"]) else nss.iloc[0:0]
        if not len(sols):
            rows.append({**keep, "n_nss": 0})
            continue
        for s in sols.to_dict("records"):
            row = {**keep, "n_nss": len(sols), "nss_solution_type": s["nss_solution_type"],
                   "nss_period": s["period"], "nss_period_error": s["period_error"],
                   "nss_eccentricity": s["eccentricity"], "nss_k1": s["semi_amplitude_primary"]}
            if pd.notna(r["t1_btjd"]) and pd.notna(r["t2_btjd"]):
                m = phase_match(r["t2_btjd"] - r["t1_btjd"], th["dt_err_d"], s["period"],
                                s["period_error"], nsig=th["nsig"], max_window=th["max_window"])
                row.update(phase_verdict=m["verdict"], phase_n=m["n"], phase_frac=m["frac"],
                           phase_window=m["window"])
            if pd.notna(r["period_d"]):
                m = period_match(r["period_d"], s["period"], s["period_error"], nsig=th["nsig"],
                                 rel_floor=th["period_rel_floor"], harmonics=th["harmonics"])
                row.update(period_verdict=m["verdict"], period_harmonic=m["harmonic"],
                           period_ratio=m["ratio"])
            rows.append(row)
    rep = pd.DataFrame(rows)
    rep["source_id"] = rep.source_id.astype("Int64")
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--schema", default="gaiadr3")
    ap.add_argument("--out", default=str(HERE / "results_dr3"))
    ap.add_argument("--prereg", default=str(HERE / "prereg_dr4.yaml"))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(open(a.prereg))
    hosts = pd.read_csv(HERE / "hosts_frozen.csv")
    cands = pd.read_csv(HERE / "candidates_frozen.csv")

    if (out / "ids.csv").exists():
        ids = pd.read_csv(out / "ids.csv", dtype={"source_id": "Int64", "dr2_source_id": "Int64"})
    else:
        ids = resolve_ids(hosts)
        ids.to_csv(out / "ids.csv", index=False)
    sid = ids.source_id.dropna().astype("int64").unique()
    print(f"hosts {len(hosts)}; with a source_id {ids.source_id.notna().sum()}; "
          f"methods {ids.id_method.value_counts().to_dict()}", flush=True)

    if (out / "source.csv").exists():
        source = pd.read_csv(out / "source.csv", dtype={"source_id": "int64"})
    else:
        source = fetch(a.schema, "gaia_source", SOURCE_COLS, sid)
        source.to_csv(out / "source.csv", index=False)
    if (out / "nss.csv").exists():
        nss = pd.read_csv(out / "nss.csv", dtype={"source_id": "int64"})
    else:
        nss = fetch(a.schema, "nss_two_body_orbit", NSS_COLS, sid)
        nss.to_csv(out / "nss.csv", index=False)

    rep = build_report(cands, ids, source, nss, cfg)
    rep.to_csv(out / "report.csv", index=False)
    host = rep.drop_duplicates("tic")
    print(f"schema {a.schema}: sources {len(source)}, NSS rows {len(nss)} on "
          f"{nss.source_id.nunique()} stars")
    print(f"hosts flagged: any {int(host.any_flag.sum())}, RUWE {int(host.ruwe_flag.sum())}, "
          f"NSS {int(host.nss_flag.sum())}, RV-variable {int(host.rv_flag.sum())}")
    for col in ("phase_verdict", "period_verdict"):
        if col in rep:
            print(col, rep[col].value_counts().to_dict())


if __name__ == "__main__":
    main()
