"""
Known-status check for the K2-singles leads (run_full2 events whose best period explains >= 2
repeats). This runs BEFORE any hand vetting, so we never re-discover a catalogued period
(lesson from TOI-5560.01).

Per event:
  - TOI / CTOI rows in the bulk ExoFOP tables (data/raw, 2026-10-11) and the TESS EB catalogue
    (Prsa et al. 2022)
  - the ExoFOP target page JSON (TOIs, CTOIs, community planet parameters): fetched live
  - Gaia DR3: the brightest source within 5.4", RUWE, non_single_star, RV variability, any
    NSS two-body orbit, any vari_eclipsing_binary solution (Mowlavi et al. 2023, with period),
    and the vari_classifier best class
  - VSX (VizieR B/vsx/vsx) within 10": name, type, period
  - SIMBAD main type within 5"
Output: k2singles/data/lead_status.csv, plus raw ExoFOP JSON in k2singles/cache/exofop/.
Usage:  python k2singles/scripts/status_check.py [min_support] [run_dir] [EPIC ...]
        (run_dir defaults to full2; with EPICs, rows are appended to the existing CSV)
"""
import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
from gaia_flags import query, rv_variable  # noqa: E402

OUT = ROOT / "k2singles" / "data" / "lead_status.csv"
XDIR = ROOT / "k2singles" / "cache" / "exofop"


def exofop(tic):
    f = XDIR / f"{tic}.json"
    if not f.exists():
        XDIR.mkdir(parents=True, exist_ok=True)
        for k in range(3):
            try:
                r = requests.get(f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json", timeout=60,
                                 headers={"User-Agent": "duotransit-research/0.1"})
                r.raise_for_status()
                f.write_text(r.text, encoding="utf-8")
                break
            except Exception:  # noqa: BLE001
                time.sleep(5 * (k + 1))
    if not f.exists():
        return dict(exofop="fetch failed")
    j = json.loads(f.read_text(encoding="utf-8"))
    tois = [f"{t.get('toi')}:{t.get('tfopwg_disp') or t.get('disp') or ''}" for t in (j.get("tois") or [])]
    ctois = [str(t.get("ctoi")) for t in (j.get("ctois") or [])]
    pp = [f"{p.get('per')}({p.get('user', {}).get('name', '') if isinstance(p.get('user'), dict) else p.get('user', '')})"
          for p in (j.get("planet_parameters") or []) if p.get("per")]
    return dict(exofop_tois=" ".join(tois), exofop_ctois=" ".join(ctois), exofop_periods=" ".join(pp)[:300])


def gaia(ra, dec):
    g = query(f"""SELECT source_id, phot_g_mean_mag, ruwe, non_single_star, rv_chisq_pvalue,
                 rv_renormalised_gof, rv_amplitude_robust, rv_nb_transits,
                 DISTANCE(POINT({ra},{dec}),POINT(ra,dec))*3600 AS sep
                 FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec),CIRCLE({ra},{dec},0.0015))
                 ORDER BY phot_g_mean_mag""")
    if g is None or not len(g):
        return dict(gaia="none")
    sid = int(g["source_id"].iloc[0])
    s = {k: v for k, v in g.iloc[0].to_dict().items() if k != "source_id"}
    s["gaia_source_id"] = str(sid)
    s["gaia_n_within_5as"] = len(g)
    s["rv_variable"] = rv_variable(s.get("rv_chisq_pvalue"), s.get("rv_renormalised_gof"), s.get("rv_nb_transits"))
    n = query(f"SELECT nss_solution_type, period, eccentricity FROM gaiadr3.nss_two_body_orbit WHERE source_id={sid}")
    if n is not None and len(n):
        s.update({f"nss_{k}": v for k, v in n.iloc[0].to_dict().items()})
    e = query(f"SELECT frequency, model_type, global_ranking FROM gaiadr3.vari_eclipsing_binary WHERE source_id={sid}")
    if e is not None and len(e):
        s["gaia_eb_period"] = 1 / float(e["frequency"].iloc[0])
        s["gaia_eb_model"] = e["model_type"].iloc[0]
    c = query(f"SELECT best_class_name, best_class_score FROM gaiadr3.vari_classifier_result WHERE source_id={sid}")
    if c is not None and len(c):
        s["gaia_vari_class"] = f"{c['best_class_name'].iloc[0]}({float(c['best_class_score'].iloc[0]):.2f})"
    return s


def vsx_simbad(ra, dec):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.simbad import Simbad
    from astroquery.vizier import Vizier
    c = SkyCoord(ra, dec, unit="deg")
    out = {}
    try:
        v = Vizier(columns=["Name", "Type", "Period", "max", "min"]).query_region(c, radius=10 * u.arcsec,
                                                                                     catalog="B/vsx/vsx")
        if v:
            t = v[0]
            out["vsx"] = "; ".join(f"{r['Name']}|{r['Type']}|P={r['Period']}" for r in t)
    except Exception as ex:  # noqa: BLE001
        out["vsx"] = f"query failed {str(ex)[:40]}"
    try:
        s = Simbad()
        s.add_votable_fields("otype")
        t = s.query_region(c, radius=5 * u.arcsec)
        if t is not None and len(t):
            cols = {k.lower(): k for k in t.colnames}
            out["simbad"] = "; ".join(f"{r[cols['main_id']]}|{r[cols['otype']]}" for r in t)
    except Exception as ex:  # noqa: BLE001
        out["simbad"] = f"query failed {str(ex)[:40]}"
    return out


def main():
    min_sup = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    run = sys.argv[2] if len(sys.argv) > 2 else "full2"
    only = [int(x) for x in sys.argv[3:]]
    s = pd.read_csv(ROOT / "k2singles" / "results" / run / "summary.csv")
    cat = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv")
    cat["event_idx"] = cat.groupby("EPIC").cumcount()
    s = s.merge(cat[["EPIC", "event_idx", "ra", "dec"]], on=["EPIC", "event_idx"])
    s = s[s.max_support >= min_sup]
    if only:
        s = s[s.EPIC.isin(only)]
    toi = pd.read_csv(ROOT / "data" / "raw" / "exofop_toi.csv", comment="#")
    ctoi = pd.read_csv(ROOT / "data" / "raw" / "exofop_ctoi.csv", comment="#")
    eb = pd.read_csv(ROOT / "data" / "raw" / "tess_eb_prsa2022.tsv", sep="\t", comment="#")
    eb_tics = set(pd.to_numeric(eb.iloc[:, 0], errors="coerce").dropna().astype(int))
    rows = []
    for r in s.itertuples():
        tic = int(r.TIC)
        d = dict(EPIC=r.EPIC, event_idx=r.event_idx, TIC=tic, best_periods=r.periods)
        d["bulk_toi"] = " ".join(map(str, toi[toi["TIC ID"] == tic]["TOI"]))
        d["bulk_ctoi"] = " ".join(map(str, ctoi[ctoi["TIC ID"] == tic]["CTOI"]))
        d["tess_eb_prsa"] = tic in eb_tics
        d.update(exofop(tic))
        d.update(gaia(r.ra, r.dec))
        d.update(vsx_simbad(r.ra, r.dec))
        rows.append(d)
        print(tic, {k: d.get(k) for k in ("bulk_toi", "bulk_ctoi", "exofop_tois", "gaia_eb_period", "vsx")}, flush=True)
    new = pd.DataFrame(rows)
    if only and OUT.exists():
        old = pd.read_csv(OUT)
        new = pd.concat([old[~old.EPIC.isin(only)], new], ignore_index=True)
    new.to_csv(OUT, index=False)


if __name__ == "__main__":
    main()
