"""
Overlap and prior-work check for the T16 Planet Hunt single-transit candidates
(Roth et al., arXiv 2604.18579; table from doi:10.7910/DVN/CWEUGW, file singleTransit_candidate_table.tab).

For each of the 411 stars: is it a TOI, a CTOI (and whose), a confirmed-planet host, a TESS eclipsing
binary (Prsa et al. 2022), or already in this project's lists? Does its transit time coincide with other
T16 singles on the same sector, camera and CCD (a sign of a systematic, not a transit)?

Inputs  : gaia_dr4/data_t16/singleTransit_candidate_table.tab
          data/raw/exofop_toi.csv, exofop_ctoi.csv, nasa_pscomppars.csv  (fetched here if absent)
          data/raw/tess_eb_prsa2022.tsv
          gaia_dr4/candidates_frozen.csv
Outputs : gaia_dr4/data_t16/t16_overlap.csv, t16_overlap_summary.txt
Usage   : python gaia_dr4/t16_overlap.py
"""
import datetime
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
T16 = ROOT / "gaia_dr4" / "data_t16"
URLS = {
    "exofop_toi.csv": "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv",
    "exofop_ctoi.csv": "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv",
    "nasa_pscomppars.csv": ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"
                            "query=select+pl_name,hostname,tic_id+from+pscomppars&format=csv"),
}


def load_t16():
    d = pd.read_csv(T16 / "singleTransit_candidate_table.tab", sep="\t", dtype={"GAIAID": "Int64"})
    return d.rename(columns={"GAIAID": "gaia_dr3", "TICID": "tic", "Sector": "sector", "Camera": "cam",
                             "CCD": "ccd", "Duration": "dur_d", "Tc": "tc_btjd", "Depth": "depth_mag"})


def coincident(d, tol_d=0.1):
    """Number of other rows on the same sector/camera/CCD with a transit time within tol_d days."""
    n = pd.Series(0, index=d.index)
    for _, g in d.groupby(["sector", "cam", "ccd"]):
        t = g.tc_btjd.to_numpy()
        n.loc[g.index] = [(abs(t - x) <= tol_d).sum() - 1 for x in t]
    return n


def fetch_raw():
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        p = RAW / name
        if not p.exists():
            r = requests.get(url, timeout=180, headers={"User-Agent": "python-requests"})
            r.raise_for_status()
            p.write_bytes(r.content)
            print(f"fetched {name} ({len(r.content)} bytes)")


def main():
    fetch_raw()
    d = load_t16()
    toi = pd.read_csv(RAW / "exofop_toi.csv")
    ctoi = pd.read_csv(RAW / "exofop_ctoi.csv")
    ps = pd.read_csv(RAW / "nasa_pscomppars.csv")
    eb = pd.read_csv(RAW / "tess_eb_prsa2022.tsv", sep="\t", comment="#", skiprows=[1, 2],
                     dtype=str).dropna(subset=["TIC"])
    eb = eb[eb.TIC.str.strip().str.isdigit()]
    eb_per = dict(zip(eb.TIC.str.strip().astype(int), pd.to_numeric(eb.Per, errors="coerce")))
    ours = pd.read_csv(ROOT / "gaia_dr4" / "candidates_frozen.csv")
    ours = ours[ours.catalogue != "t16"].groupby("tic").catalogue.agg(lambda s: ";".join(sorted(set(s))))

    tg = toi.groupby("TIC ID")
    toi_ids = tg["TOI"].agg(lambda s: ";".join(f"{v:.2f}" for v in s))
    toi_disp = tg["TFOPWG Disposition"].agg(lambda s: ";".join(sorted(set(s.dropna().astype(str)))))
    toi_per = tg["Period (days)"].agg(lambda s: ";".join(f"{v:.4f}" for v in s.dropna()))
    cg = ctoi.groupby("TIC ID")
    ctoi_user = cg["User"].agg(lambda s: ";".join(sorted(set(s.dropna().astype(str)))))
    ctoi_per = cg["Period (days)"].agg(lambda s: ";".join(f"{v:.4f}" for v in s.dropna()))
    ps_tic = set(pd.to_numeric(ps.tic_id.astype(str).str.replace("TIC ", ""), errors="coerce")
                 .dropna().astype(int))

    d["n_rows_this_tic"] = d.groupby("tic").tic.transform("size")
    d["n_coincident"] = coincident(d)
    d["toi"] = d.tic.map(toi_ids)
    d["toi_disposition"] = d.tic.map(toi_disp)
    d["toi_period_d"] = d.tic.map(toi_per)
    d["ctoi_user"] = d.tic.map(ctoi_user)
    d["ctoi_period_d"] = d.tic.map(ctoi_per)
    d["confirmed_host"] = d.tic.isin(ps_tic)
    d["tess_eb_period_d"] = d.tic.map(eb_per)
    d["in_our_lists"] = d.tic.map(ours)

    def reason(r):
        if r.confirmed_host:
            return "host has a confirmed planet (pscomppars)"
        if pd.notna(r.toi):
            return "star has a TOI"
        if pd.notna(r.tess_eb_period_d) or r.tic in eb_per:
            return "known TESS EB (Prsa+2022)"
        return ""
    d["drop_reason"] = d.apply(reason, axis=1)
    d.to_csv(T16 / "t16_overlap.csv", index=False)

    lines = [
        f"T16 single-transit overlap check, run {datetime.date.today().isoformat()}",
        f"rows {len(d)}; distinct TICs {d.tic.nunique()}; sectors {d.sector.min()}-{d.sector.max()}",
        f"depth (mag): median {d.depth_mag.median():.4f}, max {d.depth_mag.max():.4f}; "
        f"deeper than 0.05 mag: {(d.depth_mag > 0.05).sum()}",
        f"duration (h): median {24 * d.dur_d.median():.1f}, max {24 * d.dur_d.max():.1f}",
        f"star has a TOI: {d.toi.notna().sum()} (dispositions {d.toi_disposition.value_counts().to_dict()})",
        f"star has a CTOI: {d.ctoi_user.notna().sum()} (users {d.ctoi_user.value_counts().head(8).to_dict()})",
        f"  of those, CTOI with a period > 0: "
        f"{d.ctoi_period_d.fillna('').map(lambda s: any(float(x) > 0 for x in s.split(';') if x)).sum()}",
        f"confirmed-planet host: {int(d.confirmed_host.sum())}",
        f"TESS EB (Prsa+2022): {d.tic.isin(eb_per.keys()).sum()}",
        f"already in our lists: {d.in_our_lists.notna().sum()} "
        f"({d.in_our_lists.value_counts().to_dict()})",
        f"transit time shared with another T16 single on the same sector/cam/CCD (0.1 d): "
        f"{(d.n_coincident > 0).sum()}",
        f"dropped: {(d.drop_reason != '').sum()} ({d.drop_reason[d.drop_reason != ''].value_counts().to_dict()})",
        f"kept: {(d.drop_reason == '').sum()} rows on {d.tic[d.drop_reason == ''].nunique()} stars",
    ]
    (T16 / "t16_overlap_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
