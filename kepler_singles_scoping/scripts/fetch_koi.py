"""Download Kepler KOI tables (cumulative + DR25) from the NASA Exoplanet Archive TAP service."""
import requests, io, pandas as pd
from pathlib import Path
OUT = Path(__file__).resolve().parents[1] / "data" / "raw"
TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
q = {
 "koi_cumulative": "select kepid,kepoi_name,kepler_name,koi_disposition,koi_pdisposition,koi_score,koi_fpflag_nt,koi_fpflag_ss,koi_fpflag_co,koi_fpflag_ec,koi_period,koi_period_err1,koi_time0bk,koi_depth,koi_duration,koi_prad,koi_num_transits,koi_count,koi_max_mult_ev,koi_model_snr,koi_steff,koi_srad,koi_smass,ra,dec,koi_kepmag,koi_comment from cumulative",
 "koi_dr25": "select kepid,kepoi_name,koi_disposition,koi_pdisposition,koi_score,koi_period,koi_time0bk,koi_depth,koi_duration,koi_prad,koi_num_transits,koi_model_snr,koi_srad,ra,dec,koi_kepmag from q1_q17_dr25_koi",
}
for name, sql in q.items():
    r = requests.get(TAP, params={"query": sql, "format": "csv"}, timeout=300)
    r.raise_for_status()
    (OUT / f"{name}.csv").write_text(r.text, encoding="utf-8")
    d = pd.read_csv(io.StringIO(r.text))
    print(name, len(d))
