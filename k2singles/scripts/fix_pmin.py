"""Re-apply the K2-campaign lower period limit to existing pilot JSONs with the fixed k2_coverage
(only the campaign containing the epoch). Candidate periods are re-filtered from candidate_periods_all.
Run: python k2singles/scripts/fix_pmin.py"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_pilot as rp  # noqa: E402

for f in sorted(rp.OUT.glob("*.json")):
    j = json.loads(f.read_text())
    if "EPIC" not in j or j.get("status") == "error":
        continue
    tmin, tmax, iv = rp.k2_coverage(j["EPIC"], j["epoch_btjd"])
    j["k2_tmin_btjd"], j["k2_tmax_btjd"], j["k2_intervals_btjd"] = tmin, tmax, iv
    if tmin is not None:
        pmin = max(j["epoch_btjd"] - tmin, tmax - j["epoch_btjd"])
        j["P_min_from_K2"] = pmin
        cp = j.get("candidate_periods_all", j.get("candidate_periods") or [])
        j["candidate_periods_all"] = cp
        j["candidate_periods"] = [x for x in cp if x[1] >= pmin]
        j["n_candidate_periods"] = len(j["candidate_periods"])
    f.write_text(json.dumps(j, indent=1, default=float))
    print(j["TIC"], j["EPIC"], "intervals", len(iv), "Pmin", None if tmin is None else round(j["P_min_from_K2"], 1),
          "n_periods", j.get("n_candidate_periods"))
