"""
Freeze the blind validation answers. Reads ONLY toi_singles_pilot/results/validation/<TOI>.json (our pipeline
output); no ExoFOP upload, period or comment column is read.

Answer per TOI:
  resolved   one period explains the most repeats and that is >= 2 repeats (>= 3 transits)
  aliases    one or more repeats; the list = every robust_periods period (P >= 5 d) not ruled out by a flat
             window, ordered by support; 'best' = those with the top support
  no repeat  no robust repeat; the allowed period set (fine grid) is reported instead
Output: toi_singles_pilot/validation_blind.csv (+ frozen_at timestamp column).
"""
import json
import time
from pathlib import Path

import pandas as pd

PILOT = Path(__file__).resolve().parents[1]


def answer(j):
    if j.get("status") != "ok":
        return "error", [], []
    pers = j.get("periods") or []
    if not j.get("n_repeats"):
        return "no repeat", [], []
    if not pers:
        return "repeats but no consistent period", [], []
    s = max(p["n_support"] for p in pers)
    best = [p["P"] for p in pers if p["n_support"] == s]
    allp = [p["P"] for p in pers]
    if s >= 2 and len(best) == 1:
        return "resolved", best, allp
    return "aliases", best, allp


def main():
    rows = []
    for f in sorted((PILOT / "results" / "validation").glob("*.json")):
        j = json.loads(f.read_text())
        kind, best, allp = answer(j)
        reps = j.get("repeats") or []
        rows.append(dict(
            TOI=j["TOI"], TIC=j["TIC"], status=j.get("status"), n_sectors=j.get("n_sectors"),
            sectors_loaded=j.get("sectors_loaded"), event_verdict=j.get("event_verdict"),
            event_snr=j.get("event_snr"), event_T=j.get("event_T"), n_repeats=j.get("n_repeats"),
            repeats=" ".join(f"{r['T']:.3f}(s{r['sector']},{r['depth_ppm']}ppm,snr{r['snr']})" for r in reps),
            answer_type=kind, max_support=j.get("max_support"), n_best=len(best),
            best_periods=" ".join(f"{p:.4f}" for p in sorted(best)), n_periods=len(allp),
            all_periods=" ".join(f"{p:.4f}" for p in sorted(allp)),
            min_allowed_P=j.get("min_allowed_P_fine"), frac_logP_excluded=j.get("frac_logP_excluded_fine"),
            n_allowed_ranges=j.get("n_allowed_ranges"),
            note="TOI-7044.01: period is in the TOI comment (flagged; not an independent upload)"
            if j["TOI"] == "7044.01" else "",
            frozen_at=time.strftime("%Y-%m-%d %H:%M:%S")))
    d = pd.DataFrame(rows)
    out = PILOT / "validation_blind.csv"
    if out.exists():
        raise SystemExit("validation_blind.csv already frozen; delete it deliberately to re-freeze")
    d.to_csv(out, index=False)
    print(d[["TOI", "n_sectors", "event_verdict", "n_repeats", "answer_type", "n_best", "best_periods",
             "n_periods", "min_allowed_P"]].to_string())


if __name__ == "__main__":
    main()
