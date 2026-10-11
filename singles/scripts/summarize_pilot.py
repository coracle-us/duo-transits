"""
Single-transit project, step 4: classify pilot results with rules fixed BEFORE looking at the
pilot outcomes (written 2026-10-11 after only the first target had finished).

Difference image (per image; validated controls had peak SNR 8-26, offsets 1-3 arcsec):
  inconclusive   peak SNR < 5
  off-target     offset >= 3 sigma AND >= 6 arcsec (2x the validated floor)
  on-target      otherwise
Target-level source verdict: off-target if any conclusive image is off-target and none is
on-target; on-target if all conclusive images are on-target; mixed if both; inconclusive if
no conclusive image.

Period verdict:
  event-not-recovered   the catalogued event is not 'supported' in our light curves
  period-unique         repeats found and exactly one candidate period survives
  period-narrowed       repeats found, 2-10 candidate periods
  repeats-many-aliases  repeats found, > 10 candidate periods
  repeats-inconsistent  repeats found but no single period explains them all
  single                no repeat detected (report excluded fraction of log P in 5-1500 d)

Usage:  python singles/scripts/summarize_pilot.py [results_subdir]
Output: singles/results/<subdir>/classified.csv
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SNR_MIN, SIG_MIN, ARCSEC_MIN = 5.0, 3.0, 6.0


def image_verdict(d):
    if "offset_arcsec" not in d:
        return "inconclusive"
    if d["snr_peak"] < SNR_MIN:
        return "inconclusive"
    if d["offset_sigma"] >= SIG_MIN and d["offset_arcsec"] >= ARCSEC_MIN:
        return "off-target"
    return "on-target"


def dip_period(times, D, pmin=0.2, pmax=60.0):
    """Best period explaining the detected dips (event + repeats): the P on a fine grid that
    puts the most dips within D/3 of one linear ephemeris through the event. Also tests 2P
    with the dips split between phase 0 and ~0.5 (eclipsing-binary signature)."""
    t = np.sort(np.asarray(times))
    if len(t) < 3:
        return None
    t0 = t[0]
    best = (0, None)
    for P in np.exp(np.linspace(np.log(pmin), np.log(pmax), 200000)):
        r = np.abs((t - t0) / P - np.round((t - t0) / P)) * P
        n = int(np.sum(r < D / 3))
        if n > best[0] or (n == best[0] and best[1] is not None and P > best[1] * 1.5):
            best = (n, P)
    n, P = best
    ph2 = ((t - t0) / (2 * P)) % 1
    second = (np.abs(ph2 - 0.5) < 0.1).sum()
    return dict(P=round(float(P), 5), n_aligned=n, n_dips=len(t),
                frac_at_half_phase_of_2P=round(second / len(t), 2))


def main():
    sub = sys.argv[1] if len(sys.argv) > 1 else "pilot"
    out = ROOT / "singles" / "results" / sub
    rows = []
    for f in sorted(out.glob("*.json")):
        j = json.loads(f.read_text())
        imgs = j.get("diff_images", []) or []
        iv = [(d.get("kind"), d.get("sector"), image_verdict(d), d.get("offset_arcsec"),
               d.get("offset_sigma"), d.get("snr_peak")) for d in imgs]
        concl = [v for v in iv if v[2] != "inconclusive"]
        kinds = {v[2] for v in concl}
        if not concl:
            src = "inconclusive"
        elif kinds == {"on-target"}:
            src = "on-target"
        elif kinds == {"off-target"}:
            src = "off-target"
        else:
            src = "mixed"
        nc = j.get("n_candidate_periods", 0) or 0
        if j.get("status") == "error":
            per = "error"
        elif j.get("event_verdict") != "supported":
            per = "event-not-recovered"
        elif j.get("n_repeats", 0) == 0:
            per = "single"
        elif nc == 0:
            per = "repeats-inconsistent"
        elif nc == 1:
            per = "period-unique"
        elif nc <= 10:
            per = "period-narrowed"
        else:
            per = "repeats-many-aliases"
        depth = j.get("depth_ppm")
        dp = None
        if j.get("n_repeats", 0) >= 3 and per != "period-unique":
            dp = dip_period([j.get("event_T", j["epoch_btjd"])] + [r["T"] for r in j["repeats"]],
                            j["dur_h"] / 24)
        rep_depths = [r["depth_ppm"] for r in j.get("repeats", [])]
        rows.append(dict(
            TIC=j["TIC"], CTOI=j.get("CTOI"), Tmag=j.get("Tmag"), depth_ppm=depth, dur_h=j.get("dur_h"),
            event=j.get("event_verdict"), event_snr=j.get("event_snr"),
            event_depth_ppm=j.get("event_depth_ppm"), n_repeats=j.get("n_repeats"),
            repeat_depth_ratio=" ".join(f"{d / depth:.2f}" for d in rep_depths) if depth else "",
            period_verdict=per, n_candidate_periods=nc,
            candidate_periods=" ".join(f"{p:.4f}" for _, p in (j.get("candidate_periods") or [])[:8]),
            frac_logP_excluded=round(j.get("frac_logP_excluded") or 0, 2),
            min_allowed_P=j.get("min_allowed_P"), source_verdict=src,
            dip_period=json.dumps(dp) if dp else "",
            diff_images=" | ".join(f"{k}:s{s}:{v}({o}\" {sg}sig snr{sn})" for k, s, v, o, sg, sn in iv)))
    df = pd.DataFrame(rows)
    df.to_csv(out / "classified.csv", index=False)
    pd.set_option("display.width", 250, "display.max_colwidth", 90)
    print(df[["TIC", "Tmag", "depth_ppm", "event", "n_repeats", "period_verdict", "candidate_periods",
              "frac_logP_excluded", "source_verdict"]].to_string(index=False))
    print("\nperiod verdicts:", df.period_verdict.value_counts().to_dict())
    print("source verdicts:", df.source_verdict.value_counts().to_dict())
    concl = (df.source_verdict != "inconclusive").mean()
    print(f"targets with >=1 conclusive difference image: {concl:.0%}")


if __name__ == "__main__":
    main()
