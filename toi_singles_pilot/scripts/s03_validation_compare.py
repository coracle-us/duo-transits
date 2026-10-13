"""
Compare the frozen blind answers (validation_blind.csv) with the community periods.

Reference periods: toi_singles_k2/targets_merged.csv post_update_periods (community ExoFOP uploads dated
after the TOI's last update), with users/notes; TOI-7044.01's reference is parsed from its TOI comment
(flagged). Full upload rows from toi_singles_k2/exofop_uploads.csv are attached for context.

Match: |P_ours - P_ref| / P_ref < 0.3 %.
  exact        our answer is 'resolved' (or a single best alias) and it matches a reference period
  in_aliases   a reference period matches one of our listed periods (rank reported)
  allowed      no repeat found; the reference period is NOT excluded by our flat windows (not recovered)
  missed       the reference period is excluded by our data, or we resolved a different period
For every reference period, the verdicts of our grid windows at its predicted epochs (our event time +
k P_ref) are listed, to show why it was or was not recovered.
Output: toi_singles_pilot/validation_compare.csv
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

PILOT = Path(__file__).resolve().parents[1]
ROOT = PILOT.parent
TOL = 0.003

# Diagnoses written AFTER the comparison (from our own grid, light curves and the upload rows); evidence in
# the per_reference column and in results/validation/ (grids, screen_*.json), plots/vet_2300.03.png.
DIAG = {
    "2300.03": "MISSED. The 2nd transit (S86, ~BTJD 3651.5) sits in cadences with SPOC quality bit 4096 and NaN "
               "PDCSAP (the uploader notes 'SAP_FLUX only'); we use quality==0 only, so it is not in our data. Our "
               "only 'repeat' (S74, 3313.12) is a segment-start ramp (plots/vet_2300.03.png), which produced the "
               "wrong 914.41 d answer; the post-hoc shape screen (s08) rejects it, leaving 1252.78 allowed.",
    "5893.01": "MISSED, but the reference looks wrong: the RAVEN P = 7.344 d (Lafarga+2026 upload) predicts "
               "transits at its own epoch BTJD 2802.48 (S55) and at 3536.87/3551.56/3558.90 (S82) that are flat in "
               "SPOC, TESS-SPOC and QLP (|depth| < 460 ppm, < 2 sigma, vs 3543 ppm uploaded / ~7000 ppm observed). "
               "TESS data exclude 7.34 d for this event.",
    "2321.01": "ALLOWED, NOT RECOVERED. A TOI duo (S23 + S50). The S50 transit is only SNR 5.5-6.9 in our box test, "
               "below the SNR-8 repeat threshold, so no alias list was produced. RAVEN 12.567 d = 716.35/57 is "
               "consistent with both transits and is not excluded by our data.",
    "7044.01": "Reference is the TOI comment itself (s25 transits), flagged; not an independent upload.",
    "1894.01": "True period is our top-ranked of 6 (support tie). Third transit (S97) is at a segment end "
               "(no right-hand baseline) and the S3 epoch is in a gap of our SPOC_120 data; S63 dip is spurious "
               "(gap-edge). After the post-hoc screen: {369.62, 739.24}.",
}


def epoch_verdicts(g, T0, P, D, tmin, tmax):
    Tg = g.T_pred.to_numpy()
    out = []
    for k in range(int(np.ceil((tmin - T0) / P)), int(np.floor((tmax - T0) / P)) + 1):
        if k == 0:
            continue
        T = T0 + k * P
        i = np.searchsorted(Tg, T)
        best = None
        for ii in (i - 1, i):
            if 0 <= ii < len(Tg) and abs(Tg[ii] - T) <= D / 8 + 1e-9:
                best = g.iloc[ii]
        if best is None:
            continue  # no data
        dep = best.get("depth")
        out.append(f"{T:.2f}:s{int(best.sector)}:{best.verdict}"
                   + (f"({dep * 1e6:.0f}ppm,snr{best.snr:.1f})" if pd.notna(dep) else ""))
    return out


def main():
    b = pd.read_csv(PILOT / "validation_blind.csv", dtype={"TOI": str})
    m = pd.read_csv(ROOT / "toi_singles_k2" / "targets_merged.csv", dtype={"TOI": str})
    m["TOI"] = m.TOI.map(lambda x: f"{float(x):.2f}")
    m = m.set_index("TOI")
    up = pd.read_csv(ROOT / "toi_singles_k2" / "exofop_uploads.csv", dtype={"TOI": str})
    up["TOI"] = up.TOI.map(lambda x: f"{float(x):.2f}")
    rows = []
    for r in b.itertuples():
        mr = m.loc[r.TOI]
        refs = [float(x) for x in re.findall(r"\d+\.?\d*", str(mr.post_update_periods))] \
            if pd.notna(mr.post_update_periods) else []
        src = f"ExoFOP upload ({mr.post_update_users})" if refs else ""
        if r.TOI == "7044.01" and not refs:
            refs = [float(x) for x in re.findall(r"~\s*(\d+\.\d+)\s*d", str(mr.Comments))]
            src = "TOI comment (FLAG: not an independent upload)"
        j = json.loads((PILOT / "results" / "validation" / f"{r.TOI}.json").read_text())
        D = j["dur_h"] / 24
        g = pd.read_csv(PILOT / "results" / "validation" / f"{j['TIC']}_{r.TOI}_grid.csv.gz")
        ours = [float(x) for x in str(r.all_periods).split()] if pd.notna(r.all_periods) else []
        best = [float(x) for x in str(r.best_periods).split()] if pd.notna(r.best_periods) else []
        ranges = j.get("allowed_ranges") or []
        results, whys = [], []
        for P in refs:
            hit = [i for i, q in enumerate(ours) if abs(q - P) / P < TOL]
            # rank in the support-ordered list of the JSON
            jp = j.get("periods") or []
            rank = next((i + 1 for i, q in enumerate(jp) if abs(q["P"] - P) / P < TOL), None)
            allowed = any(a * (1 - 1e-4) <= P <= bb * (1 + 1e-4) for a, bb in ranges) if P >= 5 else None
            if best and len(best) == 1 and abs(best[0] - P) / P < TOL:
                res = "exact"
            elif hit:
                res = f"in_aliases(rank {rank} of {len(jp)})"
            elif r.answer_type == "no repeat" and allowed:
                res = "allowed_not_recovered"
            else:
                res = "missed"
            results.append(res)
            ev = epoch_verdicts(g, j["event_T"], P, D, g.T_pred.min(), g.T_pred.max())
            whys.append(f"P={P}: {res}; allowed_by_flat_windows={allowed}; predicted epochs in our data: "
                        + (" ".join(ev) if ev else "none"))
        order = ["exact", "in_aliases", "allowed_not_recovered", "missed"]
        overall = min(results, key=lambda x: order.index(x.split("(")[0])) if results else "no reference"
        u = up[up.TOI == r.TOI]
        scr_f = PILOT / "results" / "validation" / f"screen_{r.TOI}.json"
        scr = ""
        if scr_f.exists():
            sp = json.loads(scr_f.read_text()).get("periods_screened") or []
            top = max([p["n_support"] for p in sp], default=0)
            sb = sorted(p["P"] for p in sp if p["n_support"] == top)
            scr = " ".join(f"{p:.3f}" for p in sb[:12]) + (f" (+{len(sb) - 12})" if len(sb) > 12 else "")
            scr = scr or "no passing repeat"
            if refs:
                pr = any(abs(p - P) / P < TOL for p in [q["P"] for q in sp] for P in refs)
                scr += f" | ref in screened list: {pr}"
        rows.append(dict(TOI=r.TOI, TIC=r.TIC, reference_periods=" ".join(map(str, refs)), reference_source=src,
                         reference_notes=str(mr.post_update_notes)[:300] if pd.notna(mr.post_update_notes) else "",
                         our_answer_type=r.answer_type, our_best=r.best_periods, our_n_periods=r.n_periods,
                         our_n_repeats=r.n_repeats, our_repeats=r.repeats, our_min_allowed_P=r.min_allowed_P,
                         sectors=r.sectors_loaded, result=overall, diagnosis=DIAG.get(r.TOI, ""),
                         post_hoc_screened_top=scr, per_reference=" | ".join(whys),
                         all_uploads=" | ".join(f"{x.user} {x.date} P={x.period}" for x in u.itertuples())[:600]))
    d = pd.DataFrame(rows)
    d.to_csv(PILOT / "validation_compare.csv", index=False)
    print(d[["TOI", "reference_periods", "our_answer_type", "our_best", "our_n_periods", "result"]].to_string())
    print(d.result.str.split("(").str[0].value_counts().to_string())


if __name__ == "__main__":
    main()
