"""
Pilot results table: merges results/pilot/<TOI>.json (pipeline), gaia_flags.csv (s06), novelty_check.csv (s04)
screen_pilot.csv (s08 shape screen) and vet_notes.csv (hand-written vetting conclusions per TOI, with the files they rest on).
Output: toi_singles_pilot/pilot_results.csv
"""
import json
from pathlib import Path

import pandas as pd

PILOT = Path(__file__).resolve().parents[1]


def main():
    gaia = pd.read_csv(PILOT / "gaia_flags.csv", dtype={"TOI": str}).set_index("TOI")
    nov = pd.read_csv(PILOT / "novelty_check.csv", dtype={"TOI": str}).set_index("TOI")
    vet = pd.read_csv(PILOT / "vet_notes.csv", dtype={"TOI": str}).set_index("TOI")
    scr = pd.read_csv(PILOT / "screen_pilot.csv", dtype={"TOI": str}).set_index("TOI")
    rows = []
    for f in sorted((PILOT / "results" / "pilot").glob("*.json")):
        if f.name.startswith(("vet_", "screen_", "diffs_")):
            continue
        j = json.loads(f.read_text())
        t = j["TOI"]
        pers = j.get("periods") or []
        top = max([p["n_support"] for p in pers], default=0)
        best = [p["P"] for p in pers if p["n_support"] == top]
        reps = j.get("repeats") or []
        g = gaia.loc[t] if t in gaia.index else None
        n = nov.loc[t] if t in nov.index else None
        v = vet.loc[t] if t in vet.index else None
        sc = scr.loc[t] if t in scr.index else None
        rows.append(dict(
            TOI=t, TIC=j["TIC"], Tmag=j.get("Tmag"), depth_ppm=round(j["depth_ppm"]), dur_h=round(j["dur_h"], 2),
            n_sectors=j.get("n_sectors"), own_photometry_sectors=" ".join(map(str, j.get("own_photometry_sectors") or [])),
            siblings_masked=" ".join(j.get("siblings_masked") or []),
            event_verdict=j.get("event_verdict"), event_snr=j.get("event_snr"), re_anchored=j.get("re_anchored"),
            n_repeats=len(reps),
            repeats=" ".join(f"{r['T']:.2f}(s{r['sector']},{r['depth_ppm']}ppm,snr{r['snr']})" for r in reps),
            max_support=top, n_top_periods=len(best),
            top_periods=" ".join(f"{p:.3f}" for p in sorted(best)[:12]) + (" ..." if len(best) > 12 else ""),
            n_periods=len(pers), min_allowed_P=round(j.get("min_allowed_P_fine") or 0, 1),
            frac_logP_excluded=round(j.get("frac_logP_excluded_fine") or 0, 3), n_allowed_ranges=j.get("n_allowed_ranges"),
            screen_passing=None if sc is None else sc.passing, screen_failing=None if sc is None else sc.failing,
            screened_max_support=None if sc is None else sc.max_support,
            screened_periods=None if sc is None else sc.top_periods,
            diff_images="; ".join(f"{d.get('kind')} s{d.get('sector')}: " + (
                f"{d['offset_arcsec']}+/-{d['offset_err_arcsec']}\" snr{d['snr_peak']}" if "offset_arcsec" in d else
                str(d.get("note") or d.get("error", ""))[:40]) for d in (j.get("diff_images") or [])),
            gaia_ruwe=None if g is None else round(g.ruwe, 2), gaia_binarity_flag=None if g is None else g.binarity_flag,
            gaia_nss=None if g is None or pd.isna(g.get("nss_nss_solution_type")) else
            f"{g.nss_nss_solution_type} P={g.nss_period:.1f} e={g.nss_eccentricity:.2f}",
            exofop_checked_at=None if n is None else n.checked_at,
            exofop_community_uploads=None if n is None else n.community_uploads,
            pscomppars=None if n is None else n.pscomppars,
            pipeline_verdict=None if v is None else v.pipeline_verdict,
            checks=None if v is None else v.checks, conclusion=None if v is None else v.conclusion,
            novelty_status=None if v is None else v.novelty_status, plots=None if v is None else v.plots))
    d = pd.DataFrame(rows)
    d.to_csv(PILOT / "pilot_results.csv", index=False)
    print(d[["TOI", "n_sectors", "event_verdict", "n_repeats", "max_support", "n_top_periods", "top_periods",
             "min_allowed_P", "conclusion", "novelty_status"]].to_string())


if __name__ == "__main__":
    main()
