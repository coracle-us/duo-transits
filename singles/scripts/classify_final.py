"""
Final per-star classification of every analysed stale single-transit CTOI (pilot, batch2,
batch3, full), with the strict source-location rules agreed 2026-10-11:

Difference image (per image): conclusive if peak SNR >= 5; 'off' if offset >= 6" and >= 3 sigma.
Source verdict per star:
  off-target (firm)      >= 2 conclusive images, all off, offsets mutually consistent
                         (the nearest Gaia star to the fitted source is the same star, or
                         the offsets agree within 10")
  off-target (possible)  exactly 1 conclusive image, off; or off images disagree
  on-target              all conclusive images on (>= 1)
  mixed                  conclusive on and off images
  inconclusive           no conclusive image
For off-target, the eclipse depth the neighbour would need: depth * (1 + 10^(0.4 (G_nb - G_t))).

Star class (first match):
  error / event-not-recovered
  off-target-EB          repeats and firm off-target (or single event firmly off-target)
  EB-on-target (Gaia)    repeats, not off-target, Gaia binarity flag (non_single_star > 0 or
                         RV-variable per Katz+23 cut; RUWE alone is not enough)
  periodic-on-target     repeats, unique period, not off-target, Gaia-clean  -> vet by hand
  duotransit / narrowed  repeats with 2+ candidate periods
  repeats-inconsistent   repeats but no single period explains them
  single                 no repeat; report excluded fraction of log P (5-1500 d) and min P

Manual verdicts from the vetting notes override where listed (MANUAL below).
Output: singles/results/final_classification.csv, printed summary.
Usage:  python singles/scripts/classify_final.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "singles" / "results"
RUNS = ["pilot", "batch2", "batch3", "full"] + sorted(p.name for p in RES.glob("rerun_*") if p.is_dir())
# rerun_<date> directories (new_sector_rerun.py) come last, so the newest result per star wins
SNR_MIN, SIG_MIN, ARCSEC_MIN = 5.0, 3.0, 6.0

MANUAL = {
    147035057: "periodic 37.31497 d (3 eclipses incl. own S108); Gaia RV-variable -> likely EB; see writeups/tic147035057",
    268534931: "tentative 3rd transit S107 (own photometry, egress in gap at momentum dump) -> recheck list",
    64549889: "periodic 35.753 d on target; host A star R*=3.5 Rsun -> companion ~2.1 RJ, likely stellar; G=15.9 star at 3.6\"",
    53847727: "known TESS EB P=3.0113 d (Prsa+2022), rediscovered",
    231969683: "EB: TESS 21.978 d = Gaia SB1 21.969 d",
    230130556: "EB: TESS 22.259 d = Gaia SB1 22.245 d",
    458419988: "EB: Gaia SB1 18.355 d = dT/40",
    58946313: "on-target duotransit but depths differ at 8 sigma in the same SPOC product -> EB-like",
    67884964: "14.42 d rejected on vetting: no dips at predicted S20/S60 epochs, depth chi2 46/4, QLP 0.34x; "
              "submitter (PHT) notes the S47 dip at TJD~2589 is background -> S47 'repeats' likely systematics",
    # hand vetting of the full-run leads, 2026-10-11 (singles/results/vetting/leads_full_run.md)
    159614322: "LEAD: 2nd transit S14 (TESS-SPOC) of S41 event; depths 0.99/1.00, all products agree, on target, "
               "independent re-check agrees; Rp~0.96 RJ; 5 aliases 719.66/n, n=1,2,4,5,8; RUWE 1.51",
    229055790: "LEAD: 2nd transit S68 of PHT S1 event; depths 1.04/0.97, on target, independent re-check agrees; "
               "Rp~3.8 Re; 38 aliases 1838.95/n >= 26.7 d; CHIRON spectra on ExoFOP (2022)",
    139726266: "tentative: 2nd transit only in own S107 photometry (depth 1.07); Rp~2.8 Re on K dwarf; "
               "diff images inconclusive; G=16.1 star at 6.7\" -> recheck list",
    168702709: "tentative: S107 repeat in own photometry is partial (ingress at a gap); 2-min centroid shift 3.9 sigma -> recheck list",
    69745398: "not confirmed: S96 repeat absent in independent re-check; S28 dip shorter than the 17 h catalogued duration",
    337665591: "rejected: depths differ (chi2 12.5/1), crowded field, S66 diff image 8.4\" off at SNR 11",
    365250897: "rejected: depths differ (chi2 40/1), QLP shows no dip",
    51079186: "rejected: depths inconsistent across epochs and products (chi2 22/4)",
}


# Re-run results that hand vetting rejected (the earlier result is kept for these stars).
# Re-run changes NOT listed here enter the classification unvetted; treat them as provisional.
# The RNAAS note uses the frozen copy results/final_classification_2026-10-11_rnaas.csv.
REJECT_RERUN = {
    "rerun_2026-10-11": {
        152359413: "3 'repeats' in S106 own photometry 2-11 d apart: no single P >= 5 d; scattered-light systematics",
        228463329: "3 'repeats' in S107 own photometry 0.8-6.7 d apart: no single P >= 5 d; systematics",
        98794593: "S108 own-photometry 'repeat' sits on a post-gap ramp (+-3 ppt sector trends); diff image 19\" off",
    },
}


def image_ok(d):
    return d.get("snr_peak") is not None and d.get("offset_arcsec") is not None and d["snr_peak"] >= SNR_MIN


def is_off(d):
    return d["offset_arcsec"] >= ARCSEC_MIN and (d.get("offset_sigma") or 0) >= SIG_MIN


def source_verdict(imgs):
    concl = [d for d in imgs if image_ok(d)]
    if not concl:
        return "inconclusive", None
    off = [d for d in concl if is_off(d)]
    on = [d for d in concl if not is_off(d)]
    if off and on:
        return "mixed", None
    if on:
        return "on-target", None
    if len(off) >= 2:
        g = {(round(d.get("nearest_gaia_G") or 0, 1), round(d.get("nearest_gaia_sep") or 0)) for d in off}
        offs = [d["offset_arcsec"] for d in off]
        consistent = len(g) == 1 or (max(offs) - min(offs) <= 10)
        return ("off-target (firm)" if consistent else "off-target (possible)"), off
    return "off-target (possible)", off


def main():
    gaia = pd.read_csv(ROOT / "singles" / "data" / "gaia_binarity_flags.csv").set_index("tic")
    rows = []
    for run in RUNS:
        for f in sorted((RES / run).glob("*.json")):
            j = json.loads(f.read_text())
            tic = int(j["TIC"])
            if run.startswith("rerun_") and tic in REJECT_RERUN.get(run, {}):
                continue  # re-run result vetted and rejected: keep the earlier result
            imgs = [d for d in (j.get("diff_images") or []) if isinstance(d, dict)]
            sv, off = source_verdict(imgs)
            g = gaia.loc[tic] if tic in gaia.index else None
            gaia_bin = bool(g is not None and ((g.get("non_single_star", 0) or 0) > 0 or bool(g.get("rv_variable"))))
            nrep = j.get("n_repeats") or 0
            ncand = len(j.get("candidate_periods") or [])
            req = ""
            if off and g is not None and pd.notna(g.get("phot_g_mean_mag")):
                d0 = off[0]
                if d0.get("nearest_gaia_G"):
                    dg = d0["nearest_gaia_G"] - g["phot_g_mean_mag"]
                    req = round((j.get("depth_ppm") or 0) * 1e-6 * (1 + 10 ** (0.4 * dg)), 3)
            if j.get("status") == "error":
                cls = "error"
            elif j.get("event_verdict") not in ("supported",):
                cls = "event-not-recovered"
            elif sv == "off-target (firm)":
                cls = "off-target-EB" if nrep else "off-target (single event)"
            elif nrep and gaia_bin:
                cls = "EB-on-target (Gaia)"
            elif nrep and ncand == 1:
                cls = "periodic-on-target" if sv in ("on-target", "inconclusive") else "periodic (source " + sv + ")"
            elif nrep and ncand > 1:
                cls = "duotransit / narrowed"
            elif nrep:
                cls = "repeats-inconsistent"
            else:
                cls = "single"
            rows.append(dict(
                TIC=tic, run=run, CTOI=j.get("CTOI"), Tmag=j.get("Tmag"), depth_ppm=j.get("depth_ppm"),
                dur_h=j.get("dur_h"), n_repeats=nrep, n_candidate_periods=ncand,
                periods=" ".join(f"{p:.4f}" for _, p in (j.get("candidate_periods") or [])[:6]),
                frac_logP_excluded=round(j.get("frac_logP_excluded") or 0, 2), min_allowed_P=j.get("min_allowed_P"),
                source=sv, neighbour_required_depth=req, gaia_binary_flag=gaia_bin,
                gaia_ruwe=round(g["ruwe"], 2) if g is not None and pd.notna(g.get("ruwe")) else None,
                own_photometry_sectors=" ".join(map(str, j.get("own_photometry_sectors") or [])),
                cls=cls, manual_note=MANUAL.get(tic, "")))
    df = pd.DataFrame(rows).drop_duplicates("TIC", keep="last")
    df.to_csv(RES / "final_classification.csv", index=False)
    print(f"stars classified: {len(df)}")
    print(df.cls.value_counts().to_string())
    print("\nsource verdicts:", df.source.value_counts().to_dict())
    single = df[df.cls == "single"]
    if len(single):
        print(f"singles: median excluded log-P fraction {single.frac_logP_excluded.median():.2f}")
    keep = df[df.cls.isin(["periodic-on-target", "duotransit / narrowed"]) & ~df.gaia_binary_flag]
    print("\nGaia-clean, not off-target, with repeats (vet by hand):")
    print(keep[["TIC", "run", "Tmag", "depth_ppm", "n_repeats", "periods", "source", "gaia_ruwe", "manual_note"]].to_string(index=False))


if __name__ == "__main__":
    main()
