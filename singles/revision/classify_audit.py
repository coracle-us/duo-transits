"""
Assemble audit_34.csv: one row per star whose catalogued event was not recovered, with the
reason category assigned by hand from audit_measurements.csv, audit_sector_scan.csv,
offset_check.csv, the ExoFOP records and the diagnostic plots (plots/<TIC>_audit.png).

Categories (the first that applies, in the order a, b, f, d, c, e):
 a  no TESS data at the epoch in any product, TESScut included
 b  data exist at the epoch but the pipeline did not use them (product not loaded, PDCSAP NaN,
    or cadences quality-flagged), and the dip is visible in those unused data
 c  the dip is at the catalogued time in the pipeline data but failed the event test
    (depth < 0.5x catalogued, SNR < 5)
 d  the dip is present at a different time: an epoch error in the CTOI
 e  no dip in any reduction at or near the catalogued time
 f  other (CTOI attached to the wrong TIC, blend shared with a neighbour CTOI)

Usage:  python singles/revision/classify_audit.py
Output: singles/revision/audit_34.csv, audit_34_summary.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

REV = Path(__file__).resolve().parent
ROOT = REV.parents[1]

R = {  # TIC: (category, note)
    152283270: ("d", "Epoch uploaded as an integer day (BJD 2458336). The only feature is a 3.3 ppt dip "
                     "(0.33x catalogued) at +9.5 h (BTJD 1336.41), SNR 43 in SPOC 2-min SAP/PDCSAP and TASOC "
                     "2-min, weak in TESS-SPOC (2.2 ppt) and absent in QLP/TGLC/eleanor/own FFI photometry: an "
                     "aperture-dependent feature, not a 10 ppt transit. Would fail the depth test even at the right time."),
    209217764: ("d", "Dip at +12.05 h (BTJD 3129.45), 13.3 ppt, SNR ~100 in SPOC, TESS-SPOC, QLP, TARS and own "
                     "photometry: the UAN epoch is ~0.5 d early (half-day slip, as for 5 other UAN CTOIs)."),
    408921284: ("a", "Epoch BTJD 3124.50 falls in S66, but TESScut and MAST have this star only in S54 and S81. "
                     "No TESS data at the epoch (wrong epoch or wrong TIC); no 36 ppt dip in S54/S81 either."),
    431945694: ("c", "Dip at the epoch in all products: 5.6 ppt (0.40x catalogued 14 ppt), SNR 16 (PDCSAP); "
                     "4.5-7.0 ppt across SAP/PDCSAP/TESS-SPOC/QLP/own. Fails the 0.5x depth cut; the CTOI depth is overstated."),
    6121688: ("d", "Dip of 12-13 ppt at +12 h (BTJD 2694.16) in SPOC SAP, TESS-SPOC SAP, QLP and own photometry, "
                   "all in quality-flagged cadences (PDCSAP NaN there): UAN half-day epoch slip. Two more ~13 ppt "
                   "dips in S51 (2699.58, 2702.26) and one at 2715.79: periodic/EB-like, not a single transit."),
    130718055: ("d", "Epoch BTJD 1897.82 is in S21, which did not observe the star (S23, 46, 50, 91 only). ExoFOP "
                     "carries a 2026-09-24 correction by another user (single transit in S46, BTJD 2572.68, 4.3 ppt); "
                     "we see it at BTJD 2572.71, 3.9 ppt, SNR 55. Original epoch wrong."),
    1461133: ("f", "Truncated TIC ID: duplicate of CTOI 14611330.01 (same uploader, epoch, depth, duration; "
                   "re-uploaded under the correct TIC 2021-03-11). TIC 1461133 was not observed in S33. TIC "
                   "14611330 is in our sample and its event was recovered."),
    169251938: ("c", "Weak dip at the epoch in SPOC 2-min: 577 ppm (0.49x), SNR 4.2 (PDCSAP), 715 ppm SNR 3.8 "
                     "(SAP); 2.6 ppt in TESS-SPOC FFI from a few low points; absent in QLP/TGLC/T16/own."),
    250103631: ("c", "Shallow dip at the epoch, consistent across products: 1.2-1.3 ppt (0.28-0.31x of 4.2 ppt), "
                     "SNR 5-7 (SPOC, TESS-SPOC, own). Fails the depth cut; CTOI depth overstated ~3x."),
    255989063: ("e", "No 2 ppt dip at the epoch in any of 6 products (SPOC 2-min depth -109 ppm, SNR -0.1). "
                     "Data end 11 h after the epoch (end of S13)."),
    3181884: ("f", "Truncated TIC ID: duplicate of CTOI 318188400.01 (same epoch, depth, duration). TIC 3181884 "
                   "was not observed in S33. TIC 318188400 is in our sample (classified single, event recovered)."),
    326397318: ("e", "16-min, 1 ppt event: SPOC 2-min gives 1.2 ppt at SNR 2.0 (PDCSAP) / 2.8 (SAP); 30-min "
                     "FFI products cannot resolve it. Not significant in any reduction."),
    329064084: ("d", "Dip at +11.3 h (BTJD 2336.41), 4.8 ppt SNR 17, in all products: UAN half-day slip. Six "
                     "dips in S38 every ~4.75 d: periodic, not a single transit."),
    354032500: ("f", "Blend with the neighbour TIC 354032501 (7.5\" away, 1.7 mag brighter; also a UAN CTOI). The "
                     "event at the epoch (BTJD 3341.04) is 131 ppt in QLP but only 5.8 ppt in the (flagged) SPOC SAP and 9.0 ppt "
                     "(SNR 15, 0.29x) in own photometry; PDCSAP is NaN there. Aperture-dependent depth of a neighbour's "
                     "event. Both stars show the same ~1.6 d dips in S75."),
    354032501: ("d", "Same event as CTOI 354032500.01 (BTJD 3341.04), 0.50 d after this CTOI's epoch: UAN half-day "
                     "slip. Seen in QLP and own photometry; the 2-min cadences there are flagged. Periodic ~1.6 d dips "
                     "in S75 (shared blend)."),
    355095316: ("c", "1.15-h event present at the epoch: 6.8 ppt (1.4x) SNR 4.1 in PDCSAP, 6.4 ppt SNR 4.2 in SAP, "
                     "3.2-6.3 ppt in QLP/TESS-SPOC (QLP KSPSAP SNR 8.5). Fails SNR >= 5 in the pipeline product."),
    356240333: ("c", "420 ppm, 9.7-h event: 434 ppm (1.03x) at SNR 3.0 in PDCSAP, 0.3-0.5 ppt SNR 2.5-6 across "
                     "SPOC/TESS-SPOC/TASOC. Present but below SNR 5."),
    359819459: ("c", "Epoch uploaded as an integer day (BJD 2458558). At that time TESS-SPOC PDCSAP gives 650 ppm "
                     "(0.93x) at SNR 4.1; the best dip, 720 ppm SNR 8, is at -2.2 h (BTJD 1557.91), within the 1-d epoch "
                     "precision but beyond the pipeline's +-0.03 d search. Marginal: S9 has two other dips of similar SNR."),
    36887038: ("c", "Clear dip at the epoch in all products: 1.9 ppt (0.48x of 4.0 ppt), SNR 15.7 (PDCSAP). "
                    "Just below the 0.5x depth cut; CTOI depth overstated 2x."),
    38005437: ("c", "Dip at the epoch: PDCSAP 1.6 ppt (0.74x), SNR 3.3; SAP 3.3 ppt SNR 7.2, TESS-SPOC and own "
                    "'supported' (SNR 5.5-7.9). PDCSAP-specific suppression. Ten dips every ~2.1 d in S62: periodic."),
    425585592: ("c", "Dip at the epoch in SPOC 2-min PDCSAP (1.55 ppt, 0.89x, SNR 3.9) and TESS-SPOC (1.65 ppt, "
                     "SNR 5.8); absent in QLP, TGLC, TARS, eleanor and own photometry. Marginal."),
    427459423: ("c", "19-h, 1.1 ppt event: 516 ppm (0.47x) at SNR 6.2 in PDCSAP; 0.4 ppt SNR 3 in TESS-SPOC; "
                     "absent in QLP/CDIPS/own. Marginal, fails the depth cut."),
    4689632: ("f", "Truncated TIC ID: duplicate of CTOI 468963200.01 (same epoch, depth, duration). TIC 4689632 "
                   "was not observed in S33 (TESScut: S92 only). TIC 468963200 is in our sample (single, recovered)."),
    63700569: ("c", "343 ppm event present at the epoch in all products: 440-620 ppm (1.3-1.8x), SNR 2.0 (PDCSAP), "
                    "4.0-4.6 (SAP, TESS-SPOC, QLP, own). Below SNR 5."),
    64254425: ("d", "Dip at +11.9 h (BTJD 1493.23), 15-18 ppt, SNR 5.6-5.9 in SPOC PDCSAP, TESS-SPOC and TGLC "
                    "(T = 14.2): UAN half-day slip."),
    65290698: ("e", "Epoch in S3 data that are quality-flagged in every product (PDCSAP NaN); neither flagged nor "
                    "unflagged data (12 product/column combinations) show a 2.2 ppt, 9-h dip. The only feature nearby is a "
                    "monotonic post-gap ramp (-3 ppt recovering over ~1 d from BTJD 1396.65)."),
    66013259: ("b", "Epoch in a scattered-light stretch flagged in SPOC/TESS-SPOC/QLP (PDCSAP NaN). TARS shows 6.9 ppt "
                    "(1.03x) SNR 6.7 and own photometry 4.3 ppt SNR 4.1 at the catalogued time; dips recur every "
                    "~2.91 d in S59 on an ephemeris through the epoch: periodic."),
    66101137: ("e", "No 18.7 ppt dip at the epoch in any product (PDCSAP 3.2 ppt SNR 2.9, own 1.1 ppt). A single "
                    "deep event 12.2 d later (BTJD 2251.23; 37 ppt SPOC, 12 ppt own) is not a simple epoch typo."),
    83524788: ("c", "576 ppm event present at the epoch in SPOC and TESS-SPOC: 720-830 ppm (1.3-1.45x), SNR 4.4-4.9; "
                    "own photometry 430 ppm SNR 2.2. Below SNR 5."),
    87949211: ("c", "Dip at the epoch in all products: 2.8-4.7 ppt (0.5-0.8x of 5.7 ppt), SNR 3.1-3.9. Below SNR 5."),
    88212067: ("d", "Dip at +11.96 h (BTJD 2941.06), 4.85 ppt SNR 22-23 in all products: UAN half-day slip. Five "
                    "dips in S60: periodic."),
    88467632: ("e", "PDCSAP NaN at the epoch; the flagged SPOC SAP, TESS-SPOC SAP, CDIPS and own photometry are flat "
                    "(< 0.5 ppt vs 3.5 ppt catalogued)."),
    92709947: ("e", "Flat in all 19 product/column combinations (PDCSAP 375 ppm, SNR 3.2 vs 1.9 ppt catalogued). "
                    "The epoch matches that of CTOI 237332597.02 (S1, 33 deg away) to 2.6 min."),
    95768667: ("d", "Epoch BTJD 1818.35 is in S19, not observed. A 0.8 ppt, SNR 8 dip sits at BTJD 1918.36 (S22), "
                    "exactly +100.0 d: a one-digit typo in the uploaded epoch."),
}

LABELS = {"a": "no TESS data at the epoch", "b": "data exist but not used by the pipeline (flagged / PDC NaN)",
          "c": "dip at the epoch but failed the event test (depth < 0.5x or SNR < 5)",
          "d": "dip at a different time (epoch error in the CTOI)", "e": "no dip in any reduction",
          "f": "other: wrong (truncated) TIC or blend shared with a neighbour CTOI"}
PRIMARY = {"SPOC": "pdcsap_flux", "TESS-SPOC": "pdcsap_flux", "QLP": "sap_flux", "TGLC": "cal_psf_flux",
           "GSFC-ELEANOR-LITE": "corr_flux", "T16": "tfa1", "CDIPS": "tfa1", "TARS": "flux",
           "TASOC_cbv-lc": "flux_corr", "TESScut-own": "aper_flux"}


def main():
    meta = json.loads((REV / "audit_meta.json").read_text())
    m = pd.read_csv(REV / "audit_measurements.csv")
    off = pd.read_csv(REV / "offset_check.csv").set_index("TIC")
    cl = pd.read_csv(ROOT / "singles" / "results" / "final_classification_2026-10-11_rnaas.csv").set_index("TIC")
    rows = []
    for tic, (cat, note) in R.items():
        d = meta[str(tic)]
        x = m[m.TIC == tic]
        parts = []
        for _, r in x.iterrows():
            if PRIMARY.get(r["product"]) != r["column"] and not (r["product"] == "SPOC" and r["column"] == "sap_flux"):
                continue
            mode = "pipe" if r.pipe_verdict not in ("no-data", "unconstrained") or not np.isfinite(r.get("all_depth", np.nan)) else "all"
            dep, snr = r.get(f"{mode}_depth"), r.get(f"{mode}_snr")
            tag = f"{r['product']}{'_' + str(int(r.exptime)) if r.exptime else ''}:{r['column'].replace('_flux', '')}"
            if np.isfinite(dep):
                parts.append(f"{tag}={dep * 1e6:.0f}ppm/SNR{snr:.1f}{'(flagged)' if mode == 'all' else ''}")
            else:
                parts.append(f"{tag}=no data")
        pipe = x[x.pipeline_used]
        o = off.loc[tic] if tic in off.index else None
        rows.append(dict(
            TIC=tic, CTOI=d["CTOI"], uploader=d["user"], Tmag=round(d["Tmag"], 2), epoch_btjd=round(d["epoch_btjd"], 4),
            depth_ppm=d["depth_ppm"], dur_h=d["dur_h"], epoch_sector=d["epoch_sector"],
            observed_at_epoch=bool(d["observed_epoch_sectors"]),
            pipeline_product=";".join(f"S{k}:{v}" for k, v in d["pipeline_product"].items()) or "none",
            mast_products_at_epoch=";".join(d["products_at_epoch_sector"]),
            pipeline_event_verdict=pipe.pipe_verdict.iloc[0] if len(pipe) else "no-data",
            pipeline_depth_ppm=round(pipe.pipe_depth.iloc[0] * 1e6) if len(pipe) and np.isfinite(pipe.pipe_depth.iloc[0]) else np.nan,
            pipeline_snr=round(pipe.pipe_snr.iloc[0], 1) if len(pipe) and np.isfinite(pipe.pipe_snr.iloc[0]) else np.nan,
            dip_at_plus_half_day=(f"{o.bestd_p05:.0f}ppm/SNR{o.bestsnr_p05:.1f} at {o.bestdt_h_p05:+.1f}h"
                                  if o is not None and np.isfinite(o.get("bestsnr_p05", np.nan)) else ""),
            depth_snr_per_product=" | ".join(parts), cls_in_note=cl.loc[tic, "cls"],
            reason=cat, reason_label=LABELS[cat], notes=note,
            exofop_notes=(d.get("exofop_params") or [{}])[0].get("pnotes", "")))
    df = pd.DataFrame(rows)
    df.to_csv(REV / "audit_34.csv", index=False)
    s = df.groupby("reason").agg(n=("TIC", "size"), stars=("TIC", lambda v: " ".join(map(str, v)))).reset_index()
    s["label"] = s.reason.map(LABELS)
    s.to_csv(REV / "audit_34_summary.csv", index=False)
    print(s.to_string(index=False))
    print("uploader x reason:\n", pd.crosstab(df.uploader, df.reason, margins=True))


if __name__ == "__main__":
    main()
