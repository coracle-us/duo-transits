"""
Hand classification of the 21 catalogued K2 events that are not 'supported' at their catalogued epoch, made from
the diagnostic plots (plots/<EPIC>_0.png, plots/<EPIC>_campaign.png) and the metrics written by audit21.py,
campaign_scan.py, k2sff_apertures.py and t0_211411112.py. The judgements are in CLASSES below (one per event);
this script attaches the supporting numbers and writes audit21.csv and REPORT.txt.

Classes:
 a  event present but removed / gapped / suppressed in PDCSAP (visible in SAP, flagged cadences or another reduction)
 b  event present at a different time (catalogue epoch error)
 c  event present at the epoch but shallower than 0.5x catalogued depth (window mean) or below SNR 5
 d  no event in any available K2 product
 e  other (explained in the note)

Usage:  python k2singles/revision/audit21/classify21.py   (run audit21.py, campaign_scan.py, k2sff_apertures.py
        and t0_211411112.py first)
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import audit21 as A  # noqa: E402
import t0_211411112 as T0  # noqa: E402

warnings.filterwarnings("ignore")
TABLE1 = {211924561, 211503363, 211411112, 246194159, 211489484, 211892898, 248811085, 246089124, 246696483,
          211390677, 247967714, 212012030, 202126877, 202072917, 202071902}

# EPIC: (class, product used for the measured numbers, note)
CLASSES = {
    248873506: ("e", "K2SFF", "K2SFF-only: flat-bottomed 11% box in K2SFF (0.95x) absent from PDCSAP, mission SAP and "
                "EVEREST (<=0.06x); smallest K2SFF aperture shows 0.09x and depth falls in larger apertures, so the "
                "signal is off-centre in K2SFF pixels -> not on the target as seen by PDC/EVEREST (contaminant or "
                "K2SFF artefact). No other dip in C14."),
    246906371: ("c", "PDC_q", "Real dip in PDCSAP at the epoch (min depth 0.74x, ~18,500 ppm) but shorter/V-shaped vs "
                "the catalogued 7.35 h, so the 0.7D window mean is 0.34x. K2SFF 0.61x; K2SFF depth rises with "
                "aperture (0.43x->1.02x), possible blend."),
    201132684: ("e", "EVEREST", "Real ~0.2% event at the epoch in EVEREST (0.59x, SNR 32) and K2SFF (0.55x), constant "
                "across K2SFF apertures; PDCSAP is not gapped but carries ~1.5% uncorrected C10 roll systematics that "
                "swamp it. Second, shorter dip at -0.5 d (comment: multi-planet) also visible."),
    228786343: ("a", "EVEREST", "PDCSAP NaN for 4 of 6 core cadences without any quality flag (PDC gapped the "
                "event); full event in EVEREST (0.55x window, 0.80x min), K2SFF (0.55x), SAP; aperture-independent."),
    220186865: ("c", "PDC_q", "Event present in PDCSAP at the epoch (0.58x window, 0.82x min) but SNR 4.1 < 5 because "
                "the red-noise factor is 2.3 (stellar variability); EVEREST 0.55x SNR 39, K2SFF 0.63x, K2SC 0.44x. "
                "K2SFF depth rises with aperture (0.26x->1.12x), possible blend."),
    220152847: ("a", "EVEREST", "PDCSAP NaN for 12 of 18 core cadences without quality flags (PDC gapped the event); "
                "full event in EVEREST (0.86x), K2SFF (0.84x); remaining PDC cadences at full depth (min 0.96x)."),
    220208795: ("c", "PDC_q", "Deep V-shaped eclipse at the epoch, shorter than the catalogued 6.87 h: PDCSAP window "
                "mean 0.28x (min 0.51x); its deepest cadence is a thruster-firing cadence that is NaN in PDCSAP. "
                "EVEREST 0.80x, K2SFF 0.56x. K2SFF depth rises with aperture (0.31x->1.5x): depth is "
                "aperture-dependent, possible blend."),
    219240689: ("a", "EVEREST", "36-h flat-bottomed 2.2% transit clear in EVEREST (0.84x), K2SFF (0.80x), K2SC (0.77x) "
                "and SAP (0.95x) and constant across K2SFF apertures; PDCSAP co-trending suppressed it to 0.25x and "
                "gapped cadences (rolling band, coarse point, thruster firing)."),
    212820423: ("a", "K2SFF", "Smooth 0.68% ~21-h dip in EVEREST (min 1.0x), K2SFF (0.80x) and K2SFF FRAW (0.86x), "
                "constant across K2SFF apertures; PDCSAP (0.17x) and K2SC (0.22x, built on PDC) cut it into "
                "offset pieces. The EVEREST window test is spoiled by extreme unflagged outliers in its baseline (EVEREST min depth 1.0x)."),
    212343520: ("a", "EVEREST", "PDCSAP NaN for 7 of 9 core cadences (rolling-band, other or no flag); full eclipse in "
                "EVEREST (0.68x), K2SFF (0.68x), SAP (0.74x), aperture-independent."),
    212775521: ("e", "K2SFF", "K2SFF-only: 1% dip (0.67x) at the bottom of a roll-sawtooth ramp, absent in PDCSAP "
                "(0.02x), EVEREST (-0.01x), K2SC (-0.04x) and K2SFF FRAW; K2SFF depth depends on aperture "
                "(0.17x->1.14x). K2SFF correction artefact or off-target; no other dip in C6."),
    212715204: ("a", "EVEREST", "PDCSAP NaN for 5 of 7 core cadences (at most 1 flagged); eclipse in EVEREST and "
                "K2SFF (window 0.49x, min 0.72x; V-shaped, ~15% deep), SAP; aperture-independent."),
    212026444: ("e", "K2SFF", "K2SFF-only clean 700-ppm box (best aperture 0.99x) not reproduced across K2SFF "
                "apertures (-0.4x to 1.8x); EVEREST 0.13x (excluded), K2SC excluded; PDCSAP/SAP show a depression "
                "of about the catalogued depth that is part of longer variability (PDC 0.50x, SNR 2.2, beta 2.5). "
                "Not a secure event."),
    211490542: ("b", "PDC_all", "Catalogued epoch is the SECONDARY eclipse (~9,000 ppm, 0.11x, in every product). "
                "The eclipse with the catalogued 80,800 ppm depth is the primary at BTJD 183.813 (+7.64 d), "
                "present in PDCSAP, EVEREST, K2SFF. The repeat search used the secondary epoch with the primary depth."),
    211411112: ("a", "EVEREST", "PDCSAP NaN for 5 of 6 core cadences with quality 0 (PDC removed the eclipse core); "
                "full eclipse in SAP, EVEREST, K2SFF (0.73x window, 0.93x min). Table 1 target: see note in REPORT."),
    210843533: ("e", "K2SFF", "K2SFF-only (0.37x, two ~6% cadences), absent in PDCSAP (0.00x), EVEREST (0.01x), "
                "K2SC; K2SFF depth grows with aperture (0.37x->3.5x): off-target source (neighbour or moving "
                "object) in the larger K2SFF apertures. Faint target (Kp 16.9)."),
    203746451: ("e", "K2SFF", "No on-target event: PDCSAP corrupted at the epoch (spurious +5% excursion, 8 of 17 "
                "core cadences NaN); EVEREST 0.04x; K2SFF 0.47x with depth growing with aperture (0.37x->3.1x), "
                "i.e. K2SFF-only and off-target."),
    204086428: ("e", "K2SFF", "K2SFF-only (0.57x, SNR 16) at the bottom of a roll-sawtooth ramp; absent in PDCSAP "
                "(0.11x) and EVEREST (0.07x); K2SFF depth grows with aperture (0.24x->1.6x). Artefact or off-target."),
    201775904: ("a", "EVEREST", "PDCSAP NaN for 6 of 8 core cadences (at most 2 flagged); eclipse in EVEREST (0.73x), "
                "K2SFF (0.73x), aperture-independent."),
    202067195: ("e", "K2SFF", "K2SFF-only (0.65x; FRAW 0.57x), absent in PDCSAP (0.10x) and EVEREST (0.01x); K2SFF "
                "depth grows strongly with aperture (0.17x->3.9x): off-target source in the larger apertures (C0)."),
    202085278: ("c", "EVEREST", "Shallow broad depression at the epoch in EVEREST (0.29x window, 0.43x min, SNR 6.6, "
                "beta 2.6) and PDCSAP (0.21x, 0.38x min); K2SFF 0.59x but its depth grows with aperture "
                "(0.16x->3.3x), so the catalogued depth is aperture-dependent; likely an off-target signal (C0)."),
}


def primary_211490542():
    """Time and depth of the deep eclipse near BTJD 183.8 in PDCSAP (all cadences), EVEREST and K2SFF."""
    out = {}
    prods, _ = A.load_products(211490542, 5, 183.82, 10.296 / 24)
    for name in ("PDC_all", "EVEREST", "K2SFF"):
        t, f, flag = prods[name]
        t, f = t[~flag], f[~flag]
        w = np.abs(t - 183.82) < 0.7
        t, f = t[w], f[w]
        p0 = [183.82, 0.08, 0.4, 0.1, np.median(f), 0.0]
        lo, hi = [183.6, 0, 0.1, 0.005, 0.5, -1], [184.0, 0.5, 0.8, 0.4, 1.5, 1]
        r = least_squares(lambda p: T0.model(p, t) - f, p0, bounds=(lo, hi))
        out[name] = dict(t0=r.x[0], depth=r.x[1], t14_h=r.x[2] * 24)
    return out


def main():
    m = pd.read_csv(HERE / "audit21_metrics.csv")
    ap = pd.read_csv(HERE / "k2sff_apertures.csv")
    tg = A.targets().set_index("EPIC")
    prim = primary_211490542()
    rows = []
    for epic, (cls, prod, note) in CLASSES.items():
        r = tg.loc[epic]
        mm = m[m.EPIC == epic].set_index("product")
        avail = [p for p in A.PRODUCTS if bool(mm.loc[p, "available"])]
        pm, pdc = mm.loc[prod], mm.loc["PDC_q"]
        a = ap[(ap.EPIC == epic) & (ap.column == "FCOR")].set_index("aperture").ratio
        row = dict(EPIC=epic, event_idx=int(r.event_idx), campaign=int(r.campaign), Kp=r.Kp,
                   cat_epoch_btjd=r.epoch_btjd, cat_depth_ppm=r.depth_ppm, cat_dur_h=r.dur_h,
                   catalogue_comment=r.comment if isinstance(r.comment, str) else "",
                   check_verdict=r.verdict, cls=cls, measure_product=prod,
                   measured_depth_ppm=pm.depth * 1e6, measured_ratio=pm.ratio, measured_snr=pm.snr,
                   measured_min3_ratio=pm.min3_ratio, pdc_ratio=pdc.ratio, pdc_snr=pdc.snr,
                   pdc_core_cadences=int(pm.n_cad_core), pdc_core_finite=int(pm.n_pdc_finite_core),
                   pdc_core_flagged=int(pm.n_flagged_core),
                   k2sff_ratio_smallest_circ=a.get("CIRC_APER0", np.nan),
                   k2sff_ratio_largest_circ=a.get("CIRC_APER9", np.nan),
                   epoch_offset_d=np.nan, products_checked="+".join(avail), in_table1=epic in TABLE1,
                   table1_depends_on_event="", note=note)
        if epic == 211490542:
            p = prim["PDC_all"]
            row.update(epoch_offset_d=p["t0"] - r.epoch_btjd, measured_depth_ppm=p["depth"] * 1e6,
                       measured_ratio=p["depth"] * 1e6 / r.depth_ppm, measured_snr=np.nan,
                       measure_product="PDC_all (trapezoid fit at the offset epoch)",
                       note=note + " Primary fit: " + "; ".join(
                           f"{k} T={v['t0']:.4f} depth={v['depth'] * 1e6:.0f} ppm T14={v['t14_h']:.1f} h"
                           for k, v in prim.items()))
        if epic == 211411112:
            row["table1_depends_on_event"] = "yes (C5 eclipse time fitted to PDCSAP ingress only)"
        rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(HERE / "audit21.csv", index=False)

    eph = json.loads((A.ROOT / "k2singles" / "results" / "ephemerides" / "211411112.json").read_text())
    t0 = pd.read_csv(HERE / "t0_211411112.csv")
    c5 = t0[(t0.event == "C5")].set_index("product")
    names = {"a": "event present but removed/gapped/suppressed in PDCSAP", "b": "catalogue epoch error",
             "c": "present at epoch, shallower than 0.5x (window mean) or SNR < 5", "d": "no event in any K2 product",
             "e": "other"}
    L = ["Audit of the 21 catalogued K2 events not 'supported' at the catalogued epoch",
         "(k2singles/revision/catalogued_event_check.csv: 7 excluded, 7 ambiguous, 7 unconstrained)",
         "Products: K2 PDCSAP (masked and all cadences), K2 SAP_FLUX (masked and all), EVEREST v2, K2SFF v1 (FCOR",
         "and FRAW, all 20 apertures), K2SC where MAST has it (3 of 21). All HLSP downloads succeeded.",
         "Scripts: audit21.py, campaign_scan.py, k2sff_apertures.py, t0_211411112.py, classify21.py (this dir).",
         "Numbers: 'x' = fraction of the catalogued depth; 'window' = k2lib.test_window_k2 mean over 0.7 D;",
         "'min' = local baseline minus mean of the 3 lowest cadences inside +-D/2.", "",
         "Counts per class:"]
    for k in "abcde":
        L.append(f"  ({k}) {names[k]}: {int((d.cls == k).sum())}")
    e_sff = d[(d.cls == "e") & d.note.str.startswith(("K2SFF-only", "No on-target"))]
    L.append(f"      of (e): {len(e_sff)} K2SFF-only, aperture-dependent signals absent from PDCSAP and EVEREST "
             f"(likely off-target or K2SFF artefacts); 1 real event swamped by PDCSAP systematics (201132684)")
    L += ["", "Per event (EPIC, campaign, check verdict -> class: note):"]
    for _, x in d.iterrows():
        L.append(f"{x.EPIC} C{x.campaign} {x.check_verdict:13s} -> ({x.cls}) [comment: "
                 f"{x.catalogue_comment or '-'}] {x.note}")
    L += ["", "Table 1 EPICs among the 21: " + ", ".join(str(e) for e in d[d.in_table1].EPIC),
          f"211411112 (Sanders 617): Table 1 P = {eph['P']:.5f} +- {eph['sigma_P']:.5f} d uses the C5 eclipse "
          f"timed in PDCSAP at T = {eph['events'][0]['T']:.5f}, from the ingress only (chi2 = {eph['chi2']:.1f} "
          f"for 1 dof; O-C of that event {eph['oc_min'][0]:+.1f} min).",
          f"  Re-timed on the full eclipse: EVEREST {c5.loc['EVEREST', 't0']:.5f} +- {c5.loc['EVEREST', 'sig_t0']:.5f}, "
          f"K2SFF {c5.loc['K2SFF', 't0']:.5f} +- {c5.loc['K2SFF', 'sig_t0']:.5f} "
          f"({c5.loc['EVEREST', 'dt_vs_table1_min']:+.1f} min vs the Table 1 fit; "
          f"{c5.loc['EVEREST', 'dt_vs_C16S46_pred_min']:+.1f} min vs the C16+S46 prediction).",
          "  C16 control timed with the same code reproduces the pipeline's C16 time to <0.01 min.",
          "  P from C16+S46 alone (n=23): 62.594786 d. P from C5(EVEREST)+C16+S46: 62.594821 +- 0.000027 d, "
          "chi2 0.62 (1 dof),",
          "  vs Sandquist+2021 62.59484 +- 0.00003 d. The Table 1 value is pulled ~0.00016 d low by the "
          "ingress-only C5 time;",
          "  it is within its (chi2-inflated) error of the literature; the corrected value has a ~6x smaller formal error.",
          "No other Table 1 EPIC is among the 21, so no other Table 1 result relies on a catalogued K2 event that is "
          "missing at its catalogued epoch.", "",
          "Notes for the note/revision:",
          "  - 211490542: the repeat search (results/full3/211490542_0.json: no repeats, 69% of log P excluded) "
          "was run on the secondary-eclipse epoch with the primary depth; that result is not meaningful as stated.",
          "  - The 7 K2SFF-only events (incl. 212026444) suggest the catalogue was built on K2SFF; their TESS repeat "
          "searches test a signal that is probably not on the target.",
          "  - Even among real events, K2SFF depth rises with aperture size for 220186865, 220208795, 246906371, "
          "202085278 (blend flag only; not tested further)."]
    (HERE / "REPORT.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
