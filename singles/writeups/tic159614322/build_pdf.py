"""
Technical write-up of TIC 159614322 / CTOI 159614322.01 for the groups who can verify it.
Run:  python singles/writeups/tic159614322/{analysis,s21_check,ephemeris,secondary,make_figures,build_pdf}.py
Output: singles/writeups/tic159614322/TIC159614322_writeup.pdf
Numbers trace to analysis.json, s21_check.json, ephemeris.json, secondary.json, maps.json,
predictions.csv (this folder), independent/results.json (blind re-analysis) and
singles/results/vetting/159614322_n1.json.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (Image, KeepTogether, Paragraph, SimpleDocTemplate, Spacer, Table,
                                TableStyle)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INK, MUTED, ACCENT = colors.HexColor("#1d2b3a"), colors.HexColor("#5d6b7c"), colors.HexColor("#d9682b")
W = letter[0] - 1.4 * inch

title = ParagraphStyle("t", fontName="Helvetica-Bold", fontSize=14.5, leading=18, textColor=INK)
sub = ParagraphStyle("s", fontName="Helvetica", fontSize=8.8, leading=11.5, textColor=MUTED)
h2 = ParagraphStyle("h", fontName="Helvetica-Bold", fontSize=10.5, leading=13, textColor=ACCENT,
                    spaceBefore=7, spaceAfter=3)
body = ParagraphStyle("b", fontName="Helvetica", fontSize=8.9, leading=12.2, textColor=INK, spaceAfter=4)
cap = ParagraphStyle("c", fontName="Helvetica-Oblique", fontSize=7.8, leading=10, textColor=MUTED, spaceAfter=4)
cell = ParagraphStyle("cell", fontName="Helvetica", fontSize=7.8, leading=9.6, textColor=INK)
cellb = ParagraphStyle("cellb", parent=cell, fontName="Helvetica-Bold")

# Summary of the blind re-analysis (independent/report.md); set after reading it.
INDEPENDENT = ("A separate agent re-analysed the star with fresh code, given only the S14 and S41 transits "
               "(independent/report.md). It found the S21 dip independently and narrowed P to the same two values. "
               "It measured Tc3 = 1883.767 &plusmn; 0.007 and depths of 3.3&ndash;3.6 ppt in SAP, QLP and its own photometry. "
               "Sliding the transit model over all 462 product-days, S21 is the only dip seen in two products. "
               "Its verdict: probably a real, on-target transit; the claim is 'narrowed to two periods', not 'resolved'.")


def fig(name, caption, width):
    from PIL import Image as PILImage
    with PILImage.open(HERE / name) as im:
        w, h = im.size
    # embed as palette PNG (smaller PDF, so it can go out as an email attachment)
    jpg = HERE / "_pdf_img" / (Path(name).stem + ".png")
    jpg.parent.mkdir(exist_ok=True)
    with PILImage.open(HERE / name) as im:
        im.convert("RGB").resize((int(im.width * 0.7), int(im.height * 0.7))).quantize(colors=64).save(jpg, optimize=True)
    return KeepTogether([Image(str(jpg), width=width, height=width * h / w), Paragraph(caption, cap)])


def table(rows, widths, header=True):
    data = [[Paragraph(str(c), cellb if (header and i == 0) else cell) for c in r] for i, r in enumerate(rows)]
    t = Table(data, colWidths=widths)
    st = [("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#c8d1db")),
          ("VALIGN", (0, 0), (-1, -1), "TOP"), ("TOPPADDING", (0, 0), (-1, -1), 2),
          ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]
    if header:
        st.append(("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef2f6")))
    t.setStyle(TableStyle(st))
    return t


def main():
    an = json.load(open(HERE / "analysis.json"))
    ep = json.load(open(HERE / "ephemeris.json"))
    s21 = json.load(open(HERE / "s21_check.json"))
    sec = json.load(open(HERE / "secondary.json"))
    maps = json.load(open(HERE / "maps.json"))
    vet = json.load(open(ROOT / "singles" / "results" / "vetting" / "159614322_n1.json"))
    pred = pd.read_csv(HERE / "predictions.csv")
    g = an["gaia"]
    tr = ep["transits"]
    e8, e4 = ep["ephemerides"]["n8"], ep["ephemerides"]["n4"]
    e8_3 = ep["ephemerides_3transit"]["n8"]
    s21c = ep["S21_check"]
    infl = ep["sigma_inflation"]
    rstar = g["radius_flame"]
    rp_rj = an["best_b_P90"]["rp"] * rstar / 0.10045
    rho_star = g["mass_flame"] / rstar ** 3 * 1.41
    b8 = {r["b"]: r for r in an["b_profile"]["8"]}
    b1 = {r["b"]: r for r in an["b_profile"]["1"]}
    fa = sorted(s21["false_alarm_orbit_ends"], key=lambda r: -r["drop_ppm"])
    n_fa = len(fa)
    second = next(r for r in fa if not (r['sector'] == 21 and abs(r['t_end'] - 1883.97) < 0.2))
    q21 = next(f for f in s21["fits"] if f["product"] == "QLP")
    own = s21["pixel"]

    doc = SimpleDocTemplate(str(HERE / "TIC159614322_writeup.pdf"), pagesize=letter,
                            leftMargin=0.7 * inch, rightMargin=0.7 * inch, topMargin=0.6 * inch,
                            bottomMargin=0.6 * inch, author="Jon McBee",
                            title="TIC 159614322: a third TESS transit narrows the period to 89.965 or 179.930 d")
    s = [Paragraph("CTOI 159614322.01: a third TESS transit in Sector 21 narrows a Jupiter-sized "
                   "duotransit to P = 89.965 d or 179.930 d", title), Spacer(1, 3),
         Paragraph("Jon McBee (independent researcher), jon.mcbee@gmail.com &middot; 2026-10-11 &middot; "
                   "analysis of public TESS and Gaia data", sub)]

    s.append(Paragraph("Summary", h2))
    s.append(Paragraph(
        "CTOI 159614322.01 was submitted in October 2021 (ExoFOP user <i>vogl</i>, \"in-house single "
        "transit finder\") as a single transit in Sector 41, with no period. It has no follow-up on "
        "ExoFOP. We find two more transits in archival TESS data:", body))
    s.append(Paragraph(
        f"&bull; <b>Sector 14</b> (BTJD {tr[0]['T']:.3f}): a full transit of the same depth and duration, "
        f"in the TESS-SPOC and QLP FFI light curves. It leaves P = {ep['dT']:.2f}/n.", body))
    s.append(Paragraph(
        f"&bull; <b>Sector 21</b> (BTJD {tr[1]['T']:.3f}): ingress and about 9 h of the transit bottom, "
        "up to the orbit-1 downlink gap. The egress is not covered. Our automated pipeline therefore "
        "left this window unclassified, and we examined it by hand.", body))
    s.append(Paragraph(
        "<b>Periods.</b> Flat TESS data exclude every n except 1, 2, 4, 5 and 8. The S21 event falls "
        "exactly where n = 4 and n = 8 predict a transit, and at phase 0.25 or 0.5 of the others.", body))
    s.append(Paragraph(
        f"<b>If S21 is a transit</b>, the period is <b>{e8['P']:.3f} d or {e4['P']:.3f} d</b>.", body))
    s.append(Paragraph(
        f"<b>The companion.</b> Its radius is about {rp_rj:.2f} R<sub>Jup</sub>, and the host is a "
        f"Gaia RV-stable G subgiant (Teff {g['teff_gspphot']:.0f} K, R* {rstar:.2f} R<sub>sun</sub>).", body))
    s.append(Paragraph(
        "<b>Next test.</b> Both periods predict a transit on <b>2026-12-28 at "
        f"{pred.iloc[0].mid_UTC[11:]} UTC &plusmn; {pred.iloc[0].sigma_plan_min} min</b>. The formal error is "
        f"&plusmn; {pred.iloc[0].sigma_min} min, inflated &times;{infl:.1f}; see the ephemeris note. The target is "
        "circumpolar from mid-northern latitudes.", body))
    s.append(table([
        ["Ephemeris (BJD_TDB)", "P = 89.965 d (n = 8)", "P = 179.930 d (n = 4)"],
        ["P (d)", f"{e8['P']:.4f} &plusmn; {e8['P_err']:.4f}", f"{e4['P']:.4f} &plusmn; {e4['P_err']:.4f}"],
        ["T0", f"{e8['T0'] + 2457000:.4f} &plusmn; {e8['T0_err']:.4f}", f"{e4['T0'] + 2457000:.4f} &plusmn; {e4['T0_err']:.4f}"],
        ["Cycles of S14 / S21 / S41", "0 / 2 / 8", "0 / 1 / 4"],
        ["S21 O&minus;C (min)", " / ".join(f"{r['O_minus_C_min']:+.0f}" for r in s21c["reductions"])
         + " (TESS-SPOC SAP / QLP / TESS-SPOC PDC)", "same"]],
        [1.7 * inch, (W - 1.7 * inch) / 2, (W - 1.7 * inch) / 2]))
    s.append(Spacer(1, 2))
    s.append(Paragraph(
        "<b>How the ephemeris was made.</b> P and T0 come from the two full transits (S14, S41) only. "
        "The partial S21 transit is used only to select n = 4 or 8, and as a consistency check. Transit times (BTJD), "
        "each fitted with one transit shape plus a local linear baseline: "
        + "; ".join(f"S{r['sector']} {r['T']:.4f} &plusmn; {r['T_err']:.4f} ({r['product']} {'SAP' if r['column'] == 'sap_flux' else 'PDC'})" for r in tr)
        + f". Against the S14&ndash;S41 ephemeris, S21 arrives {min(r['O_minus_C_min'] for r in s21c['reductions']):+.0f} to "
        f"{max(r['O_minus_C_min'] for r in s21c['reductions']):+.0f} min late "
        f"({min(r['sigma'] for r in s21c['reductions']):.1f}&ndash;{max(r['sigma'] for r in s21c['reductions']):.1f}&sigma;), "
        "possibly because it is a partial transit. A fit through all three transits has &chi;<super>2</super> = "
        f"{e8_3['chi2']:.1f} for 1 dof. To allow for underestimated timing errors, predicted transit times carry a "
        f"planning uncertainty of the formal error &times; &radic;&chi;<super>2</super> = &times;{infl:.1f}.", cap))
    s.append(fig("fig_transits.png", "Figure 1. The three transits, with one transit shape (orange) "
                 "scaled to each event's fitted depth. Sector 21 is shown in three reductions: TESS-SPOC SAP, "
                 "QLP, and our own background-subtracted TESScut photometry (2-px radius). The data end at the "
                 "downlink gap just before the predicted egress.", W))
    s.append(fig("fig_aliases.png", f"Figure 2. Every period P = {ep['dT']:.2f}/n (n = 1&ndash;12) against TESS coverage "
                 "(grey: sectors observed). A red cross marks a predicted transit on flat data, which excludes that n. "
                 "Only n = 1, 2, 4, 5 and 8 survive. n = 4 and 8 predict the S21 event; n = 1, 2 and 5 do not.", W))

    s.append(Paragraph("Is the Sector 21 event real?", h2))
    s21rows = [["Check", "Result"]]
    s21rows.append(["Independent reductions",
                    "TESS-SPOC SAP: depth scale "
                    + f"{tr[1]['depth_scale']:.2f} &plusmn; {tr[1]['depth_scale_err']:.2f}. QLP: {q21['depth_scale']:.2f} &plusmn; {q21['depth_scale_err']:.2f}. "
                    + "Own TESScut photometry with a local background: 4.3, 3.3, 3.5, 3.2, 2.8 and 2.7 ppt for apertures of radius 1&ndash;4 px "
                    + "(4&ndash;52 pixels). Each is &ge; 8&sigma;; S14/S41 depth is about 3.3 ppt."])
    s21rows.append(["Timing and shape", f"Sharp ingress at BTJD &asymp; 1883.57, which n = 4/8 predict at 1883.58. "
                    f"&Delta;&chi;<super>2</super> (transit vs. no transit, same baseline freedom): {tr[1] and next(f for f in s21['fits'] if f['product']=='TESS-SPOC' and f['column']=='sap_flux')['dchi2_vs_linear']:.0f} "
                    f"(TESS-SPOC SAP), {q21['dchi2_vs_linear']:.0f} (QLP); against a quadratic baseline "
                    f"{next(f for f in s21['fits'] if f['product']=='TESS-SPOC' and f['column']=='sap_flux')['dchi2_vs_quadratic']:.0f} and {q21['dchi2_vs_quadratic']:.0f}."])
    s21rows.append(["Background / scattered light", "The TESScut background is flat at 74.2&ndash;74.9 e/s/px from the ingress to 1883.89. "
                    "It rises (scattered light) only in the last 4 cadences before the gap."])
    s21rows.append(["Pointing", "No step at the ingress. The centroid drifts smoothly by &lt; 0.01 px across the window."])
    s21rows.append(["Systematics at orbit ends", f"The same test at the end of all {n_fa} orbits observed for this star "
                    f"(all products) finds S21 the largest by far. The largest drop at any other orbit end is {second['drop_ppm']:.0f} ppm "
                    f"(S{second['sector']} {second['product']}, {second['snr']:.1f}&sigma;)."])
    s21rows.append(["Source location", f"Transit-amplitude maps (Fig. 3) put the source on the target in S14 "
                    f"({maps['S14']['offset_arcsec']:.1f}\" &plusmn; {maps['S14']['offset_err_arcsec']:.1f}\") and S41 "
                    f"({maps['S41']['offset_arcsec']:.1f}\" &plusmn; {maps['S41']['offset_err_arcsec']:.1f}\"). "
                    f"In S21 the fit gives {maps['S21']['offset_arcsec']:.1f}\" &plusmn; {maps['S21']['offset_err_arcsec']:.1f}\" "
                    f"(SNR {maps['S21']['peak_snr']:.1f}; 43 frames, no post-egress baseline). No catalogued star can host the "
                    "signal there: making 3 ppt needs G &lt; 17.9, and the nearest Gaia star that bright is 39.8\" away (G = 14.5). "
                    "The aperture test (depth largest in the smallest aperture) also points at the target."])
    if INDEPENDENT:
        s21rows.append(["Blind re-analysis", INDEPENDENT])
    s.append(table(s21rows, [1.45 * inch, W - 1.45 * inch]))
    s.append(Spacer(1, 4))
    s.append(fig("fig_maps.png", "Figure 3. Per-pixel transit-amplitude maps (each pixel: quadratic baseline + transit "
                 "shape). Blue circle: target; green cross: fitted source position; white squares: Gaia stars with G &lt; 16.", 0.75 * W))

    s.append(Paragraph("Star and companion", h2))
    s.append(table([
        ["Quantity", "Value", "Source"],
        ["TIC / Gaia DR3", f"159614322 / {g['source_id']}", "ExoFOP, Gaia DR3"],
        ["RA, Dec (J2000)", "220.00977, +65.48404", "TIC"],
        ["Tmag / G / BP&minus;RP", f"11.18 / {g['phot_g_mean_mag']:.2f} / {g['bp_rp']:.2f}", "TIC, Gaia DR3"],
        ["Distance", f"{1000 / g['parallax']:.0f} pc (parallax {g['parallax']:.3f} &plusmn; {g['parallax_error']:.3f} mas)", "Gaia DR3"],
        ["Teff / log g / [M/H]", f"{g['teff_gspphot']:.0f} K / {g['logg_gspphot']:.2f} / {g['mh_gspphot']:+.2f} (GSP-Phot)", "Gaia DR3"],
        ["R* / M* / age", f"{rstar:.2f} R<sub>sun</sub> / {g['mass_flame']:.2f} M<sub>sun</sub> / {g['age_flame']:.1f} Gyr "
         f"(&rho;* &asymp; {rho_star:.2f} g cm<super>&minus;3</super>)", "Gaia DR3 FLAME"],
        ["Binarity", f"RUWE {g['ruwe']:.2f} (mildly elevated); ipd_frac_multi_peak {g['ipd_frac_multi_peak']:.0f}; non_single_star 0; "
         f"RV {g['radial_velocity']:.1f} &plusmn; {g['radial_velocity_error']:.1f} km/s, not variable "
         f"(rv_chisq_pvalue {g['rv_chisq_pvalue']:.2f}, {g['rv_nb_transits']:.0f} transits)", "Gaia DR3"],
        ["Companion", f"R<sub>p</sub>/R* = {an['best_b_P90']['rp']:.4f}, so R<sub>p</sub> &asymp; {rp_rj:.2f} R<sub>Jup</sub>; "
         f"T14 &asymp; 9.7 h; depth &asymp; 3.3 ppt", "this work"]],
        [1.5 * inch, 4.1 * inch, W - 5.6 * inch]))

    s.append(Paragraph("Vetting of the S14 and S41 transits", h2))
    sbp = vet["scale_by_product"]
    s.append(table([
        ["Check", "Result"],
        ["Depth consistency", f"Per-epoch depth scale S14 0.99 &plusmn; 0.05, S41 1.00 &plusmn; 0.03 (&chi;<super>2</super> {vet['depth_consistency_chi2']:.2f}/1)"],
        ["Products", "Depth scale: " + ", ".join(f"{k.replace(':', ' ')} {v[0]:.2f}" for k, v in sbp.items())
         + ". A separate detrending code re-detects both events at SNR 19."],
        ["Centroids (SPOC)", f"No in-transit shift: {vet['mom_centr1_shift_pix'][2]:.1f}&sigma; and {vet['mom_centr2_shift_pix'][2]:.1f}&sigma;"],
        ["Data quality", "No momentum dumps or segment edges near either transit"],
        ["Secondary eclipse", f"P = 89.965 d: {sec['n8']['mean_depth_ppm']:.0f} &plusmn; {sec['n8']['mean_err_ppm']:.0f} ppm at phase 0.5 "
         "(S49, S75). P = 179.930 d: phase 0.5 never observed."],
        ["Transit shape vs. star", f"Circular orbit, b = 0: &rho;<sub>circ</sub> = {b8[0.0]['rho_circ']:.2f} g cm<super>&minus;3</super> at P = 89.97 d. "
         f"The FLAME density ({rho_star:.2f}) needs b &asymp; 0.7 (&Delta;&chi;<super>2</super> &asymp; {b8[0.7]['chi2'] - b8[0.0]['chi2']:.0f}) "
         "or e of at least ~0.3. At 179.93 d it needs b &asymp; 0.85 or e of at least ~0.45. So the duration mildly favours the shorter period "
         "for a circular orbit; eccentricity relaxes this."],
        ["Catalogues", "Not a TOI; not in the TESS EB catalogue (Pr&scaron;a et al. 2022); no confirmed planet; ExoFOP has "
         "no imaging, spectra or notes. A web search found no other mention of this TIC."]],
        [1.45 * inch, W - 1.45 * inch]))

    s.append(Paragraph("Predicted transits", h2))
    p8 = pred[pred.alias.str.contains("n=8")]
    p4c = set(pred[pred.alias.str.contains("n=4")].T_btjd.round(3))
    rows = [["Mid-transit (UTC)", "&sigma; plan (formal), min", "Predicted by", "Fraction observable* Boston / Kitt Peak / La Palma"]]
    for _, r in p8.head(7).iterrows():
        both = round(r["T_btjd"], 3) in p4c
        rows.append([r["mid_UTC"], f"{r['sigma_plan_min']} ({r['sigma_min']})", "both periods" if both else "89.965 d only",
                     f"{r['vis_Boston, MA']:.0%} / {r['vis_Kitt Peak, AZ']:.0%} / {r['vis_La Palma']:.0%}"])
    s.append(table(rows, [1.3 * inch, 1.2 * inch, 1.1 * inch, W - 3.6 * inch]))
    s.append(Paragraph("* Fraction of the transit &plusmn; 1 h observable with the target above 20&deg; and the Sun below "
                       "&minus;12&deg;. Transits on the 179.930-d ephemeris are every other one in the list. "
                       "tess-point lists Sectors 120&ndash;122 and 126&ndash;128 as the next TESS visits (2027&ndash;28).", cap))

    s.append(Paragraph("Suggested follow-up", h2))
    s.append(Paragraph(
        "(1) <b>Photometry on 2026-12-28.</b> A detection confirms the S21 transit and both "
        "periods; a flat light curve would point back to the 720/360/144-d aliases. The transit "
        "is long (~9.7 h) and the window is &plusmn;45 min, so two longitudes and a &ge;0.4 m telescope help. "
        "From New England the star is at only 17&ndash;20&deg; at ingress (observing_plan_providence.md).", body))
    s.append(Paragraph(
        "(2) <b>2027-03-28.</b> Only the 89.965-d period predicts a transit on this date, so it "
        "separates the two periods.", body))
    s.append(Paragraph(
        "(3) <b>RVs.</b> At V &asymp; 11.9, a Jupiter-mass companion at 90&ndash;180 d gives K &asymp; "
        "34&ndash;42 m/s. A stellar companion would show km/s variations.", body))
    s.append(Paragraph(
        "(4) <b>High-resolution imaging.</b> Gaia RUWE is 1.51, so imaging is worth doing.", body))
    s.append(Spacer(1, 6))
    s.append(Paragraph(
        "Methods: the repeat search and alias pipeline were validated blind on 26 TESS duotransit hosts "
        "with known periods (true period kept 26/26, no wrong unique answers). Code and intermediate files "
        "are available on request. Data: NASA TESS (SPOC, TESS-SPOC, QLP via MAST; TESScut) and ESA Gaia DR3. "
        "This analysis was done with substantial help from an AI system (Anthropic's Claude); the author has "
        "checked the results against the saved data and figures.", sub))
    doc.build(s)
    print(HERE / "TIC159614322_writeup.pdf")


if __name__ == "__main__":
    main()
