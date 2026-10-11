"""
Technical write-up of TIC 147035057 / CTOI 147035057.01 for the observers who can verify it.
Run:  python singles/writeups/tic147035057/make_figures.py && python singles/writeups/tic147035057/build_pdf.py
Output: singles/writeups/tic147035057/TIC147035057_writeup.pdf
Numbers trace to singles/results/vetting/147035057_* and singles/results/pilot/147035057*.
"""
import json
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (Image, KeepTogether, Paragraph, SimpleDocTemplate, Spacer, Table,
                                TableStyle)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
VET = ROOT / "singles" / "results" / "vetting"
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


def fig(name, caption, width):
    from PIL import Image as PILImage
    with PILImage.open(HERE / name) as im:
        w, h = im.size
    return KeepTogether([Image(str(HERE / name), width=width, height=width * h / w), Paragraph(caption, cap)])


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
    e = json.load(open(VET / "147035057_ephemeris.json"))
    nxt = pd.read_csv(VET / "147035057_next_eclipses.csv")
    T0b, P = e["T0_btjd"], e["P"]
    doc = SimpleDocTemplate(str(HERE / "TIC147035057_writeup.pdf"), pagesize=letter,
                            leftMargin=0.7 * inch, rightMargin=0.7 * inch, topMargin=0.6 * inch,
                            bottomMargin=0.6 * inch, title="TIC 147035057: 37.315-d period from three TESS eclipses",
                            author="Jon McBee")
    s = [Paragraph("CTOI 147035057.01 (TYC 5963-01913-1): three TESS eclipses fix the period at "
                   "37.315 d, but Gaia flags the host as RV-variable", title), Spacer(1, 3),
         Paragraph("Jon McBee (independent researcher), jon.mcbee@gmail.com &middot; 2026-10-11 &middot; "
                   "analysis of public TESS and Gaia data", sub)]

    s.append(Paragraph("Summary", h2))
    s.append(Paragraph(
        "CTOI 147035057.01 was submitted in February 2021 (ExoFOP user <i>steuer</i>, "
        "\"in-house single transit finder\") as a single transit in Sector 33. It has no follow-up "
        "on ExoFOP. TESS shows two more eclipses of the same depth (~3.2 ppt) and duration (~4.5 h):", body))
    s.append(Paragraph(
        "&bull; <b>Sector 7:</b> an earlier eclipse at BTJD 1515.065, in the TESS-SPOC and QLP "
        "light curves.", body))
    s.append(Paragraph(
        "&bull; <b>Sector 108:</b> a new eclipse at BTJD 4313.689 (2026-09-30). We found it with our "
        "own aperture photometry on the TESScut FFIs, because S108 has no SPOC/QLP light curve yet.", body))
    s.append(Paragraph(
        "The three events fit a linear ephemeris (O&minus;C &le; 1.3 min, &chi;<super>2</super> = 0.97 "
        "for 1 dof), and every shorter alias is excluded by flat TESS data. Of the four aliases that "
        "predicted an S108 eclipse, only P = 37.315 d shows one. Difference images put the source on "
        "the target.", body))
    s.append(Paragraph(
        "<b>Caveat on the companion:</b> Gaia DR3 reports large RV variability for this star "
        "(rv_chisq_pvalue = 0, rv_renormalised_gof = 99.8, rv_amplitude_robust = 49.8 km/s over "
        "27 epochs), and we see no secondary eclipse at phase 0.5. The companion may be a low-mass "
        "star on an eccentric orbit, or a substellar/planetary body in a binary system. RV "
        "monitoring will decide.", body))
    s.append(table([
        ["Ephemeris (BJD_TDB)", "Value"],
        ["P", f"{P:.5f} &plusmn; {e['P_err']:.5f} d"],
        ["T0", f"{T0b + 2457000:.4f} &plusmn; {e['T0_err']:.4f}  (BTJD {T0b:.4f})"],
        ["Eclipse times (BTJD)", " ; ".join(f"{t:.4f} &plusmn; {er:.4f} (cycle {int(k)})"
                                           for t, er, k in zip(e["epochs"], e["epoch_err"], e["cycles"]))],
        ["Depth / duration", "3.19, 3.28 ppt (S7, S33); 2.8 ppt (S108, own photometry); T14 &asymp; 4.5 h"],
        ["Next eclipse well placed from Paranal",
         f"{nxt[nxt.observable_Paranal].iloc[0].mid_UTC} UTC &plusmn; {nxt[nxt.observable_Paranal].iloc[0].sigma_min} min "
         f"(alt {nxt[nxt.observable_Paranal].iloc[0].alt_at_Paranal_deg}&deg;)"]],
        [1.7 * inch, W - 1.7 * inch]))
    s.append(Spacer(1, 4))
    s.append(fig("fig_eclipses.png", "Figure 1. The three eclipses, detrended locally, with one common transit-shape "
                 "model (orange). Grey: individual cadences; dark: 30-min bins (or all points for S7's 30-min data).", W))
    s.append(fig("fig_s108.png", "Figure 2. Sector 108 (own TESScut aperture photometry, 200 s). Vertical lines: "
                 "predicted eclipse times for every alias still allowed before S108. Only P = 37.315 d shows an eclipse; "
                 "P = 50.64 d falls on flat data (excluded), and P = 54.54 and 59.08 d fall on a gap edge and a ramp.", W))

    s.append(Paragraph("Star", h2))
    s.append(table([
        ["Identifier / quantity", "Value", "Source"],
        ["TIC / TYC / Gaia DR3", "147035057 / 5963-01913-1 / 2947647054662200192", "ExoFOP"],
        ["RA, Dec (J2000)", "104.27440, &minus;16.47848", "TIC"],
        ["Tmag / G / BP&minus;RP", "8.49 / 8.83 / 0.62", "TIC, Gaia DR3"],
        ["Distance", "166 pc (parallax 6.023 &plusmn; 0.017 mas)", "Gaia DR3"],
        ["Teff / log g", "6400 K / 3.90 (GSP-Phot); 6480 K (TIC)", "Gaia DR3, TIC"],
        ["R* / M* / age", "2.00 R<sub>sun</sub> / 1.44 M<sub>sun</sub> / 2.4 Gyr", "Gaia DR3 FLAME"],
        ["Astrometry", "RUWE 1.12, ipd_frac_multi_peak 0, no NSS orbit", "Gaia DR3"],
        ["Radial velocity", "&minus;42.2 &plusmn; 4.2 km/s; rv_chisq_pvalue 0; rv_renormalised_gof 99.8; "
                            "rv_amplitude_robust 49.8 km/s (27 epochs); vbroad 25.9 km/s", "Gaia DR3"],
        ["Companion size if on primary", "&asymp; 12.5 R<sub>Earth</sub> (1.1 R<sub>Jup</sub>) for 3.2 ppt; larger if grazing", "this work"]],
        [1.6 * inch, 3.9 * inch, W - 5.5 * inch]))

    s.append(Paragraph("Vetting", h2))
    s.append(table([
        ["Check", "Result"],
        ["Depth consistency", "S7/S33 scale factors 0.997 &plusmn; 0.015 / 1.003 &plusmn; 0.015 (&chi;<super>2</super> 0.08/1); "
         "S108 0.89 in our photometry"],
        ["Independent reductions", "QLP, TESS-SPOC and SPOC 2-min depths agree within 4%; a separate detrending code gives 2.97 / 3.09 ppt"],
        ["Difference images", "S33: 2.9\" &plusmn; 1.0\" from target (SNR 15); S108: 2.2\" &plusmn; 0.8\" (SNR 17); "
         "S7: 6.3\" &plusmn; 1.1\" (0.3 px, Cycle-1 30-min FFI). Validated on TOI-5542 b and TOI-2589 b (on target within 1&ndash;3\")"],
        ["Aperture test (S108)", "Depth 2.77 / 2.78 / 2.69 ppt for 3 / 9 / 21-pixel apertures; background flat (165.7 vs 165.8 e/s)"],
        ["Centroids (SPOC)", "No in-transit shift: 1.1&sigma; and 0.8&sigma;"],
        ["Data quality", "No momentum dumps, quality flags or segment edges near any eclipse"],
        ["Neighbours", "Brightest within 40\": two G = 13.5 stars at 17\" and 24\". Either would need a ~25% eclipse, "
         "and the difference images do not favour them"],
        ["Secondary eclipse", "None at phase 0.5 (&minus;0.05 &plusmn; 0.07 ppt). No significant dip at other phases; "
         "the largest features are ~0.2 ppt, at the noise level"],
        ["Catalogues", "Not in the TESS EB catalogue (Pr&scaron;a et al. 2022), not a TOI, no confirmed planet; "
         "not in Hawthorn et al. (2024) despite fitting its S1&ndash;13 + S27&ndash;39 duotransit pattern"]],
        [1.45 * inch, W - 1.45 * inch]))
    s.append(Spacer(1, 4))
    s.append(fig("fig_diffimg.png", "Figure 3. Difference images (out-of-eclipse minus in-eclipse) at the S33 and S108 "
                 "eclipses. Blue circle: target; white squares: Gaia stars with G &lt; 15.", 0.62 * W))

    s.append(Paragraph("Method (brief)", h2))
    s.append(Paragraph(
        "<b>Repeat search:</b> repeats were found by testing a box of the catalogued depth and duration on a "
        "D/4 grid over all SPOC 2-min, TESS-SPOC and QLP light curves. Each window has its own "
        "segment-aware polynomial baseline.", body))
    s.append(Paragraph(
        "<b>Validation:</b> run blind from the first transit only, on 26 TESS duotransit hosts with "
        "known periods, the method kept the true period every time (26/26) and recovered it uniquely "
        "for 19, with no wrong unique answers.", body))
    s.append(Paragraph(
        "<b>Alias test:</b> aliases P = &Delta;T/n were tested against flat-window exclusions in all data.", body))
    s.append(Paragraph(
        "<b>Sector 108 photometry:</b> TESScut cutout, background-subtracted aperture photometry. As a "
        "check, the same code recovers the known S33 eclipse at 3.21 &plusmn; 0.07 ppt.", body))
    s.append(Paragraph(
        "<b>Mid-times:</b> each eclipse was fitted with a common batman transit shape and a linear "
        "ephemeris by weighted least squares. Code and intermediate files are available on request.", body))

    s.append(Paragraph("Suggested follow-up", h2))
    s.append(Paragraph(
        "(1) <b>RV monitoring.</b> A few epochs over ~40 days would show whether the Gaia RV "
        "variability follows the 37.3-d period, with a stellar-mass semi-amplitude of tens of km/s.", body))
    s.append(Paragraph(
        "(2) <b>Photometry of the 2027-01-20 eclipse.</b> It is well placed for southern facilities, "
        "and a depth comparison in two bands tests for colour dependence.", body))
    nx = nxt.head(6)
    s.append(table([["Cycle", "Mid-eclipse (UTC)", "&sigma; (min)", "Alt. at Paranal", "Sun alt."]] +
                   [[int(r.cycle), r.mid_UTC, r.sigma_min, f"{r.alt_at_Paranal_deg}&deg;", f"{r.sun_alt_deg}&deg;"]
                    for r in nx.itertuples()],
                   [0.6 * inch, 1.6 * inch, 0.8 * inch, 1.1 * inch, 0.9 * inch]))
    s.append(Spacer(1, 6))
    s.append(Paragraph(
        "Data: NASA TESS (SPOC, TESS-SPOC, QLP via MAST; FFI cutouts via TESScut) and ESA Gaia DR3. "
        "This analysis was carried out with substantial assistance from an AI system (Anthropic's Claude); "
        "the author has checked the results against the saved data and figures.", sub))
    doc.build(s)
    print(HERE / "TIC147035057_writeup.pdf")


if __name__ == "__main__":
    main()
