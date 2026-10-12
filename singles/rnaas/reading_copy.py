"""
Reading copy of the RNAAS note: rnaas.tex rendered as a plain, readable PDF, with no LaTeX
needed. It lets the author read every sentence before submission (the AAS AI policy requires
that), with the AI disclosure statement highlighted.

The layout is NOT the journal's. The final PDF comes from Overleaf (aastex7.cls).
Run:    python singles/rnaas/reading_copy.py
Output: singles/rnaas/rnaas_reading_copy.pdf
"""
import re
from pathlib import Path

from pylatexenc.latex2text import LatexNodes2Text
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

from reportlab.lib.fonts import addMapping
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# Arial (Windows) covers accented names (Sipőcz) and ⊕; fall back to Helvetica elsewhere
try:
    for name, f in (("Arial", "arial.ttf"), ("Arial-Bold", "arialbd.ttf"), ("Arial-Italic", "ariali.ttf")):
        pdfmetrics.registerFont(TTFont(name, f"C:/Windows/Fonts/{f}"))
    addMapping("Arial", 0, 0, "Arial")
    addMapping("Arial", 1, 0, "Arial-Bold")
    addMapping("Arial", 0, 1, "Arial-Italic")
    addMapping("Arial", 1, 1, "Arial-Bold")
    FONT, FONT_B, FONT_I = "Arial", "Arial-Bold", "Arial-Italic"
except Exception:  # noqa: BLE001
    FONT, FONT_B, FONT_I = "Helvetica", "Helvetica-Bold", "Helvetica-Oblique"

HERE = Path(__file__).resolve().parent
TEX = (HERE / "rnaas.tex").read_text(encoding="utf-8")
BS = chr(92)
INK, MUTED, ACCENT = colors.HexColor("#1d2b3a"), colors.HexColor("#5d6b7c"), colors.HexColor("#d9682b")

ttl = ParagraphStyle("t", fontName=FONT_B, fontSize=14, leading=18, textColor=INK, spaceAfter=4)
sub = ParagraphStyle("s", fontName=FONT, fontSize=9, leading=12, textColor=MUTED, spaceAfter=8)
h2 = ParagraphStyle("h", fontName=FONT_B, fontSize=11, leading=14, textColor=ACCENT, spaceBefore=8, spaceAfter=3)
body = ParagraphStyle("b", fontName=FONT, fontSize=9.6, leading=13.2, textColor=INK, spaceAfter=5)
small = ParagraphStyle("sm", fontName=FONT, fontSize=8, leading=10.5, textColor=INK, spaceAfter=2)
cell = ParagraphStyle("c", fontName=FONT, fontSize=7.8, leading=9.6, textColor=INK)
note = ParagraphStyle("n", fontName=FONT_I, fontSize=8.5, leading=11, textColor=MUTED, spaceAfter=6)

L2T = LatexNodes2Text(math_mode="text", strict_latex_spaces=True)


def txt(s):
    """LaTeX fragment -> plain text -> ReportLab-safe markup."""
    s = re.sub(r"\\citep\[([^\]]*)\]\[\]\{([^}]*)\}", r"(\1 \2)", s)
    s = re.sub(r"\\citep\[([^\]]*)\]\{([^}]*)\}", r"(\1 \2)", s)
    s = re.sub(r"\\citep\{([^}]*)\}", r"(\1)", s)
    s = re.sub(r"\\citet\{([^}]*)\}", r"\1", s)
    s = re.sub(r"\\ref\{[^}]*\}", "1", s)
    s = re.sub(r"\\url\{([^}]*)\}", r"\1", s)
    s = re.sub(r"\\textit\{([^}]*)\}", r"@@I\1@@/I", s)
    s = re.sub(r"\\texttt\{([^}]*)\}", lambda m: m.group(1).replace(BS + "_", "@@US@@"), s)
    s = re.sub(r"\$\^\{?([ab,]+)\}?\$", r"(\1)", s)
    out = L2T.latex_to_text(s)
    out = re.sub(r"\s+", " ", out).strip()
    out = out.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    out = out.replace("\u2299", "sun").replace("\u2273", "&gt;~").replace("\u2295", "Earth").replace("\u2208", " in ")
    out = re.sub(r"([A-Za-z])_\{?([A-Za-z0-9\u2295]+)\}?", r"\1<sub>\2</sub>", out)
    out = re.sub(r"\^\{?(-?[0-9]+)\}?", r"<super>\1</super>", out)
    out = out.replace("@@US@@", "_")
    # citation keys like Prsa2022 -> Prsa 2022
    out = re.sub(r"\b([A-Z][A-Za-z]+)(\d{4})\b", r"\1 \2", out)
    return out.replace("@@I", "<i>").replace("@@/I", "</i>")


def between(a, b, s=TEX):
    i = s.index(a) + len(a)
    return s[i:s.index(b, i)]


def main():
    story = []
    title = between(BS + "title{", "}\n")
    story.append(Paragraph(txt(title), ttl))
    story.append(Paragraph("Jon McBee (ORCID 0009-0000-7524-4470) &middot; Independent researcher, USA &middot; "
                           "jon.mcbee@gmail.com", sub))
    story.append(Paragraph("READING COPY generated from rnaas.tex. The words are the submission text; the "
                           "layout is not the journal's.", note))
    story.append(Paragraph("Abstract", h2))
    story.append(Paragraph(txt(between(BS + "begin{abstract}", BS + "end{abstract}")), body))
    kw = re.findall(r"\\uat\{([^}]*)\}\{(\d+)\}", TEX)
    story.append(Paragraph("<b>Keywords:</b> " + "; ".join(f"{a} ({b})" for a, b in kw), small))

    main_txt = TEX[TEX.index(BS + "section{Sample}"):TEX.index(BS + "begin{acknowledgments}")]
    table_tex = re.search(r"\\begin\{deluxetable\*\}.*?\\end\{deluxetable\*\}", main_txt, re.S).group(0)
    main_txt = main_txt.replace(table_tex, "@@TABLE@@")
    parts = re.split(r"\\section\{([^}]*)\}", main_txt)
    for k in range(1, len(parts), 2):
        story.append(Paragraph(txt(parts[k]), h2))
        sec = parts[k + 1]
        for block in re.split(r"\n\s*\n", sec):
            if "@@TABLE@@" in block:
                before, after = block.split("@@TABLE@@")
                if before.strip():
                    story.append(Paragraph(txt(before), body))
                story.append(make_table(table_tex))
                if after.strip():
                    story.append(Paragraph(txt(after), body))
                continue
            if BS + "begin{itemize}" in block:
                pre, rest = block.split(BS + "begin{itemize}", 1)
                items, post = rest.split(BS + "end{itemize}", 1)
                if pre.strip():
                    story.append(Paragraph(txt(pre), body))
                for it in items.split(BS + "item")[1:]:
                    story.append(Paragraph("&bull; " + txt(it), body))
                if post.strip():
                    story.append(Paragraph(txt(post), body))
                continue
            if block.strip():
                story.append(Paragraph(txt(block), body))

    story.append(Paragraph("Acknowledgments", h2))
    story.append(Paragraph(txt(between(BS + "begin{acknowledgments}", BS + "end{acknowledgments}")), body))
    story.append(Paragraph("<b>Facilities:</b> " + txt(between(BS + "facilities{", "}\n")), small))
    story.append(Paragraph("<b>Software:</b> " + txt(between(BS + "software{", "}\n")), small))

    ai = between(BS + "section*{AI disclosure statement}", BS + "begin{thebibliography}")
    box = Table([[Paragraph("<b>AI disclosure statement</b>", ParagraphStyle("x", parent=h2, spaceBefore=0))],
                 [Paragraph(txt(ai), body)]], colWidths=[letter[0] - 1.6 * inch])
    box.setStyle(TableStyle([("BOX", (0, 0), (-1, -1), 1.2, ACCENT),
                             ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#fdf3ec")),
                             ("LEFTPADDING", (0, 0), (-1, -1), 8), ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                             ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5)]))
    story.append(Spacer(1, 8))
    story.append(box)

    story.append(Paragraph("References", h2))
    for b in re.findall(r"\\bibitem\[[^\]]*\]\{[^}]*\}\s*(.*)", TEX):
        story.append(Paragraph(txt(b), small))

    doc = SimpleDocTemplate(str(HERE / "rnaas_reading_copy.pdf"), pagesize=letter, leftMargin=0.8 * inch,
                            rightMargin=0.8 * inch, topMargin=0.7 * inch, bottomMargin=0.7 * inch,
                            title="RNAAS note: reading copy", author="Jon McBee")
    doc.build(story)
    print(HERE / "rnaas_reading_copy.pdf")


def make_table(t):
    cap = txt(re.search(r"\\tablecaption\{(.*?)\\label", t, re.S).group(1))
    heads = [txt(h) for h in re.findall(r"\\colhead\{([^}]*)\}", t)]
    data = t.split(BS + "startdata")[1].split(BS + "enddata")[0]
    rows = [[Paragraph(txt(c), cell) for c in r.split("&")] for r in data.strip().split(BS + BS) if r.strip()]
    comm = txt(re.search(r"\\tablecomments\{(.*)\}\s*\\end\{deluxetable", t, re.S).group(1).replace("$^a$", "(a) ").replace("$^b$", " (b) "))
    tb = Table([[Paragraph(f"<b>{h}</b>", cell) for h in heads]] + rows,
               colWidths=[0.85 * inch, 1.35 * inch, 2.1 * inch, 2.6 * inch])
    tb.setStyle(TableStyle([("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#c8d1db")),
                            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef2f6")),
                            ("VALIGN", (0, 0), (-1, -1), "TOP")]))
    return Table([[Paragraph(f"<b>Table 1.</b> {cap}", small)], [tb], [Paragraph(comm, note)]],
                 colWidths=[letter[0] - 1.6 * inch])


if __name__ == "__main__":
    main()
