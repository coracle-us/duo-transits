"""
Word count as the AAS quality check does it for RNAAS (learned from the AAS82207R1 bounce, 2026-10-13): everything
after the abstract except the reference list counts against the 1350-word limit, i.e. body + table + acknowledgments
+ facilities/software + AI statement. The abstract (limit 150) is counted separately.
Calibration: on the AAS82207R1 text (tag aas82207R1-submitted) AAS reported 1706 words; this script's estimate is
printed for comparison.
Update (AAS84041 bounce, 2026-11-24): AAS now says the 1350 words include the references as well ("text,
acknowledgments, figure/table captions and references"). On the AAS84041 text (tag aas84041-submitted) AAS reported
the main text as 300 over the limit; this script gives 1662 with references (1243 without). The "with references"
figure is the one to keep under 1350 (it also counts the \bibitem keys, about 80 words that are not in the PDF).
Usage:  python singles/rnaas/aas_wordcount.py FILE.tex [FILE.tex ...]
"""
import re
import sys
from pathlib import Path


def words(s):
    s = re.sub(r"(?<!\\)%.*", "", s)                                   # comments
    s = re.sub(r"\\(begin|end)\{[^}]*\}(\{[^}]*\})*", " ", s)           # environment markers (+ column spec)
    s = re.sub(r"\\(label|ref|cite[tp]?|citet|citep|url|dataset)\[[^\]]*\]", r"\\\1", s)
    s = re.sub(r"\\(label|tablenum)\{[^}]*\}", " ", s)
    s = re.sub(r"\\(dataset)\{[^}]*\}", " ", s)                         # the visible text is the [..] doi label
    s = re.sub(r"\\[a-zA-Z]+\*?", " ", s)                               # command names (arguments are kept)
    s = re.sub(r"[{}$&~^_]|\\\\|\\[,;!]", " ", s)
    return len(re.findall(r"[A-Za-z0-9]", " ".join(re.findall(r"\S+", s))) and re.findall(r"\S*[A-Za-z0-9]\S*", s))


def count(path):
    t = Path(path).read_text(encoding="utf-8")
    a0, a1 = t.index("\\begin{abstract}"), t.index("\\end{abstract}")
    end = t.index("\\begin{thebibliography}") if "\\begin{thebibliography}" in t else t.index("\\end{document}")
    abstract = words(t[a0 + len("\\begin{abstract}"):a1])
    rest = words(t[a1:end])
    body = words(t[a1:t.index("\\begin{deluxetable")]) if "\\begin{deluxetable" in t else rest
    full = words(t[a1:t.index("\\end{document}")])
    print(f"{path}: abstract {abstract} (limit 150); text without references {rest}; "
          f"AAS-counted text with references {full} (limit 1350)")
    return abstract, full


if __name__ == "__main__":
    for p in sys.argv[1:]:
        count(p)
