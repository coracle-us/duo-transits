"""
Sync the generated table rows into k2singles/rnaas/rnaas.tex and count words (RNAAS limits:
1,350 words of text and a 150-word abstract; the table, acknowledgments and references excluded).
Usage:  python k2singles/scripts/rnaas_wordcount.py
"""
import re
from pathlib import Path

RN = Path(__file__).resolve().parents[1] / "rnaas"


def words(s):
    s = re.sub(r"(?<!\\)%.*", "", s)   # strip comments, not escaped \%
    s = re.sub(r"\\cite[tp]\{[^}]*\}", "CITE", s)
    s = re.sub(r"\$[^$]*\$", "X", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", "", s)
    s = re.sub(r"[{}~]", " ", s)
    return len([w for w in s.split() if re.search(r"[A-Za-z0-9]", w)])


def main():
    p = RN / "rnaas.tex"
    tex = p.read_text()
    rows = (RN / "table_rows.tex").read_text()
    i = tex.index("\\startdata\n") + len("\\startdata\n")
    j = tex.index("\\enddata")
    tex = tex[:i] + rows + tex[j:]
    p.write_text(tex)
    ab = tex.split("\\begin{abstract}")[1].split("\\end{abstract}")[0]
    body = tex.split("\\end{abstract}")[1]
    body = re.sub(r"\\begin\{deluxetable\*\}.*?\\end\{deluxetable\*\}", "", body, flags=re.S)
    body = body.split("\\begin{acknowledgments}")[0]
    body = re.sub(r"\\keywords\{.*?\n", "", body)
    print(f"abstract {words(ab)} words (limit 150); body {words(body)} words (limit 1350)")


if __name__ == "__main__":
    main()
