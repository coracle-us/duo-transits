"""
Word count for eb_depth_evolution/rnaas/rnaas.tex (RNAAS limits: 1,350 words of text and a 150-word abstract;
the figure, caption, acknowledgments, AI statement and references are excluded).
Usage:  python eb_depth_evolution/rnaas/wordcount.py
"""
import re
from pathlib import Path

TEX = Path(__file__).resolve().parent / "rnaas.tex"


def words(s):
    s = re.sub(r"(?<!\\)%.*", "", s)   # strip comments, not escaped \%
    s = re.sub(r"\\cite(t|p|alt)\{[^}]*\}", "CITE", s)
    s = re.sub(r"\$[^$]*\$", "X", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", "", s)
    s = re.sub(r"[{}~]", " ", s)
    return len([w for w in s.split() if re.search(r"[A-Za-z0-9]", w)])


def main():
    tex = TEX.read_text()
    ab = tex.split("\\begin{abstract}")[1].split("\\end{abstract}")[0]
    body = tex.split("\\end{abstract}")[1]
    body = re.sub(r"\\begin\{figure\}.*?\\end\{figure\}", "", body, flags=re.S)
    body = body.split("\\begin{acknowledgments}")[0]
    body = re.sub(r"\\keywords\{.*?\n", "", body)
    print(f"abstract {words(ab)} words (limit 150); body {words(body)} words (limit 1350)")


if __name__ == "__main__":
    main()
