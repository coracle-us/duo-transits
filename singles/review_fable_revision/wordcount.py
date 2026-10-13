"""Word count and style checks for singles/rnaas/rnaas.tex (review only)."""
import re
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "singles/rnaas/rnaas.tex"
tex = open(path, encoding="utf8").read()


def strip(s):
    s = re.sub(r"%.*", "", s)
    s = re.sub(r"\\(citep|citet|cite|uat|url|label|ref|dataset)(\[[^\]]*\])*\{[^}]*\}", " CITE ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", s)
    s = re.sub(r"[{}$\\~]", " ", s)
    return s


def between(a, b):
    m = re.search(re.escape(a) + r"(.*?)" + re.escape(b), tex, re.S)
    return m.group(1) if m else ""


abs_ = between("\\begin{abstract}", "\\end{abstract}")
body = between("\\end{abstract}", "\\begin{deluxetable*}")
concl = between("\\end{deluxetable*}", "\\begin{acknowledgments}")
table = between("\\begin{deluxetable*}", "\\end{deluxetable*}")
ack = between("\\begin{acknowledgments}", "\\end{acknowledgments}")
ai = between("AI disclosure statement}", "\\begin{thebibliography}")
title = re.search(r"\\title\{(.*?)\}", tex).group(1)
cnt = lambda s: len(strip(s).split())  # noqa: E731
print("title", cnt(title), "| abstract", cnt(abs_), "| sections before table", cnt(body), "| conclusion", cnt(concl),
      "| body total", cnt(body) + cnt(concl), "| table", cnt(table), "| ack", cnt(ack), "| AI", cnt(ai))
main = re.sub(r"\\begin\{thebibliography\}.*", "", tex, flags=re.S)
print("we/our/us:", re.findall(r"\b([Ww]e|[Oo]ur|[Uu]s)\b", main))
print("duotransit lines:")
for l in tex.splitlines():
    if "uotransit" in l:
        print("   ", l.strip()[:170])
print("I/my count", len(re.findall(r"\b(I|my|My)\b", main)))
print("author/his:", re.findall(r"\b([Tt]he author|his)\b", main))
