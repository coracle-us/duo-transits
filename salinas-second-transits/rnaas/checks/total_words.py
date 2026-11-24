"""
The RNAAS 1,500-word total (https://journals.aas.org/research-note-preparation-guidelines/, read 2026-11-24:
title, headers, captions and references count; 150 of the 1,500 are reserved for the abstract). Uses the
word rule of singles/rnaas/aas_wordcount.py. Prints the parts, the total with the table body (as the AAS82207R1
quality check counted it) and without it (as the guidelines page says: table contents do not count).
Usage:  python salinas-second-transits/rnaas/checks/total_words.py
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "singles" / "rnaas"))
from aas_wordcount import words  # noqa: E402

t = (HERE.parent / "rnaas.tex").read_text(encoding="utf-8")
t = re.sub(r"(?<!\\)%.*", "", t)
title = words(re.search(r"\\title\{(.*?)\}\n", t).group(1))
author = words(re.search(r"\\author.*?\n\\affiliation\{(.*?)\}", t, re.S).group(0).replace("\\author[orcid=0009-0000-7524-4470]", ""))
a0, a1 = t.index("\\begin{abstract}"), t.index("\\end{abstract}")
abstract = words(t[a0 + len("\\begin{abstract}"):a1])
b0, b1 = t.index("\\begin{thebibliography}"), t.index("\\end{thebibliography}")
rest = words(t[a1:b0])
refs = words(t[b0:b1])
tab = t[t.index("\\startdata"):t.index("\\enddata")]
table_body = words(tab)
parts = dict(title=title, author_affil=author, abstract=abstract, body_to_refs=rest, references=refs)
tot = title + abstract + rest + refs
print(parts)
print(f"total (title+abstract+body incl. table, captions, ack, software, AI statement + refs): {tot}; "
      f"without the table body ({table_body}): {tot - table_body}; with author/affiliation: {tot + author}; limit 1500")
