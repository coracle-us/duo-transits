"""
Pre-submission gate for Research Notes of the AAS. Required before any RNAAS upload or resubmission.

AAS's initial quality check counts the abstract (limit 150) and the main text including acknowledgments, figure and
table captions **and references** (limit 1350). Both AAS84041 and AAS84033 bounced on 2026-11-24 because our old
counter left the references out. The 2026-11-25 audit (rnaas/audit-2026-11-25.md on branch rnaas-audit-2026-11-25)
found three more blind spots, closed here:
  1. references kept in a .bib file were not counted (the .bbl is now expanded into the count);
  2. body text between the keywords and the first \\section was counted but not caught, and AASTeX 7 typesets it
     ABOVE the title (AAS84041's page 1 went to AAS that way); the gate now fails on it;
  3. nobody looked at page 1: the gate now renders page 1 of the built PDF to PNG and fails if body text comes
     before the ABSTRACT heading on that page.

Checks (exit status 1 if any fails):
  - abstract <= 150 words;
  - main text with references <= 1300 by our count (a 50-word margin under AAS's 1350);
  - no body text before the first \\section (.tex);
  - with --pdf (default: <stem>.pdf next to the .tex, if present): page 1 rendered to <stem>_page1.png, and at most
    60 words on page 1 before "ABSTRACT" (title, author, affiliation and the draft header fit easily).
It also prints the AAS range: our count x 0.993-1.022, calibrated on the two AAS quality-check counts of
2026-11-24 (AAS84041: AAS 1650 vs ours 1662; AAS84033: AAS 1640 vs ours 1605).

Usage:  python scripts/rnaas_gate.py FILE.tex|FILE.md [--pdf FILE.pdf] [...]
Paste the output into the commit message of the version you upload, and look at the page-1 PNG.

.tex: counted with singles/rnaas/aas_wordcount.py's rules, after expanding \\input{} and \\bibliography{} (from the
.bbl, which must exist: build first); abstract = the abstract environment; main text = everything from
\\end{abstract} to \\end{document} (tables, captions, acknowledgments, facilities/software, AI statement and
references; \\bibitem keys too, which slightly overcounts).
.md drafts: abstract = the "Abstract" section; main text = every later section, references included.
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "singles" / "rnaas"))
from aas_wordcount import words  # noqa: E402

ABSTRACT_MAX, TEXT_MAX, AAS_LIMIT = 150, 1300, 1350
AAS_RATIO = (0.993, 1.022)
PAGE1_MAX_WORDS_BEFORE_ABSTRACT = 60
TOK = re.compile(r"\S*[A-Za-z0-9]\S*")


def expand_inputs(text, base):
    def sub(m):
        f = base / m.group(1)
        f = f if f.suffix else f.with_suffix(".tex")
        return expand_inputs(f.read_text(encoding="utf-8"), f.parent) if f.exists() else m.group(0)
    return re.sub(r"\\input\{([^}]*)\}", sub, text)


def bbl_as_text(bbl):
    """A .bbl reduced to what the reference list prints: no preamble, no \\bibitem[label]{key}, macros kept as words."""
    t = bbl.read_text(encoding="utf-8")
    t = t[t.index("\\bibitem"):] if "\\bibitem" in t else ""
    t = t.replace("\\end{thebibliography}", "")
    out, i = [], 0
    while True:                                   # drop each \bibitem[...]{key}; the [...] label nests braces
        j = t.find("\\bibitem", i)
        out.append(t[i:j if j >= 0 else len(t)])
        if j < 0:
            break
        k = j + len("\\bibitem")
        if t[k:k + 1] == "[":
            depth = 0
            for k in range(k, len(t)):
                depth += {"[": 1, "]": -1}.get(t[k], 0)
                if depth == 0:
                    break
            k += 1
        k = t.index("}", k) + 1                   # the {key}
        i = k
    t = " ".join(out)
    t = re.sub(r"\\bibinfo\{[^}]*\}", " ", t)                         # \bibinfo{field} names; keep the {value}
    t = re.sub(r"\\dodoi\{([^}]*)\}", r" doi: \1 ", t)                   # printed as "doi: 10.xxx"
    t = re.sub(r"\\(doarXiv|doeprint)\{([^}]*)\}", r" \2 ", t)           # printed as a URL
    t = re.sub(r"\\(apj|apjs|apjl|aj|aap|mnras|pasp|nat|araa|icarus|psj)\b", " J ", t)   # journal macros = 1 word
    return t


def expand_bibliography(text, base):
    m = re.search(r"\\bibliography\{([^}]*)\}", text)
    if not m:
        return text, None
    stem = m.group(1).split(",")[0].strip()
    bbl = base / f"{Path(stem).stem}.bbl"
    if not bbl.exists():
        sys.exit(f"rnaas_gate: {bbl} not found; build the note (bibtex) first so the references can be counted")
    return text[:m.start()] + bbl_as_text(bbl) + text[m.end():], bbl.name


def count_tex(path):
    t = expand_inputs(path.read_text(encoding="utf-8"), path.parent)
    t, bbl = expand_bibliography(t, path.parent)
    a0, a1 = t.index("\\begin{abstract}"), t.index("\\end{abstract}")
    end = t.index("\\end{document}")
    s = t.find("\\section", a1)
    pre = t[a1 + len("\\end{abstract}"):s if s >= 0 else end]
    pre = re.sub(r"(?m)(?<!\\)%.*$", "", pre)
    pre = re.sub(r"\\keywords\{.*\}", "", pre)
    pre_words = words(pre)
    return words(t[a0 + len("\\begin{abstract}"):a1]), words(t[a1:end]), pre_words, bbl


def md_words(s):
    s = re.sub(r"<!--.*?-->", " ", s, flags=re.S)
    s = re.sub(r"!\[[^\]]*\]\([^)]*\)", " ", s)              # images
    s = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", s)           # links keep their text
    s = re.sub(r"[#*_`|>]", " ", s)
    return len(re.findall(r"\S*[A-Za-z0-9]\S*", s))


def count_md(path):
    t = path.read_text(encoding="utf-8")
    heads = [(m.start(), m.group(1).strip()) for m in re.finditer(r"(?m)^#{1,6}\s*(.+)$", t)]
    i = next(k for k, (_, h) in enumerate(heads) if re.fullmatch(r"(?i)abstract", h))
    a_start = heads[i][0] + t[heads[i][0]:].index("\n")
    a_end = heads[i + 1][0] if i + 1 < len(heads) else len(t)
    return md_words(t[a_start:a_end]), md_words(t[a_end:]), 0, None


def page1(pdf):
    """Render page 1 to PNG; return (png path, words on page 1 before ABSTRACT) or (None, None) without poppler."""
    if not (shutil.which("pdftoppm") and shutil.which("pdftotext")):
        return None, None
    stem = pdf.with_name(pdf.stem + "_page1")
    subprocess.run(["pdftoppm", "-r", "80", "-png", "-singlefile", "-f", "1", "-l", "1", str(pdf), str(stem)],
                   check=True)
    txt = subprocess.run(["pdftotext", "-f", "1", "-l", "1", str(pdf), "-"], capture_output=True, text=True).stdout
    i = txt.find("ABSTRACT")
    return stem.with_suffix(".png"), len(TOK.findall(txt[:i] if i >= 0 else txt))


def main(argv):
    paths, pdfs, k = [], {}, 0
    while k < len(argv):
        if argv[k] == "--pdf":
            pdfs[len(paths) - 1] = Path(argv[k + 1])
            k += 2
        else:
            paths.append(Path(argv[k]))
            k += 1
    ok = True
    for n, p in enumerate(paths):
        ab, text, pre, bbl = (count_md if p.suffix == ".md" else count_tex)(p)
        lo, hi = round(text * AAS_RATIO[0]), round(text * AAS_RATIO[1])
        problems = []
        if ab > ABSTRACT_MAX:
            problems.append(f"abstract over {ABSTRACT_MAX}")
        if text > TEXT_MAX:
            problems.append(f"main text over {TEXT_MAX}")
        if pre:
            problems.append(f"{pre} words of body text before the first \\section (AASTeX 7 sets them above the title)")
        pdf = pdfs.get(n) or (p.with_suffix(".pdf") if p.suffix == ".tex" and p.with_suffix(".pdf").exists() else None)
        shot = ""
        if pdf:
            png, before = page1(pdf)
            if png is None:
                shot = "; page 1 not rendered (poppler missing)"
            else:
                shot = f"; page 1: {png.name}, {before} words before ABSTRACT"
                if before > PAGE1_MAX_WORDS_BEFORE_ABSTRACT:
                    problems.append(f"page 1 has {before} words before ABSTRACT: body text above the title?")
        elif p.suffix == ".tex":
            shot = "; no PDF found, page 1 not checked"
        good = not problems
        ok &= good
        refs = f" (references from {bbl})" if bbl else ""
        print(f"rnaas_gate {'PASS' if good else 'FAIL'}: {p}: abstract {ab}/{ABSTRACT_MAX}; "
              f"main text with references {text}/{TEXT_MAX}{refs}; AAS range {lo}-{hi} (limit {AAS_LIMIT})" + shot
              + ("" if good else "; FAIL: " + "; ".join(problems)))
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1:]))
