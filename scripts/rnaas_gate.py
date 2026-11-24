"""
Pre-submission word-count gate for Research Notes of the AAS. Required before any RNAAS upload or resubmission.

AAS's initial quality check counts the abstract (limit 150) and the main text including acknowledgments, figure and
table captions **and references** (limit 1350). Both AAS84041 and AAS84033 bounced on 2026-11-24 because our old
counter left the references out. This gate fails if the abstract is over 150 words or the main text, references
included, is over 1300 words (a 50-word margin: AAS's count ran 35 words above ours on AAS84033, about 2%).

Usage:  python scripts/rnaas_gate.py FILE.tex|FILE.md [...]      exit status 0 = pass, 1 = fail
Paste the output into the commit message of the version you upload.

.tex: counted with singles/rnaas/aas_wordcount.py's rules, after expanding \\input{}; abstract = the abstract
environment; main text = everything from \\end{abstract} to \\end{document} (tables, captions, acknowledgments,
facilities/software, AI statement and references; \\bibitem keys too, which slightly overcounts).
.md drafts: abstract = the "Abstract" section; main text = every later section, references included.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "singles" / "rnaas"))
from aas_wordcount import words  # noqa: E402

ABSTRACT_MAX, TEXT_MAX, AAS_LIMIT = 150, 1300, 1350


def expand_inputs(text, base):
    def sub(m):
        f = base / m.group(1)
        f = f if f.suffix else f.with_suffix(".tex")
        return expand_inputs(f.read_text(encoding="utf-8"), f.parent) if f.exists() else m.group(0)
    return re.sub(r"\\input\{([^}]*)\}", sub, text)


def count_tex(path):
    t = expand_inputs(path.read_text(encoding="utf-8"), path.parent)
    a0, a1 = t.index("\\begin{abstract}"), t.index("\\end{abstract}")
    return words(t[a0 + len("\\begin{abstract}"):a1]), words(t[a1:t.index("\\end{document}")])


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
    return md_words(t[a_start:a_end]), md_words(t[a_end:])


def main(paths):
    ok = True
    for p in map(Path, paths):
        ab, text = (count_md if p.suffix == ".md" else count_tex)(p)
        good = ab <= ABSTRACT_MAX and text <= TEXT_MAX
        ok &= good
        print(f"rnaas_gate {'PASS' if good else 'FAIL'}: {p}: abstract {ab}/{ABSTRACT_MAX}; "
              f"main text with references {text}/{TEXT_MAX} (AAS limit {AAS_LIMIT})")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1:]))
