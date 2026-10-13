"""
RNAAS style for single-author notes is first person singular (editor comment 1 on AAS82207). This converts
'we/We/our/Our/ours' to 'I/I/my/My/mine' as whole words in a LaTeX file, skipping comment lines. 'us' is left
alone; it only appears in URLs. It prints each changed line for review.
Usage:  python singles/revision/first_person_singular.py FILE.tex [FILE.tex ...]
"""
import re
import sys
from pathlib import Path

MAP = [(r"\bWe\b", "I"), (r"\bwe\b", "I"), (r"\bOur\b", "My"), (r"\bour\b", "my"), (r"\bours\b", "mine")]


def convert(path):
    lines = Path(path).read_text().split("\n")
    n = 0
    for k, line in enumerate(lines):
        if line.lstrip().startswith("%"):
            continue
        new = line
        for pat, rep in MAP:
            new = re.sub(pat, rep, new)
        if new != line:
            n += 1
            print(f"{path}:{k + 1}: {new.strip()[:150]}")
            lines[k] = new
    Path(path).write_text("\n".join(lines))
    print(f"{path}: {n} lines changed")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        convert(p)
