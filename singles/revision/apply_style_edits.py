"""
Editor comments 1-3 on AAS82207: first person singular, no 'duotransit' (except in the precise sense of
Hawthorn et al. 2024), and an explicit description of the owner exclusion. Also the minor light-travel fix
(06:42 -> 06:40 UTC). Each replacement must match exactly once.
Usage:  python singles/revision/apply_style_edits.py
"""
from pathlib import Path

P = Path(__file__).resolve().parents[1] / "rnaas" / "rnaas.tex"

R = [
    ("Mostly Eclipsing Binaries, Two New Duotransits}",
     "Mostly Eclipsing Binaries, Two Candidates with a Second Transit}"),
    ("We re-examined 306 such stars against all TESS data through Sector 108. We used a validated\n"
     "repeat search, our own photometry for sectors not yet processed by the TESS pipelines,",
     "I re-examined 306 such stars against all TESS data through Sector 108. I used a validated\n"
     "repeat search, my own photometry for sectors not yet processed by the TESS pipelines,"),
    ("event we recover, 79 show repeat events.", "event I recover, 79 show repeat events."),
    ("orbits reproduce our periods)", "orbits reproduce the TESS periods)"),
    ("Two on-target duotransits survive.", "Two on-target candidates with two transits each survive."),
    ("2020--2022. We took the 450 open CTOIs flagged as single or duo transits, excluding those of\n"
     "groups still working their candidates. We then removed",
     "2020--2022. Of 794 open CTOIs flagged as single or duo transits, I excluded the 344 uploaded by\n"
     "groups still following up their own single transits, identified by the ExoFOP uploader: members of\n"
     "the TESS Single Transit Planet Candidate working group (223), the NGTS duotransit program\n"
     "\\citep[55;][]{Hawthorn2024}, and \\citet[66]{Salinas2025}. From the remaining 450 I then removed"),
    ("Another 22 fall outside the range we analyzed:", "Another 22 fall outside the range I analyzed:"),
    ("100\\%. We analyzed the remaining 305 stars", "100\\%. I analyzed the remaining 305 stars"),
    ("We used SPOC 2-min, TESS-SPOC", "I used SPOC 2-min, TESS-SPOC"),
    ("without pipeline products, we performed our own aperture photometry",
     "without pipeline products, I performed my own aperture photometry"),
    ("rejected. We validated the method blind, starting from the first transit only, on 26 TESS\n"
     "duotransit hosts with known periods \\citep[including][]{Hawthorn2024}.",
     "rejected. I validated the method blind, starting from the first transit only, on 26 TESS\n"
     "stars with two transits and a known period, including duotransit candidates of \\citet{Hawthorn2024}."),
    ("We localized sources with TESScut", "I localized sources with TESScut"),
    ("do not show the catalogued event in our reductions,", "do not show the catalogued event in my reductions,"),
    ("Three Gaia SB1 orbits reproduce our TESS periods", "Three Gaia SB1 orbits reproduce the TESS periods"),
    ("seen so far only in our own Sector~107 photometry.", "seen so far only in my own Sector~107 photometry."),
    ("(Sectors 7, 33, and 108, the last from our own\nphotometry)", "(Sectors 7, 33, and 108, the last from my own\nphotometry)"),
    ("It appears in TESS-SPOC, QLP, and our own\nphotometry", "It appears in TESS-SPOC, QLP, and my own\nphotometry"),
    ("$^a$Our own TESScut photometry.", "$^a$My own TESScut photometry."),
    ("on-target duotransits with alias lists short enough to target,",
     "on-target candidates with two transits and alias lists short enough to target,"),
    ("We also used ExoFOP", "I also used ExoFOP"),
    ("2026-12-28 06:42\\,UTC", "2026-12-28 06:40\\,UTC"),
]

BIB_ANCHOR = "\\bibitem[Pr{\\v s}a et al.(2022)]{Prsa2022}"
SALINAS = ("\\bibitem[Salinas et al.(2025)]{Salinas2025} Salinas, H., Brahm, R., Olmschenk, G., et al. 2025, "
           "MNRAS, 538, 2031\n")


def main():
    s = P.read_text()
    for old, new in R:
        n = s.count(old)
        assert n == 1, f"{n} matches for: {old[:60]!r}"
        s = s.replace(old, new)
    i = s.index(BIB_ANCHOR)
    j = s.index("\n", i) + 1
    s = s[:j] + SALINAS + s[j:]
    P.write_text(s)
    print(f"applied {len(R)} replacements + Salinas2025 bibitem")


if __name__ == "__main__":
    main()
