# arXiv posting of AAS82207 (accepted 2026-11-14)

`McBee_AAS82207_arXiv.tar.gz` is the source to upload. It holds:
- `rnaas.tex`, byte for byte the accepted text (tag `aas82207R2-submitted`);
- `aastex701.cls`, AASTeX 7.0.1 (2025/05/09) from CTAN, LPPL.

Rebuild it with `build_arxiv.sh /path/to/aastex701.cls OUTDIR`. On 2026-11-14 the package compiled with
tectonic 0.17.0 with no errors and no undefined references, giving 3 pages that match the AAS build. The
only warning is about a comment in the tectonic bundle's own `lineno.sty`. The source is pure ASCII, has
no figures and has no .bib file, so arXiv's pdflatex needs nothing else.

## Form fields

| Field | Value |
|---|---|
| Licence | CC BY 4.0 |
| Primary category | astro-ph.EP |
| Cross-list | astro-ph.SR |
| Title | Revisiting Stale TESS Single-transit Community Candidates: Mostly Eclipsing Binaries, Two Candidates with a Second Transit |
| Authors | Jon McBee |
| Abstract | `abstract_plaintext.txt` (the manuscript abstract with the TeX removed) |
| Comments | Accepted for publication in Research Notes of the AAS. 3 pages, 1 table |
| Journal-ref, DOI | Blank for now. Add both, with no new version needed, once RNAAS publishes. |

## Notes

- arXiv compiles the source again, so the header reads "Draft version" with arXiv's compile date. That
  is AASTeX's default and is common for RNAAS preprints. The text is unchanged.
- Line 1 of `rnaas.tex` is a source comment ("RNAAS submission draft ... Compile on Overleaf"). Anyone
  who downloads the source can see it. It doesn't appear in the PDF.
- A first submission to astro-ph from an account without an institutional address may need an endorsement.
