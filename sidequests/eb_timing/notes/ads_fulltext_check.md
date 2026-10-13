# ADS full-text search: the four southern candidates

Done 2026-10-12 in the ADS web UI (ui.adsabs.harvard.edu) with `full:"…"` queries. **Repeated 2026-10-13 through the ADS API**
(`scripts/ads_check.py` → `notes/ads_api_check.json`, with highlights). The API run adds the Gaia DR3 source IDs, 2MASS,
GSC, CD and CPD names, and more positive controls. Its results are identical to the UI's; the extra identifiers add no hits. The full-text index covers abstracts, bodies and, for many journals,
arXiv full text.

## Positive control, and a coverage caveat

`full:"424461577"` returns **1** paper: Marcadon & Prša 2024. That TIC is also tabulated in Mitnyan+2024 and Borkovits+2026,
so **ADS full text misses identifiers that appear only in tables** (or in journals it does not index in full, e.g.
MDPI *Universe*). A null result here is therefore not proof of absence. It complements the LaTeX-source grep of
27 triple/ETV papers (`scripts/known_triples.py`) and the SIMBAD/MSC/WDS check (`ncvz_unreported_check.md`).

### Controls (API run, 2026-10-13)

| Control | Hits | Notes |
|---|---|---|
| TIC 424461577 (tabulated in Mitnyan+2024, Borkovits+2026, Marcadon & Prša 2024) | 1 | Only Marcadon & Prša |
| TIC 229771234 (same three papers) | 2 | Marcadon & Prša plus a δ Sct catalogue; misses Mitnyan and Borkovits |
| TIC 167692429 (named in title/text of Borkovits+2020) | 12 | Discovery paper and its follow-ups all found |
| TIC 220397947 (same) | 11 | Same |

**Interpretation:** a star that is the *subject* of a paper (named in its title, abstract or text) is reliably found.
A star that is only a *row in a large table* is mostly missed. So the search rules out any paper about our candidates,
but cannot rule out a catalogue row, for which the LaTeX-source grep of 27 papers is the complementary check.

## Results

| Candidate | Query | Hits | What they are |
|---|---|---|---|
| TIC 32000625 (TYC 9155-22-1) | `full:"32000625" OR full:"TYC 9155-22-1" OR full:"TYC 9155-0022-1"` | 4 (2 relevant) | Two **conference abstracts** from CATS 2023 (Comets, Asteroids, exoplanetary Transits and variable Stars, 1–2 July 2023), J. Phys. Stud. 27(3), 3998-2: **Larkin & Pyshna 2023** (2023JPhSt..2739982L) and **Dzygunenko & Tvardovskyi 2023** (2023JPhSt..2739982D). Both classify it as EA, refine P and T0 from TESS minima (MAVKA), and **plot an O−C curve**. Neither abstract mentions a cyclic O−C, a third body or an LTTE orbit. The other 2 hits are seagrass papers that match the number by coincidence. |
| TIC 279323859 (UCAC4 155-008982) | `full:"279323859" OR full:"UCAC4 155-008982"` | 0 | — |
| TIC 278861826 (TYC 8544-1497-1) | `full:"278861826" OR full:"TYC 8544-1497-1"` (plus GSC/Gaia/2MASS in the API run) | 1 | Chakraborty, Wheeler & Kipping 2020 (MNRAS 499, 4011), "weirddetector" periodic signals in TESS Year 1. Table entry only (no text highlight). |
| TIC 294327312 (HD 56910) | `full:"HD 56910" OR full:"HD56910" OR full:"HIP 35095" OR full:"HIP35095" OR full:"294327312"` (plus CD/CPD/Gaia/2MASS) | 1 | Mellon et al. 2019 (ApJS 244, 15), bright southern variables in bRing. Named once in a list of stars whose catalogued period was confirmed by inspecting the folded light curve. No timing work. |

## Other checks in the same pass

- **Web search** for each identifier (standard engine): no ETV or triple reports, and no AAVSO/BAV/IBVS/OEJV timing papers.
- **Czavalinga 2023 PhD dissertation** (Szeged; Borkovits and Mitnyan were consultants; doktori.bibl.u-szeged.hu/11937):
  none of the four IDs appear, and southern-CVZ plans are not mentioned.
- **Szeged diploma thesis record 136154** turned up in a web search for TIC 279323859/278861826. The repository returned
  an EPrints system error at both URL forms, so its content is **unverified**. It may be a fuzzy search-engine match. Recheck
  later, because a Szeged student thesis on southern ETVs would be a direct overlap.
- Vokrouhlický 2026 (arXiv:2603.21799, A&A) is LTTE theory and names no stars.

## Verdict

- **TIC 279323859, 278861826 and HD 56910:** no prior timing or third-body work found.
- **TIC 32000625:** O−C diagrams from TESS were made in 2023 by two student projects (conference abstracts), but no
  third body or orbit was reported. Our LTTE orbit is still new. **Both abstracts must be cited** in any note, e.g.
  "eclipse times were previously analysed for period refinement by Larkin & Pyshna (2023) and Dzygunenko &
  Tvardovskyi (2023)".
- Open item: the Szeged thesis 136154.
