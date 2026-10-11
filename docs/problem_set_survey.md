# Next problem set: survey and plan

*Started 2026-10-11. Living document: update it as options are tested or ruled out.*

## 1. Why we're looking

The duotransit project (Hawthorn et al. 2024 catalogue, `CLAUDE.md` §1–8) reached a natural
end. The catalogue is heavily picked over:

| Stage | Count | What it means |
|---|---|---|
| Already solved by others when we started | 26 of 85 | Papers, an NGTS PhD thesis, and ExoFOP community users already mining TESS Sectors 87–103 |
| Leads where the new data gave a single period | 9 | — |
| Signal comes from a neighbouring star | 4 of 9 | Shown with difference imaging |
| Source ambiguous | 1 of 9 | — |
| On-target periods that passed every check | 3 | 207783865 (24.839 d), 38138512 (30.836 d), 153838604 (56.079 d) |

The results went to the catalogue authors by email
(`results/outreach/email_to_warwick_duotransit_team.md`); the vetting record is
`results/vetting/vetting_log.csv`.

**Constraint learned:** Jon cannot get an ExoFOP account, because accounts require an
institutional affiliation. Results therefore leave the project in one of two ways:
- an RNAAS note, which needs no affiliation;
- email to the groups that own the candidates.

## 2. Assets we can reuse

| Asset | Script | Validation |
|---|---|---|
| Alias pipeline | `scripts/alias_pipeline.py` | 26/26 known periods keep the true alias; 13 recovered uniquely from TESS alone |
| Light-curve vetting | `scripts/vet_candidate.py` | Joint fit, odd/even, secondary eclipse, SAP/PDCSAP/QLP comparison |
| Difference-image localisation | `scripts/diff_image.py` | Confirmed planets TOI-5542 b and TOI-2589 b locate on target within ~1–3″ |
| Independent period check | `scripts/independent_bls.py` | Shares no code with the pipeline |
| Literature and ExoFOP triage | `fetch_exofop.py`, `triage_*.py`, `build_status.py` | — |
| Data loader | `scripts/fetch_data.py` | Fetches SPOC, TESS-SPOC and QLP light curves |

**Data lessons:**
- QLP raw flux in the long Sectors 97/98, and in Sector 103, has jumps at the start of each
  data segment. The pipeline now handles this with per-segment baselines plus trimming of
  QLP segment starts.
- TESS-SPOC full-frame light curves stop at Sector 72, so recent sectors rely on QLP or our
  own TESScut photometry.

## 3. Survey (background agent, 2026-10-11; spot-checked)

The survey's working files are in the session scratchpad. The numbers below were re-derived
where noted.

| Rank | Problem set | Size | How picked over | Our edge | Effort |
|---|---|---|---|---|---|
| 1 | Stale single-transit CTOIs | ~450 | Low: 36/45 sampled have no period upload | Repeat-transit search + aliases + difference imaging | Pilot 2–3 weeks; full ~2 months |
| 2 | Source localisation of single-event candidates | ~550–800 | Low: no systematic catalogue found | Validated difference-image tool | 4–6 weeks |
| 3 | Open long-period / duotransit TOIs | 163 | High: 53 already have a period on ExoFOP | Little; use as a validation set | — |
| 4 | Our own northern duotransit catalogue | Probably 50–100 | No catalogue exists, but PHT, TSTPC and Tuson pipelines already catch these | Same pipeline at scale | Months; heavy downloads |
| 5 | K2 single transits re-observed by TESS | 164 events | Not done systematically | Alias pipeline | 3–4 weeks; low yield |
| — | Kepler long-period singles (Herman et al.) | 15 | — | **Dropped:** too few TESS windows | — |

**Notes on the size column:**
- **Row 1** covers never-promoted single/duo CTOIs from groups no longer actively working them,
  mostly uploaded 2020–22. Re-derived count: 450 CTOIs on 435 TICs. The agent's filter gave 472.
- **Row 2** covers 76 open P=0 TOIs plus rows 1 and 3.
- **Row 3** covers 87 TOIs with P > 200 d and 76 with P = 0.

### Evidence and sources

1. **Stale single-transit CTOIs** (from local `data/raw/exofop_ctoi.csv`).
   - *Who uploaded them:* the Single Transit Finder group (steuer, vogl: ~215), Planet Hunters
     TESS (eisner: ~50), ExoMAST (amarino: 52), and others (lipponen, deen, gsun).
   - *Target quality:* median Tmag ≈ 10.8 and depth ≈ 3,200 ppm, so the signal-to-noise is good.
   - *Coverage:* TESScut lists sectors up to 108, and sampled targets have 1–3 sectors after
     their CTOI date.
   - *Staleness:* some entries are out of date. The agent reports that TIC 1042868, still an
     unpromoted CTOI, is now the confirmed planet HD 88986 b (not yet verified by us).
2. **Localisation gap.** SPOC runs its Data Validation step, including the difference-image
   centroid test, only on threshold crossings with two or more transits (per arXiv:2103.12538
   and arXiv:1909.09094 as summarised by the agent). Single events get no pixel-level vetting
   from the pipeline. TESS-Localize (Higgins & Bell, arXiv:2204.06020) is applied case by case.
   On Hawthorn, 4 of our 9 leads were off-target, which suggests the off-target fraction is large.
3. **Open TOIs are picked over:**
   - ExoFOP users (mason, slee, lhill, vincent) post solutions from new sectors continually.
   - Bass & Fabrycky checked TOIs with P > 27 d through S62.
   - Tschudi (arXiv:2607.23781) resolved four single/dual-transit M-dwarf TOIs.
   - RAVEN (arXiv:2603.22597) posts periods.
4. **No northern duotransit catalogue.** arXiv full-text search for "duotransit" returns only
   arXiv:2310.17268 and arXiv:2402.07800, and the review arXiv:2604.09254 lists none.
5. **K2 has an amateur RNAAS precedent:** LaCourse & Jacobs (arXiv:1802.05823). One K2 event
   plus one TESS event about 2,500 d apart still leaves hundreds of possible aliases.

**Do not duplicate:** the TESS Single Transit Planet Candidate (TSTPC) group (Skinner, Dragomir,
Essack) runs a MonoTools all-sector period search on its own ~220 CTOIs (2026 poster:
https://vmnexsci.ipac.caltech.edu/workshop/2026/posters/Poster_BrettSkinner_33.pdf). Exclude
their candidates, or contact them about collaborating.

## 4. Current recommendation: combine options 1 and 2

> "Revisiting ~450 stale TESS single-transit community candidates with four more years of
> TESS data."

Every candidate gets one of three outcomes:
- a period or narrowed alias set;
- an off-target flag from difference imaging;
- "still single, no repeat transit in N new sectors".

**Expected yield (be realistic):** with 1–3 new sectors per target, only a minority will show a
second transit. Most of the output will be vetting (off-target and false-positive flags) rather
than periods. Null rows are still informative and publishable in an RNAAS table.

**Validation set:** the community-solved open TOIs (row 3) give dozens of known answers.

## 5. Pilot plan (not started)

1. **`scripts/build_ctoi_singles.py`:** filter `exofop_ctoi.csv` for candidates that are:
   - not promoted, and not FP/FA/CP/KP;
   - single or duo (P = 0/NaN, or notes saying single/duo);
   - not uploaded by essack, dragomir, skinner, hawthorn or salinas.

   Then cross-match with `nasa_pscomppars.csv` and the TOI table, fetch the ExoFOP JSON per
   target and drop any with community period uploads, and keep targets with at least one
   TESS sector after the event.
2. **Single-epoch mode for the alias pipeline:** with one known event there is no dT, so search
   for repeat events across all sectors with a duration-matched transit template. Any
   corroborated repeat gives the dT; then compute aliases as before.
3. **Difference imaging at every event,** including the original single transit.
4. **Pilot on ~20 bright targets.** Measure the yield (repeat transits found, off-target
   fraction) before committing to the full set.
5. **Decision point:** if the pilot finds less than about 10% repeats and the off-target
   fraction is uninformative, revisit options 2 or 5.

## 6. Open questions

- Is HD 88986 b really the TIC 1042868 CTOI? How many other stale CTOIs are already solved
  under a different name?
- Should we contact the TSTPC group before starting, to avoid overlap and possibly get ExoFOP
  access through them?
- What is the detection threshold for "no repeat transit" in single-epoch mode, i.e.
  injection-recovery on real light curves (the flare side quest's approach)?
- RNAAS policy on AI assistance: check the current AAS journals policy before any submission
  (`CLAUDE.md` §6).
