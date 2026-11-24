# Strategy: a living all-sky catalogue of TESS long-period alias sets

*Draft for Jon's approval, 2026-11-22. Nothing here has been started. Statuses are taken from
`docs/PROJECT_FACTORY_BRIEF.md` (snapshot 2026-11-16) and `origin/main` on 2026-11-22.*

## 1. Intent

Stop publishing one short note per catalogue. Build one uniform, calibrated, regularly updated
table of allowed periods for every published TESS single-transit and duotransit candidate, and
describe it in one peer-reviewed paper.

**Why.**
- Catalogue-picking has a low ceiling: one to four objects per catalogue, three scoping studies
  ended NO-GO, and ExoFOP community uploaders reach the easy cases first.
- The project's reusable assets are its method and its verification protocol, not any one target
  list. Five workstreams already ran nearly the same analysis under slightly different rules.
- Jon cannot get an ExoFOP account. A citable, versioned table is a distribution channel that needs
  no account.
- Jon's sign-off is the bottleneck (one unsubmitted note, one pending arXiv upload, three unsent
  emails). Fewer, larger outputs spend that time better.

## 2. Goals, non-goals, anti-goals

**Goals**
1. One frozen method applied to every catalogue (§5.1), giving a per-candidate table: allowed
   periods, a probability for each, the data used, vetting flags.
2. Each exclusion carries a measured error rate from injections and from confirmed planets.
3. A versioned public release (Zenodo DOI), regenerated when new sectors arrive.
4. A peer-reviewed paper (AJ, or PASP/RASTI if the scope ends up methods-heavy).
5. A Gaia DR4 layer: stellar companions flagged, and periods taken from Gaia orbits where they exist.
6. Next-transit windows that discriminate between surviving aliases, published for ground observers.

**Non-goals**
- Discovering new candidates (the blind search in §7 is a later, separate decision).
- Confirming planets. The product is ephemerides and flags.
- Replacing ExoFOP or MonoTools.

**Anti-goals** (failure even if it runs)
- A table whose exclusions are wrong for real planets at a rate nobody measured.
- Claiming someone else's candidate or period as ours. Every row cites its source catalogue and any
  existing community period.
- A "living" product that is stale after two sectors because regeneration needs manual work.

## 3. Constraints

- The hard rules in `CLAUDE.md` §5 and the factory brief §3 and §3a apply unchanged: verification
  over volume, novelty check before any claim, pre-registration, traceability, read-only toward the
  outside world. Only Jon submits, emails or posts.
- Several catalogues belong to groups still working them (TSTPC, Hawthorn/Warwick, Salinas). The
  CTOI scope excluded them for that reason. Including them here needs their knowledge, and ideally
  their agreement (§5.6).
- Data limits already found: TESS-SPOC stops at S72; QLP S100+ timestamps are off by up to about
  6 minutes; own TESScut photometry is diluted by 10-40% and is good for times, not depths.
- Shared checkout: one worktree per branch.

## 4. Approaches considered

| | Approach | Trade-off |
|---|---|---|
| A | Keep going: one RNAAS note per catalogue | Cheap and proven (AAS82207 accepted). Yield per catalogue is small and falling; not peer reviewed; each note re-litigates the method. |
| **B** | **Unified calibrated catalogue + peer-reviewed paper (chosen)** | Reuses everything already run. Needs peer review, a homogeneous rerun and probably a collaborator. Risk: someone already maintains one (gate G0). |
| C | Own blind search of the continuous viewing zones, with completeness, for an occurrence rate | Highest scientific payoff. Competes directly with funded groups and needs the calibration from B first. Kept as the stretch bet (§7). |
| D | Leave exoplanets (Be stars, eclipsing-binary depth evolution, ML vetting) | Each is a fresh start that discards the validated pipeline. Be stars and the EB object continue as side lines, not the flagship. |

Jon chose B on 2026-11-22.

## 5. Design

### 5.1 Unify
Target list: the union of Hawthorn et al. 2024 (85), the stale CTOIs (368 on 359 stars), Salinas
et al. (173 rows, 169 stars), LaCourse & Jacobs K2 events (164), the open single-transit TOIs
(about 156), plus the TSTPC and Planet Hunters TESS lists. Order of a thousand stars; the exact
count is the first deliverable.

One schema, one config, one code path, built from the existing shared modules (`scripts/alias_pipeline`,
`salinas-second-transits/st2core.py`, `tesscut_phot.py`, `k2singles/scripts/k2lib.py`). Earlier
workstream results are not copied in; they are rerun and compared, and differences are explained.

### 5.2 Calibrate (the peer-reviewable core)
What exists: `salinas-second-transits/results/round2/EPOCH_RULE_TEST_2026-11-01.md`. On 83 confirmed
planets, **this project's own former rule** excluded the true period for 46 (55%) when the reference
depth is measured as the pipeline measures it, and 68 (82%) against the TOI table's fitted depth.
The half-depth rule cuts that to 3 (4%), at the cost of many more surviving false aliases at low
SNR (median 102 below epoch SNR 6).

What that does and does not show: it is a result about our rule. Whether published alias lists
(Hawthorn's, MonoTools output) share the weakness is untested, and is the first thing to test,
because it decides whether the paper has a methods finding of general interest or only a
description of a careful pipeline.

To build: replace the in-or-out verdict with a probability per alias, from injection-recovery as a
function of epoch SNR, coverage, cadence and reduction. Validate on confirmed planets held out from
tuning, and on the community-solved candidates.

### 5.3 Release
A table with a fixed schema and a `numbers.md`-style provenance file, on Zenodo with a DOI per
version. Regeneration is one command, run when a sector lands on MAST. The release is public and
permanent, so each one is Jon's action.

### 5.4 Gaia DR4 (release 2026-12-02)
Pre-register before release day, then cross-match every host:
- non-single-star solutions, RV amplitude and RUWE, to flag stellar companions at scale;
- orbital periods from Gaia where a solution exists, compared with the alias set.

Not yet checked: what DR4 publishes per star (epoch astrometry, epoch RVs, magnitude limits). The
pre-registration is written against the DR4 documentation, not against this paragraph.
`singles/TODO_gaia_dr4.md` (two stars) becomes a special case of this sweep.

### 5.5 Ground windows
For each candidate with two to about five surviving aliases and a depth a small telescope can
reach, compute the upcoming windows that separate them and publish the list with the release.
TIC 159614322 (2026-12-28 06:42 UTC) is the first case. How to reach observers (AAVSO exoplanet
section, Exoplanet Watch, ExoClock, Planet Hunters TESS) is Jon's call per channel.

### 5.6 Collaboration
An affiliated co-author solves ExoFOP access and helps at review. Natural openings, all drafts for
Jon to send or not:
- the Salinas group: the four periodic stars and the epoch-rule test;
- the MIT QLP team: the S100+ clock offset;
- Warwick: already contacted about TIC 147035057.

## 6. Phases and gates

| Phase | When | What | Gate to continue |
|---|---|---|---|
| G0 | 1 day, first | Novelty check: does a uniform, updating alias catalogue already exist (TSTPC, MonoTools-based, Planet Hunters, NGTS)? Read the DR4 release documentation. | If one exists and is maintained: narrow to the calibration paper plus the DR4 layer, or contribute to theirs. |
| 1 | by 2026-12-01 | DR4 pre-registration; host list frozen; cross-match code tested on DR3. | Jon approves the pre-registration. |
| 2 | 2026-12-02 onward | DR4 sweep. | Reported whatever it finds. |
| 3 | Dec-Jan | Test whether published alias lists share the old-rule weakness (§5.2). Build the injection calibration. | Decides the paper's framing. |
| 4 | Jan-Feb | Unified rerun; comparison with the five earlier workstreams. | Agreement, or every difference explained. |
| 5 | Mar | Zenodo v1; paper draft; windows list. | Jon releases and submits. |
| 6 | ongoing | Regenerate per sector. | Automatic; Jon releases. |

**Kill or shrink criteria**
- G0 finds the niche occupied: shrink as above.
- Phase 3 shows the calibrated probabilities change few published alias sets: publish the DR4 layer
  and the unified table as a data paper, and drop the methods claim.
- No collaborator and a first referee report that turns on affiliation or access: fall back to
  Zenodo plus an RNAAS pointer.

## 6a. G0 result (2026-11-22)

Six web searches and three documents read. Not exhaustive: no ADS full-text search, and the second
half of the preprint below was read through a summary, not verbatim.

**Verdict: proceed, with the scope changes listed at the end.** No maintained, all-sky, updating
alias catalogue was found. The niche is partly occupied.

| Who | What | Overlap with this plan |
|---|---|---|
| Pagés Navarrete, arXiv:2611.03019 (submitted to MNRAS 2026-11-02, v2 11-17; independent researcher) | 187 TOIs with P = 0 or P ≥ 200 d (94 duo, 93 mono), all TESS data through S107. Per-alias weights, 300 duo-transit injections, difference images, Gaia DR3 flags. 28 periods resolved, 8 fixed by exclusion. | **High for TOIs.** Calibrated alias weights for TOI duos are no longer new. Not covered: non-TOI candidates, mono-transit injections ("left for future work"), a test on confirmed planets, DR4, a maintained release. |
| Skinner, Dragomir, Essack & TSTPC (2026 workshop poster, Univ. of New Mexico) | "PeriodSearch" pipeline: re-searches all sectors for their own candidates and runs MonoTools alias probabilities. 62 candidates so far (11 mono, 11 duo, 39 multi); 5 unsolved targets sent to ExoFOP. | Their own list only, in progress, no paper found. Confirms TSTPC candidates should stay out unless they agree. |
| T16 Planet Hunt, arXiv:2604.18579 | 411 single-transit candidates from Cycle 1, published without period constraints. | **New input**, larger than every list in §5.1. Not yet checked for overlap with the others or for later work on it. |
| Independent GitHub efforts (FBLS, two dated-record repos) | Blind searches and per-candidate records, not peer reviewed. | Low. Shows the hobbyist space is active. |

**Gaia DR4, from ESA's content page (marked provisional, last updated 2026-09-14)**
- Archive reopens with DR4 on 2026-12-02 at 12:00 CET. Data span 2014-07-25 to 2020-01-20.
- Epoch astrometry for all sources in the main astrometric processing; per-transit radial
  velocities (`rvs_epoch_parameters_single`, 6.9 billion rows); `nss_two_body_orbit` with 5.9 million
  orbits (DR3 had far fewer); a `gaia_exoplanet_candidates` table. Magnitude limits and period ranges
  are not stated.
- **The Gaia ESA Archive goes offline on 2026-11-25 at 12:00 CET.** Any DR3 test of the cross-match
  code against that archive has to run before then, or use a mirror (VizieR, MAST).

**Scope changes that follow**
1. The paper's claim moves from "calibrated alias probabilities" to: the candidates nobody has
   treated uniformly (CTOIs, Hawthorn, Salinas, K2, T16 singles), mono-transit calibration, the
   confirmed-planet test, the DR4 layer, and a maintained release.
2. The 187 TOIs in arXiv:2611.03019 are a comparison set, not a claim. Rerun them, compare with
   the published table (`table2_full.csv`, arXiv ancillary), cite, and report disagreements.
3. Add T16's 411 singles to §5.1 after an overlap and prior-work check.
4. Phase 3's first question widens: does the old-rule weakness appear in MonoTools-style
   exclusion and in the 2611.03019 rule, not only in ours?
5. Pagés Navarrete is a possible collaborator as well as a competitor. Whether to write is Jon's call.

## 6b. First phase 3 result (2026-11-22)

`calibration/README.md`. Hawthorn et al.'s published alias limits kept all 26 periods that were later
solved (22 certainly in the list; 4 need the full online list). The exclusion weakness found on
2026-11-01 was therefore this project's own, not a property of published catalogues.

**Consequence, by the shrink criterion in section 6:** the paper carries no claim that published alias
lists are unreliable. What remains for the methods part is the single-transit calibration nobody has
done (`calibration/mono/prereg_mono_calibration.yaml`, a draft), and an empirical check of the period
prior: on the 26 solved candidates a flat prior is rejected and P^-5/3 fits best, with a selection
caveat.

## 7. Stretch bet (not approved by this document)

A blind single-transit search of the continuous viewing zones with measured completeness, aimed at
an occurrence rate for giants at 100-1000 d. It depends on the §5.2 calibration and competes with
funded groups. Decide after phase 4.

## 8. Existing work under this strategy

| Item | Proposal |
|---|---|
| AAS82207 (accepted) | arXiv upload as planned; becomes a cited input. |
| AAS82225 revision (K2) | Finish it; the blind search built for it feeds §5.2. |
| EB note (EPIC 202126877) | Submit as is. Independent of the flagship. |
| Salinas second transits | No separate note. Results go to the Salinas group (§5.6) and into the catalogue. |
| Boyajian's star | Publish the Kepler survivors with the window-3 null, or drop. No new TESS data before Aug 2027. |
| Be-star pulsations (PR #8) | Off-thesis. Own repo, own decision. |
| Robovetter loop, new catalogue scoping studies | Parked / stopped. |

## 9. Dates that do not wait

- **2026-12-02** Gaia DR4. The pre-registration must be approved before it.
- **2026-12-28** TIC 159614322 transit. The Planet Hunters TESS email is drafted and unsent.
- 2026-12-08, 12-30, 2027-01-21, 02-11: EPIC 202126877 predicted eclipses.

## 10. Calls made in this draft (Claude's, for Jon to correct)

- Journal: AJ first.
- Earlier results are rerun, not imported.
- Salinas results go by email, not as a note.
- Be stars leave this repo.
- The stretch bet waits for phase 4.
- The catalogue includes candidates owned by active groups, subject to §5.6.

## 11. Open for Jon

1. Is including the owner groups' candidates acceptable before they have answered?
2. Is a public Zenodo release under your name and Coracle the outlet you want?
3. Boyajian: publish the null or drop?
