# Approvals: the catalogue run (`prereg_catalogue.yaml`)

| Version | Approved by | Date | SHA-256 of `prereg_catalogue.yaml` | Note |
|---|---|---|---|---|
| 1.0 | Jon McBee | 2026-11-23 | `1177568b08db3431af02d488cab5bcd7fbee8db7306ebd6d80dab73adfb75d72` | Step 1 of 2. "Approved", in the Claude Code session, after draft 0.1 (PR 28) was put to him with four calls to check (K2 rows left out; two blind agents do the first vetting; Hawthorn tables through the two catalogued times; rows with a stated period kept as candidates), the four validation arms and the known limits. The text is draft 0.1 with only the version and status lines changed. The approval also grants two downloads: MAST light curves of the scored stars, and the arXiv source of 2611.03019. Nothing had been run and no candidate scored. The freeze is a separate approval. |
| 1.0 | Jon McBee | 2026-11-23 | `1177568b08db3431af02d488cab5bcd7fbee8db7306ebd6d80dab73adfb75d72` | Step 2 of 2 (the freeze). "Agreed", in the Claude Code session, to the freeze as proposed in `FREEZE_PROPOSED.md` (commit 0fdfb9f) with its 13 choices, plus a 14th put to him the same day: results stay out of git until he releases them, and only their SHA-256 hashes are committed. Frozen files and hashes: `FREEZE.json` (SHA-256 `f2caf9e2deaaad551d10b5266fe6dac16d6dd694a9c007031bc4357b0746b58c`). No candidate's light curve had been read. He also agreed: PR 28 is squash-merged and the run is made from `main`; no force-push. |

## Corrections that do not change a frozen file

- `prereg_catalogue.yaml`, line 3, says the draft was "commit 3 of PR 28". It was the first commit of
  PR 28, 7f3944c. The file is frozen by hash and is left as it is (Jon, 2026-11-23).
- Its line 12 still reads "step 2 (the freeze) not yet approved", for the same reason. Step 2 is the
  row above.
- Commit 1735b23 of PR 28 put the LaTeX source of arXiv:2611.03019 in the repository; 60d0753 removed
  it. PR 28 is squash-merged so that `main` never holds it. The old commit remains reachable through
  the pull request's own reference on GitHub.


## Deviations approved after the freeze

- 2026-11-23, Jon ("Agreed, split the vetters into batches"): each of the two vetters is several fresh
  agents, one per batch of sheets, with the same instructions (`VETTING.md`); no agent sees a verdict
  of the other vetter. The pre-registration says "two fresh agents". Reason: about 1,900 sheets do not
  fit in one agent's context.

## Changes of procedure after the freeze that change no frozen file (not put to Jon for approval; reported to him)

- 2026-11-23/24: `sheets_runner.py` runs the frozen sheet generators star by star, each in its own
  process with a time limit, retries and checkpoints, because `make_sheets.py` hung on MAST.
- 2026-11-24, about 11:35 CST: MAST's filtered TIC service (`Catalogs.query_criteria`, used by
  `sheets.neighbours`) stopped answering. From then the runner supplies the same list of stars from
  MAST's TIC cone service, filtered by Tmag locally. Sheets made before and after differ only if the
  two services differ at an exact magnitude boundary.
- The vetter-check sheets were vetted in two parts (616, then the rest), because two stars with
  dozens of events each took hours longer.
- An independent re-derivation of the search stage (`catalogue/independent/`) was run before
  vetting finished, not after the tables. It read no sheet and no verdict.
- 2026-11-26: one vetter-2 batch of the vetter check (`b003`) had 39 verdict lines for 40 sheets.
  The missing sheet (8303517e81d7, a false event) was vetted by one more fresh agent, which saw
  only that sheet and `VETTING.md`; its line was appended. Both vetters reject it.
- 2026-11-26: vetter check (V1) scored by `merge_verdicts.py` on all 682 sheets: true events
  accepted 332 of 341 (vetter 1) and 333 of 341 (vetter 2). Criterion (at least 90% each) met, so
  the candidate sheets (1,199, in 30 batches per vetter) go to vetting.
