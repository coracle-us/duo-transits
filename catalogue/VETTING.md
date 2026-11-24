# Vetting instructions (given to each vetter, word for word)

You are vetting dips in TESS light curves. You work alone. Do not look for, or ask about, the star's
name, any period, any other analysis, or another vetter's verdict. Read only the files in the sheets
folder you are given. Each sheet is a PNG (`<sheet_id>.png`) with a JSON of the same numbers.

Each sheet shows a REFERENCE dip (a catalogued transit-like event on the star) and an EVENT (another
dip on the same star that an automatic search found). Your question is only this:

> Is the EVENT a real dip in the light of this target, of the same depth and duration as the
> REFERENCE, and not an instrumental or sector-wide feature?

You are not asked whether the same body caused both dips.

## The panels

- A1-A4: reference and event, raw and detrended, on one scale. The orange band is one reference
  duration; the dashed line is the reference depth.
- B1-B2: the event's times in other light-curve products of the same sector, if any exist.
- C: the simple-aperture flux (SAP) beside the flux the search used.
- D: the event's whole sector, with flagged cadences and momentum dumps marked, and the distance of
  the event from the sector's edges and gaps.
- E: the same times on neighbouring stars of similar brightness on the same camera and CCD.
- F: a difference image, if one exists. Most sheets say "no difference image".
- G: the numbers.

## Verdicts

Give exactly one verdict per sheet.

- `accept`: the dip is there in the raw and the detrended flux, looks like a transit (flat or rounded
  bottom, ingress and egress, roughly symmetric), matches the reference in depth and duration, and
  nothing below applies.
- `reject`, with one reason code (the main one):
  - `not_in_other_product`: another product covers the event's times with noise low enough to show
    it, and shows no dip.
  - `not_in_sap`: the SAP flux shows no dip, or the dip is made by the systematics correction.
  - `flagged_or_dump`: the dip coincides with flagged cadences or a momentum dump and looks caused by them.
  - `edge_or_gap`: the dip is a ramp or a step at a sector edge or data gap.
  - `on_neighbours`: neighbouring stars show a dip at the same times.
  - `not_transit_shaped`: a ramp, a step, one or two outlying points, a slow wave, a strongly lopsided dip.
  - `detrending_artefact`: there is no dip in the raw flux; the detrending made it.
  - `does_not_match_reference`: the depth or the duration is clearly unlike the reference's.
  - `off_target`: a difference image puts the signal away from the target.
- `cannot_tell`: the sheet does not let you decide. Say in the note what is missing.

A missing panel (no other product, no SAP, no neighbours, no difference image) is not a reason to
reject, and not by itself a reason for `cannot_tell`.
"Event coincides with another catalogued transit on this star: True" is information only. Judge
the dip as any other.
If the REFERENCE is marked "NOT MEASURED", the reference panels may show no visible dip. Judge the
event against the reference depth and duration given in panel G.

## What to write

One CSV, `verdicts_<your label>.csv`, with a header and one line per sheet, every sheet in the folder:

```
sheet_id,verdict,reason,note
```

`reason` is empty unless the verdict is `reject`. `note` is one short sentence, no commas.
Look at every sheet yourself. Do not sample, and do not write a rule to decide from the JSON alone.
