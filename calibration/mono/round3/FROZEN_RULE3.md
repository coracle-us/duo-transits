# Round 3: the rule proposed for the freeze

**Status: proposed. Not frozen until Jon approves step 2 in `../APPROVALS.md`.** Chosen on
already-seen stars only (`TUNING_LOG.md`). In code: `harness3.FROZEN3` (a `rule3.Rule3`).

Everything not listed here is as in round 2's frozen rule (`../round2/FROZEN_RULE.md`): light-curve
choice, binning, detrending, the reference, the event search, the box definitions, the
single-transit exclusion (R1-R3), its probability table (R5) and the event vetting (R4). Symbols
(`dur`, `depth0`, `sig`, box, clean, covered, `thr`, measured) are defined there.

Round 3 changes two things: one more veto on second events (T5) and the alias table (T1).

## Second events (T5)

An event passes if it passes round 2's R4 **and** its `sec_box_snr >= 9`.

`sec_box_snr` of an event: take the detrended data of the event's own sector (t, f). Put boxes of
length `dur` at centres `numpy.arange(t.min() + dur/2, t.max() - dur/2, dur/2)`. Keep the boxes with
coverage `frac >= 0.8` whose centre is more than `dur` from the event time. If fewer than 10 are
left, `sec_box_snr` is undefined and the event fails. Otherwise
`sec_box_snr = event depth / (1.4826 * median(|d - median(d)|))` over the kept box depths d
(`core.box_stats`); if that scatter is zero the event fails. The event depth is the fitted depth
of the event (`depth` in the pipeline's event record).

The alias event is chosen as in round 2 among passing events: the highest SNR among clean events,
else the highest SNR.

## Alias table (T1)

Two transit times: `t_cat` (the catalogued time of the reference, not the pipeline's fit) and
`tc` of the alias event; `t1 < t2` are these two in order, `dT = t2 - t1`. Each time is given an
error of 15 min, so `sig_dt = hypot(15, 15) min = 21.2 min`, whatever the SNR.

For n = 1, 2, ... while `P_n = dT/n >= 1 d` (at most 3000):

- Epochs: `t1 + k P_n` for every integer k from `floor((tmin - t1)/P_n) - 1` to
  `ceil((tmax - t1)/P_n) + 1`, except k = 0 and k = n. `tmin`, `tmax` are the first and last
  time of the detrended data.
- Window of epoch k, in steps of `dur/8`: `h_k = w0 + ceil(2 * far_k * (sig_dt / n) / (dur/8) - 1e-9)`
  with `far_k = min(|k|, |k - n|)` and `w0` round 2's tolerance in steps (4, or 8 on a star with
  more than one planet candidate).
- If `h_k > 16` (the window is wider than +-2 durations) the epoch is **unconstrained**: it neither
  excludes nor supports.
- Otherwise the epoch **excludes** if (a) the box at the predicted time is clean and satisfies
  round 2's box rule (`-3 err <= depth < thr`, reference measured), and (b) no covered box at
  offsets j = -h_k .. h_k steps has `depth >= thr` (its own threshold).
- The epoch **supports** if the box at the predicted time is clean and `|depth - depth0| < 3
  hypot(err, derr)`.
- In `thr` and in the support test the depth error `derr` is the alias event's own depth error,
  as in round 2.

Status of alias n: `excluded` if any epoch excludes; else `supported` if any epoch supports; else
`allowed`. An `allowed` alias with `P_n` below the minimum period from the duration
(`core.pmin_from_duration(dur, mstar, rstar, depth0)`, when finite) becomes `allowed_below_pmin`.

## The call

As round 2: `n_ok` = number of aliases `allowed` or `supported`; `resolved` if 1, `narrowed` if
more, `all excluded` if 0, `no repeat` if no event passes. The truth is in the alias set if
`n = round(dT / P_true)` is in the table, `|dT - n P_true| < dur/2`, and alias n is not `excluded`.

## Other definitions used by the score

- An event is *at the truth* as in round 2 (within `dur/2` of a true transit other than the
  reference). `event_true`: the alias event is at the truth.
- Arm B, own events removed: as round 2 (an event that is not at the truth and lies within `dur`
  of a null-run event is removed before vetting).
- `frac_alias_excluded`: the share of the table's aliases that are `excluded`.
  `frac_alias_excluded_r2`: the same for round 2's table (h_k = w0 for every epoch, no cap) through
  `t_cat` and the same alias event.
- Arm B test trials: the catalogued depth and duration given to the pipeline are the injected
  values times `max(0.2, 1 + 0.1 g)`, g standard normal; one generator per star, numpy seed
  `[20261204, tic]`, two draws per trial in draw order (depth, then duration).

## Not adopted (tried in tuning)

Using the pipeline's fitted reference time; a timing error scaled by duration/SNR; windows capped at
1, 4 or no limit; other growth factors; a veto on events that fit no common period (it removed
real repeats); a higher SNR for events when the reference is weak; changes to the minimum-period
limit (no effect). T2 (refitting the period) and the two-product check of T5 were not tried.
T4 does not enter any criterion and was left as it is.
