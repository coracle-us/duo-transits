Dear Dr Lintott,

Thank you for your decision on AAS82225, "Orbital Periods for K2 Single-transit and Single-eclipse Events from K2
Revisits and TESS", and for inviting a revision. The revised note does both things you asked for.

**1. A thorough search of the existing literature.**

For each of the 163 stars I searched:
- ADS, by full text and by object, under every identifier (EPIC, TIC, Gaia DR3, 2MASS, and HD, TYC or BD where one
  exists);
- SIMBAD bibliographies;
- 78 VizieR catalogues: variable-star and eclipsing-binary catalogues, K2 rotation and planet-candidate catalogues,
  and the catalogue of every paper found for the stars in the table;
- the Gaia DR3 orbit catalogue;
- the SPOC data-validation reports, where a SPOC threshold-crossing event (TCE) exists.

Each source for the 15 stars in the table was classified twice, independently: as giving a period, a period
estimate, another kind of study, or only a catalogue listing. The two classifications disagreed on one period, and
reading the paper settled it. The rule, fixed before the search, is that a period is called new only if no source
gives one.

The result:
- **8 of the 15 periods are new.** Seven eclipsing binaries have no prior period, and EPIC 211503363 had only a
  period estimate.
- **5 agree with published periods,** which the note now cites.
- **2 conflict with SPOC TCEs:** TOI-5560.01 and EPIC 211489484. The note explains why my data exclude those TCEs and
  does not call these two periods new.

The table is now grouped by these verdicts. The note also states the search's known gaps: an all-catalogue VizieR
query could not run, and one known paper lists a star just outside the search radius.

**2. A blind search for transits.**

The submitted note slid a box with the catalogued depth and duration along the data, and searched only the events
with high expected signal-to-noise. The revised search:
- uses no catalogued depth, duration or epoch for detection, and has no sample cut;
- covers all 164 events on 163 stars, every K2 campaign and every TESS sector;
- takes its detection thresholds from the same search run on inverted light curves.

It recovers 163 of the 164 catalogued events and gives the reported period for 11 of the 15 systems. The note
explains the four misses: two shallow transit sets fall below the single-event threshold, and two eccentric binaries
have their primaries and secondaries grouped into one family. A linked-event test added after the first full run,
and labelled as post hoc, recovers two of the four.

Blind signals on 47 other stars are listed in the repository. They are not claimed, since the control predicts that
many are noise or stellar variability.

**Other changes.** The title now names the blind search. The audit of the catalogued events is replaced by the blind
recovery. Some descriptive detail was cut to stay within the word limit.

The search rules were written before the search ran. Every change made after the first results is dated and
explained in the repository (`k2singles/revision2/PLAN.md`), and every number in the note traces to a script there
(https://github.com/coracle-us/duo-transits).

Sincerely,
Jon McBee
