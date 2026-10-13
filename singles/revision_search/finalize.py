"""
Final tables and REPORT.txt for the unconstrained search (RNAAS AAS82207 revision, Q5).

Inputs : extra_events_auto.csv, data/per_star_summary_auto.csv, control_false_alarms.csv,
         recovery_check.csv, data/control_events_all.csv (classify.py), manual_vetting.csv (by eye,
         from the plots made by plot_events.py), data/eb_per_star/*.json (EB positive controls),
         data/per_star/*.json
Outputs: extra_events.csv, per_star_summary.csv, eb_controls.csv, REPORT.txt

Manual layer: each automatic (a)-(c) event, and (d) events only where a manual row says
apply_to_d=yes, takes the first matching row of manual_vetting.csv (TIC, time window, minimum
depth). verdict real / possible / systematic; systematic -> category d; final_category, when
given, replaces the automatic letter. Class changes are suggested only from 'real' events.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def manual_row(mv, r):
    for m in mv[mv.TIC == r.TIC].itertuples():
        if not (m.T_min <= r.T_btjd <= m.T_max) or r.depth_fit_ppm < m.min_depth_ppm:
            continue
        if r.auto_category == "d" and m.apply_to_d != "yes":
            continue
        return m
    return None


def eb_controls():
    rows = []
    for f in sorted((HERE / "data" / "eb_per_star").glob("*.json")):
        j = json.loads(f.read_text())
        ev = j["events"]
        P = j["prsa_P"]
        ref = j["event_T"] if j["event_T"] else max(ev, key=lambda e: e["depth"])["T"]
        ph = np.array([((e["T"] - ref) / P + 0.25) % 1 - 0.25 for e in ev])
        dep = np.array([e["depth"] for e in ev])
        nred = np.array([e["n_reductions"] for e in ev])
        used = np.zeros(len(ev), bool)
        clusters = []
        for i in np.argsort(-dep):
            if used[i]:
                continue
            m = (np.abs(((ph - ph[i]) + 0.5) % 1 - 0.5) < 0.04) & ~used
            used |= m
            if m.sum() >= 2:
                clusters.append((round(float(np.median(ph[m] % 1)), 3), int(m.sum()), round(float(np.median(dep[m])) * 1e6),
                                 int((nred[m] >= 2).sum())))
        clusters.sort(key=lambda c: -c[2])
        sec = [c for c in clusters if min(c[0], 1 - c[0]) > 0.05]
        rows.append(dict(TIC=j["TIC"], prsa_P_d=P, n_sectors=j["n_star_sectors"], n_events=len(ev),
                         n_control=len(j["control_events"]),
                         clusters_phase_n_depthppm_n2red=" | ".join(f"{c[0]}:{c[1]}:{c[2]}:{c[3]}" for c in clusters),
                         secondary_found=bool(sec),
                         secondary=f"phase {sec[0][0]}, {sec[0][1]} eclipses, {sec[0][2]} ppm" if sec else ""))
    return pd.DataFrame(rows)


def main():
    ex = pd.read_csv(HERE / "extra_events_auto.csv").rename(columns={"category": "auto_category",
                                                                    "category_reason": "auto_reason"})
    mv = pd.read_csv(HERE / "manual_vetting.csv")
    # dump flag (v1 rule comparison)
    dump = {}
    for f in (HERE / "data" / "per_star").glob("*.json"):
        j = json.loads(f.read_text())
        for e in j.get("events", []):
            dump[(j["TIC"], round(e["T"], 4))] = e.get("dump_in_window", False)
    ex["dump_in_window"] = [dump.get((r.TIC, round(r.T_btjd, 4)), False) for r in ex.itertuples()]
    verdict, cat, note = [], [], []
    for r in ex.itertuples():
        m = manual_row(mv, r) if (r.auto_category != "d" or (mv.TIC == r.TIC).any()) else None
        if m is None:
            verdict.append("auto-d" if r.auto_category == "d" else "unvetted")
            cat.append(r.auto_category)
            note.append("")
            continue
        verdict.append(m.verdict)
        if m.verdict == "systematic":
            cat.append("d")
        else:
            cat.append(m.final_category if isinstance(m.final_category, str) and m.final_category else r.auto_category)
        note.append(m.note)
    ex["manual_verdict"], ex["category"], ex["manual_note"] = verdict, cat, note
    assert not (ex.manual_verdict == "unvetted").any(), ex[ex.manual_verdict == "unvetted"][["TIC", "T_btjd"]]
    cols = ["TIC", "CTOI", "note_cls", "T_btjd", "sector", "product", "depth_ppm", "dur_h", "durs_detected_h",
            "depth_fit_ppm", "dur_fit_h", "snr", "beta", "in_two_reductions", "n_reductions", "n_alt_tested",
            "independent_confirm", "alt_products", "n_other_stars_same_time", "dump_in_window", "ref_depth_ppm",
            "ref_dur_h", "ref_source", "depth_ratio", "dur_ratio", "known_period_d", "phase", "auto_category",
            "auto_reason", "manual_verdict", "category", "manual_note"]
    ex = ex[cols].sort_values(["TIC", "T_btjd"])
    ex.to_csv(HERE / "extra_events.csv", index=False)

    ps = pd.read_csv(HERE / "data" / "per_star_summary_auto.csv").rename(columns={"suggested_class": "auto_suggestion"})
    star_note = {
        147035057: "EB confirmed: 41-57 ppt primary eclipses at phase 0.239 of P = 37.31494 d (S6, S87); catalogued 3.2 ppt event is the secondary",
        188835620: "EB, P = 22.217 d: 58 ppt and 11.8 ppt eclipses (phase 0.733); was repeats-inconsistent",
        75577433: "EB, P = 3.2799 d (near-equal eclipses at phase 0 and 0.496); period new",
        83384234: "periodic, P ~ 4.331 d (EB-like, deeper in QLP); was repeats-inconsistent",
        342004815: "off-target EB, P = 39.311 d with secondary at phase 0.644; period new",
        258799245: "EB confirmed: 38 ppt primary at phase 0.753 of P = 30.133 d; catalogued 2.4 ppt event is the secondary",
        230130556: "EB confirmed: 1.06 ppt secondary at phase 0.575",
        270450687: "EB: 1.5 ppt secondary-like eclipse",
        162031311: "EB-like: a 5 ppt eclipse twice the catalogued depth",
        250492521: "single -> EB/blend candidate: a second, V-shaped 1.3 ppt 10 h event; no common period P = dT/n",
        354032500: "(event not recovered) 6 ppt eclipses with P = 1.6242 d",
        354032501: "(event not recovered) 6 ppt eclipses with P = 1.6242 d",
        38005437: "(event not recovered) 2.2 ppt dips every ~1.06 d: short-period variable/EB",
    }
    rows = []
    for p in ps.itertuples():
        e = ex[ex.TIC == p.TIC]
        n = {c: int(((e.category == c) & (e.manual_verdict == "real")).sum()) for c in "abc"}
        npos = int((e.manual_verdict == "possible").sum())
        sugg = star_note.get(p.TIC)
        if sugg is None:
            sugg = p.note_cls + (" (unchanged; possible extra dip, unconfirmed)" if npos else " (unchanged)")
        rows.append(dict(TIC=p.TIC, CTOI=p.CTOI, note_cls=p.note_cls, n_star_sectors=p.n_star_sectors,
                         event_recovered_by_free_search=p.event_recovered_by_free_search,
                         n_known_repeats=p.n_known_repeats, n_known_repeats_recovered=p.n_known_repeats_recovered,
                         n_extra_detected=p.n_extra, n_extra_two_reductions=int(e.in_two_reductions.sum()),
                         n_a_real=n["a"], n_b_real=n["b"], n_c_real=n["c"], n_possible=npos,
                         n_d=int((e.category == "d").sum()), n_control=p.n_control,
                         n_control_passing=p.n_control_passing, auto_suggestion=p.auto_suggestion,
                         suggested_class=sugg))
    ps2 = pd.DataFrame(rows)
    ps2.to_csv(HERE / "per_star_summary.csv", index=False)

    eb = eb_controls()
    eb.to_csv(HERE / "eb_controls.csv", index=False)

    # ------------------------------------------------ numbers
    fa = pd.read_csv(HERE / "control_false_alarms.csv")
    cev = pd.read_csv(HERE / "data" / "control_events_all.csv")
    rec = pd.read_csv(HERE / "recovery_check.csv")
    S = set(ps2[~ps2.note_cls.isin(["event-not-recovered", "error"])].TIC)
    NR = set(ps2[ps2.note_cls == "event-not-recovered"].TIC)
    e272, c272 = ex[ex.TIC.isin(S)], cev[cev.TIC.isin(S)]
    ss272 = int(fa[fa.TIC.isin(S)].n_star_sectors.sum())
    ss34 = int(fa[fa.TIC.isin(NR)].n_star_sectors.sum())
    dips_f = e272[(e272.n_reductions >= 2) & (e272.n_other_stars_same_time < 2)]
    bumps_f = c272[c272["filter"] == "passes"]
    bins = [8, 10, 15, 20, 50, np.inf]
    lab = ["8-10", "10-15", "15-20", "20-50", ">=50"]
    hd, hb = np.histogram(dips_f.snr, bins)[0], np.histogram(bumps_f.snr, bins)[0]
    real272 = e272[e272.manual_verdict == "real"]
    poss272 = e272[e272.manual_verdict == "possible"]
    single = ps2[(ps2.note_cls == "single")]
    r272 = rec[rec.TIC.isin(S)]
    fail = r272[~r272.event_recovered].free_fit_fail_reason.fillna("no fit").value_counts().to_dict()
    nb = lambda d, c: int(((d.category == c)).sum())  # noqa: E731
    L = []
    w = L.append
    w("UNCONSTRAINED DIP SEARCH - RNAAS AAS82207 revision (editor's Q5)")
    w("=" * 70)
    w("Scripts: search_any.py (search + inverted control), classify.py (automatic rules),")
    w("plot_events.py + manual_vetting.csv (by-eye vetting), period_check.py (ephemerides),")
    w("finalize.py (this report). Run order: search_any.py; search_any.py --ebs; classify.py; finalize.py.")
    w("")
    w("1. WHAT THE NOTE'S MATCHED SEARCH IS FOR, AND WHAT IT MISSES")
    w("-" * 70)
    w("The repeat search asks one question: does the CATALOGUED event recur? Only a recurrence of the")
    w("same event (same object, same eclipse) can be used in P = dT/n; the aliases and the flat-data")
    w("exclusions are all computed for that one event. Fixing depth (0.5-2x) and duration to the")
    w("catalogued values is what makes a detection a repeat of THAT event and keeps the number of")
    w("trials (and false repeats) low. A dip of a different depth is not a repeat: a secondary eclipse")
    w("sits at an unknown phase (eccentric orbits), and a different body has its own period, so")
    w("neither can enter P = dT/n. By construction the matched search therefore misses (i) EB")
    w("eclipses of a different depth (secondaries, or a deeper primary when the catalogued event was")
    w("the secondary), (ii) transits of other bodies, (iii) repeats of a different duration, and")
    w("(iv) repeats where the CTOI depth itself is wrong. The search below covers those.")
    w("")
    w("2. METHOD (rules fixed in the search_any.py / classify.py docstrings before results)")
    w("-" * 70)
    w("Same light curves as the note: SPOC 2-min, TESS-SPOC, QLP and our TESScut photometry (OWN,")
    w("S>=100 without products); best product per sector for the search; every detection re-measured")
    w("in every other product. Free depth; durations 1, 2, 4, 8, 16, 24 h plus the catalogued D.")
    w("Detection: SNR >= 8 including a red-noise factor (measurable), coverage >= 0.8, baseline on")
    w("both sides in the same segment (>= max(0.25 d, D/2)), >= 0.25 d from segment/sector/orbit")
    w("edges, momentum-dump test, median depth >= half the mean depth, depth >= 3x the 1-h scatter.")
    w("Control: the identical search on 2 - flux. Automatic (d) filter: seen on >= 2 other stars")
    w("within 0.5 d (common-mode), not confirmed in another reduction (SNR >= 4, depth 0.5-2x), no")
    w("other reduction available, or a 2-reduction control brightening on the same star at >= the SNR.")
    w("Changes after first results (disclosed): (1) dump rule v1 -> v2 (v1 vetoed 37 of 272 genuine")
    w("catalogued events, up to SNR 76; v2 vetoes only dumps in the event core whose two sides do not")
    w("both show the dip); (2) the depth/duration used to classify an event comes from a fine")
    w("duration refit (0.5-30 h), and the catalogued-event reference from the free fit at the")
    w("catalogued time and duration, because the coarse 6-step grid split single eclipse families;")
    w("(3) (b) uses durations within 1/3-3x. Detection thresholds were never changed.")
    w(f"Extra events with a dump inside the event window (would fail the v1 rule): "
      f"{int(e272.dump_in_window.sum())} of {len(e272)} in the 272; "
      f"{int(real272.dump_in_window.sum())} of {len(real272)} vetted-real.")
    w("")
    w("3. SANITY CHECKS")
    w("-" * 70)
    w(f"Catalogued event re-detected (all 272, not just 20): {int(r272.event_recovered.sum())}/272. Failures "
      f"(free fit at catalogued T, D): {fail}. 'snr' = SNR < 8 here; the note counted an event as recovered")
    w("at SNR >= 5, so these are the weak catalogued events.")
    w(f"Known repeats re-detected: {int(r272.n_known_repeats_recovered.sum())}/{int(r272.n_known_repeats.sum())}.")
    w("Prsa+2022 EB positive controls (secondary found at its own phase, any depth):")
    for e in eb.itertuples():
        w(f"  TIC {e.TIC}: P = {e.prsa_P_d:.4f} d, {e.n_events} events, clusters (phase:n:depth ppm:n 2-red) "
          f"{e.clusters_phase_n_depthppm_n2red}; secondary {'FOUND: ' + e.secondary if e.secondary_found else 'not found'}")
    w(f"  -> secondaries found in {int(eb.secondary_found.sum())}/{len(eb)} EBs (in 396720998 the second family is")
    w("     seen only in QLP S97/98 and is as deep as the primary, so the Prsa period may be 2x off there).")
    w("")
    w("4. FALSE ALARMS (inverted light curves, same stars, same rules)")
    w("-" * 70)
    w(f"272 stars, {ss272} star-sectors. Inverted-control detections: {len(c272)} raw "
      f"({len(c272) / ss272:.3f} per star-sector); after the 2-reduction and common-mode filters: "
      f"{len(bumps_f)} ({len(bumps_f) / ss272:.3f} per star-sector) on {bumps_f.TIC.nunique()} stars.")
    w(f"Real extra dips (not catalogued event / known repeats): {len(e272)} raw ({len(e272) / ss272:.3f} per "
      f"star-sector); same filters: {len(dips_f)} ({len(dips_f) / ss272:.3f}) on {dips_f.TIC.nunique()} stars.")
    w("SNR distribution after those filters (dips vs inverted):")
    for a, b, c in zip(lab, hd, hb):
        w(f"  SNR {a:>6}: dips {b:4d}   inverted {c:4d}")
    w("Below SNR ~20 extra dips are no more common than inverted 'dips' (variability, flares,")
    w("systematics); the excess is at SNR >= 20 and comes from a handful of eclipsing binaries.")
    w(f"The 34 'event not recovered' stars: {ss34} star-sectors; inverted {int(cev[cev.TIC.isin(NR)].shape[0])} raw.")
    w("")
    w("5. EXTRA EVENTS IN THE 272 STARS")
    w("-" * 70)
    w(f"Stars with >= 1 extra detection: {int((ps2[ps2.TIC.isin(S)].n_extra_detected > 0).sum())}; events {len(e272)}.")
    w(f"Automatic: (d) {nb(e272.assign(category=e272.auto_category), 'd')}, (a) {nb(e272.assign(category=e272.auto_category), 'a')}, "
      f"(b) {nb(e272.assign(category=e272.auto_category), 'b')}, (c) {nb(e272.assign(category=e272.auto_category), 'c')}.")
    w(f"  (d) reasons: {e272[e272.auto_category == 'd'].auto_reason.value_counts().to_dict()}")
    w(f"After by-eye vetting of every non-(d) event: real {len(real272)} on {real272.TIC.nunique()} stars "
      f"[(a) {nb(real272, 'a')}, (b) {nb(real272, 'b')}, (c) {nb(real272, 'c')}], possible {len(poss272)} on "
      f"{poss272.TIC.nunique()} stars, systematic {int((e272.manual_verdict == 'systematic').sum())}.")
    st = ps2[ps2.TIC.isin(S)]
    w(f"Stars (of 272) by automatic category (any event): (a) {int((st.TIC.isin(e272[e272.auto_category == 'a'].TIC)).sum())}, "
      f"(b) {int((st.TIC.isin(e272[e272.auto_category == 'b'].TIC)).sum())}, (c) {int((st.TIC.isin(e272[e272.auto_category == 'c'].TIC)).sum())}, "
      f"(d) only {int(((st.n_extra_detected > 0) & ~st.TIC.isin(e272[e272.auto_category != 'd'].TIC)).sum())}.")
    w(f"Stars by vetted-real category: (a) {int((st.n_a_real > 0).sum())}, (b) {int((st.n_b_real > 0).sum())}, "
      f"(c) {int((st.n_c_real > 0).sum())}; any real {int(((st.n_a_real + st.n_b_real + st.n_c_real) > 0).sum())}; "
      f"possible only {int(((st.n_possible > 0) & ((st.n_a_real + st.n_b_real + st.n_c_real) == 0)).sum())}.")
    w("Stars with real extra events, by note class:")
    for cls, g in ps2[ps2.TIC.isin(S) & ((ps2.n_a_real + ps2.n_b_real + ps2.n_c_real) > 0)].groupby("note_cls"):
        w(f"  {cls}: {len(g)} -> " + "; ".join(f"{t}" for t in g.TIC))
    w("")
    w("6. CHANGES TO THE NOTE'S CLASSIFICATION")
    w("-" * 70)
    sr = single[(single.n_a_real + single.n_b_real + single.n_c_real) > 0]
    w(f"'single' stars (193): {int(((single.n_extra_detected > 0)).sum())} have >= 1 extra detection; "
      f"{int(e272[e272.note_cls.eq('single') & e272.auto_category.ne('d')].TIC.nunique())} have one surviving the "
      f"automatic filters; "
      f"after vetting {len(sr)} real ({', '.join(map(str, sr.TIC))}), "
      f"{int((single.n_possible > 0).sum())} possible ({', '.join(map(str, single[single.n_possible > 0].TIC))}).")
    w("  -> EB/blend candidate: 250492521 (second, differently shaped event; no common period).")
    w("  -> no 'single' star gains a period.")
    for t, s in [(r.TIC, r.suggested_class) for r in ps2.itertuples() if r.TIC in S and r.note_cls != "single"
                 and (r.n_a_real + r.n_b_real + r.n_c_real) > 0]:
        w(f"  {t}: {s}")
    w("Periods (period_check.py; data/period_check.csv): 188835620 22.217 d; 75577433 3.2799 d;")
    w("83384234 ~4.331 d; 342004815 39.311 d; 354032500/1 1.6242 d (34-set). 147035057 deep eclipses")
    w("at phase 0.2394 of 37.31494 d.")
    w("")
    w("7. THE 34 'EVENT NOT RECOVERED' STARS (separately)")
    w("-" * 70)
    e34 = ex[ex.TIC.isin(NR)]
    w(f"Catalogued event detected here: {int(rec[rec.TIC.isin(NR)].event_recovered.sum())}/34. Extra events {len(e34)}: "
      f"real {int((e34.manual_verdict == 'real').sum())} on {e34[e34.manual_verdict == 'real'].TIC.nunique()} stars "
      f"(38005437 short-period variable; 354032500/354032501 P = 1.6242 d), possible "
      f"{int((e34.manual_verdict == 'possible').sum())} on {e34[e34.manual_verdict == 'possible'].TIC.nunique()} stars, rest (d).")
    w("")
    w("Files: extra_events.csv, per_star_summary.csv, control_false_alarms.csv, recovery_check.csv,")
    w("eb_controls.csv, manual_vetting.csv; plots in data/plots/ (gitignored).")
    (HERE / "REPORT.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
