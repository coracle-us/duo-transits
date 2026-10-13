"""
Interpret the unconstrained dip search (search_any.py) for the RNAAS AAS82207 revision (Q5).

Inputs : singles/revision_search/data/per_star/<TIC>.json (306 stars)
Outputs: singles/revision_search/extra_events.csv, per_star_summary.csv, control_false_alarms.csv,
         recovery_check.csv, and the numbers used by REPORT.txt (printed; report_numbers.json)

RULES (fixed before the per-event results were looked at):
  known       an event is the catalogued event or a known repeat if its window overlaps
              [T_known +/- D_cat/2] with 0.1 d slack; known repeats = the note's repeat-search
              detections (singles/results/<run>/<TIC>.json)
  extra       any other detected event
  reference   our own measurement of the catalogued event (free fit at the catalogued time and
              duration if SNR >= 5, else the merged detection, else the CTOI values)
  (d) systematic / unverifiable, first matching reason:
      common-mode      >= 2 other stars have an extra dip within 0.5 d in the same sector
      not-confirmed    another reduction covers the event and does not confirm it
      single-reduction no other reduction covers it (OWN-only or QLP-only sectors)
      noisy-star       the star's inverted control has a 2-reduction brightening with SNR >= this dip
  otherwise, with r_d = depth / reference depth and r_D = duration / reference duration (depth and
  duration from the fine-grid refit, 0.5-30 h, of each event and of the recovered catalogued event):
  (b) primary-like     0.5 <= r_d <= 2 and 1/3 <= r_D <= 3 (same shape as the catalogued event);
                       if a unique period is known and the event is off its ephemeris -> (c)
  (a) secondary-like   depth outside 0.5-2x the reference, duration within 1/3-3x; if a unique
                       period is known, the event must not fall on the primary ephemeris
  (c) other            duration outside 1/3-3x of the reference (different body, blend, or
                       variability), or a primary-like event off a known unique ephemeris
The control (inverted light curves) gets exactly the same known/extra and (d) filters.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PER = HERE / "data" / "per_star"
COMMON_DT, COMMON_N = 0.5, 2


def load():
    stars = []
    for f in sorted(PER.glob("*.json")):
        j = json.loads(f.read_text())
        if j.get("status") == "ok":
            stars.append(j)
        else:
            print("not ok:", j["TIC"], j.get("status"), j.get("error", "")[-200:])
    return stars


def is_known(ev, j):
    Dc = j["dur_h"] / 24
    for Tk in [j["event_T"]] + j["known_repeats"]:
        if abs(ev["T"] - Tk) < (ev["D"] + Dc) / 2 + 0.1:
            return True
    return False


def unique_period(j):
    p = [float(x) for x in j.get("periods", "").split()] if j.get("periods") else []
    return p[0] if len(p) == 1 else None


def reference(j):
    """Our measurement of the catalogued event: free-depth fit at the catalogued time and
    duration if it has SNR >= 5 (the note's own 'supported' level); else the merged detection at
    the event (fine-grid refit); else the CTOI values. (The fine-grid box refit alone can
    drift to long durations next to gaps, e.g. TIC 188835620, so it is the second choice.)"""
    Dc = j["dur_h"] / 24
    c = j.get("catalogued_event_free")
    if c and c["depth_ppm"] > 0 and c["snr"] >= 5:
        return c["depth_ppm"] * 1e-6, Dc, "free fit at catalogued time and duration"
    for ev in j["events"]:
        if abs(ev["T"] - j["event_T"]) < (ev["D"] + Dc) / 2 + 0.1:
            return ev.get("depth_fit", ev["depth"]), ev.get("D_fit", ev["D"]), "recovered event (refit)"
    return j["depth_ppm"] * 1e-6, Dc, "CTOI values"


def collect(stars, key):
    rows = []
    for j in stars:
        for ev in j[key]:
            rows.append(dict(TIC=j["TIC"], T=ev["T"], sector=ev["sector"], known=is_known(ev, j) if key == "events" else False,
                             ev=ev))
    return rows


def common_mode_counts(rows):
    """For each row, number of OTHER stars with a non-known row within COMMON_DT in the same sector."""
    df = pd.DataFrame([dict(TIC=r["TIC"], T=r["T"], sector=r["sector"]) for r in rows if not r["known"]])
    out = []
    for r in rows:
        if df.empty:
            out.append(0)
            continue
        m = (df.sector == r["sector"]) & ((df["T"] - r["T"]).abs() < COMMON_DT) & (df.TIC != r["TIC"])
        out.append(int(df[m].TIC.nunique()))
    return out


def alt_str(ev):
    return "; ".join(f"{k}:{v.get('status')}" + (f"({v['depth_ppm']:.0f}ppm,SNR{v['snr']})" if "snr" in v else "")
                     for k, v in ev.get("alt", {}).items())


def systematic_reason(ev, cm, ctl_conf_snr):
    if cm >= COMMON_N:
        return "common-mode"
    if ev["n_alt_tested"] >= 1 and ev["n_reductions"] < 2:
        return "not-confirmed"
    if ev["n_alt_tested"] == 0:
        return "single-reduction"
    if ctl_conf_snr is not None and ctl_conf_snr >= ev["snr"]:
        return "noisy-star"
    return None


def categorise(ev, j, ref):
    d_ref, D_ref, _ = ref
    dep, D = ev.get("depth_fit", ev["depth"]), ev.get("D_fit", ev["D"])
    r_d = dep / d_ref if d_ref > 0 else np.nan
    r_D = D / D_ref
    P = unique_period(j)
    phase = None
    on_eph = None
    if P:
        phase = ((ev["T"] - j["event_T"]) / P) % 1.0
        on_eph = min(phase, 1 - phase) * P < (D + D_ref) / 2 + 0.1
    if 0.5 <= r_d <= 2 and 1 / 3 <= r_D <= 3:
        if P and not on_eph:
            return "c", "primary-like shape but off the known period's ephemeris", r_d, r_D, phase
        return "b", "same depth (0.5-2x) and duration (1/3-3x) as the catalogued event" + (" (on known ephemeris)" if P else ""), r_d, r_D, phase
    if 1 / 3 <= r_D <= 3:
        if P and on_eph:
            return "b", "on the known primary ephemeris, different depth", r_d, r_D, phase
        return "a", ("deeper" if r_d > 2 else "shallower") + " than the catalogued event (outside 0.5-2x), similar duration", r_d, r_D, phase
    return "c", "duration differs by > 3x", r_d, r_D, phase


def main():
    stars = load()
    real = collect(stars, "events")
    ctl = collect(stars, "control_events")
    cm_real = common_mode_counts(real)
    # control common-mode is computed among control events (same rule)
    cm_ctl = common_mode_counts(ctl)
    by_tic = {j["TIC"]: j for j in stars}

    # control: per-star max SNR of a brightening confirmed in >= 2 reductions and not common-mode
    ctl_conf = {}
    for r, cm in zip(ctl, cm_ctl):
        ev = r["ev"]
        if ev["n_reductions"] >= 2 and cm < COMMON_N:
            ctl_conf[r["TIC"]] = max(ctl_conf.get(r["TIC"], 0), ev["snr"])

    # ------------------------------------------------ extra events (real)
    ex_rows = []
    for r, cm in zip(real, cm_real):
        if r["known"]:
            continue
        j, ev = by_tic[r["TIC"]], r["ev"]
        ref = reference(j)
        sysr = systematic_reason(ev, cm, ctl_conf.get(r["TIC"]))
        cat, why, r_d, r_D, phase = categorise(ev, j, ref)
        if sysr:
            cat, why = "d", sysr
        ex_rows.append(dict(
            TIC=j["TIC"], CTOI=j["CTOI"], note_cls=j["cls"], T_btjd=round(ev["T"], 4), sector=ev["sector"],
            product=ev["product"], depth_ppm=round(ev["depth"] * 1e6, 1), dur_h=round(ev["D"] * 24, 2),
            durs_detected_h=" ".join(str(x) for x in sorted(ev["durs_detected"])), snr=round(ev["snr"], 1),
            depth_fit_ppm=round(ev.get("depth_fit", np.nan) * 1e6, 1), dur_fit_h=round(ev.get("D_fit", np.nan) * 24, 2),
            snr_fit=round(ev.get("snr_fit", np.nan), 1),
            beta=round(ev["beta"], 2), n_reductions=ev["n_reductions"], n_alt_tested=ev["n_alt_tested"],
            independent_confirm=ev["independent_confirm"], in_two_reductions=ev["n_reductions"] >= 2,
            alt_products=alt_str(ev), n_other_stars_same_time=cm,
            ref_depth_ppm=round(ref[0] * 1e6, 1), ref_dur_h=round(ref[1] * 24, 2), ref_source=ref[2],
            depth_ratio=round(r_d, 3), dur_ratio=round(r_D, 3),
            known_period_d=unique_period(j), phase=None if phase is None else round(phase, 4),
            category=cat, category_reason=why))
    ex = pd.DataFrame(ex_rows)

    # ------------------------------------------------ control rows (same filters)
    c_rows = []
    for r, cm in zip(ctl, cm_ctl):
        ev = r["ev"]
        sysr = None
        if cm >= COMMON_N:
            sysr = "common-mode"
        elif ev["n_alt_tested"] >= 1 and ev["n_reductions"] < 2:
            sysr = "not-confirmed"
        elif ev["n_alt_tested"] == 0:
            sysr = "single-reduction"
        c_rows.append(dict(TIC=r["TIC"], T_btjd=round(ev["T"], 4), sector=ev["sector"], product=ev["product"],
                           height_ppm=round(ev["depth"] * 1e6, 1), dur_h=round(ev["D"] * 24, 2), snr=round(ev["snr"], 1),
                           n_reductions=ev["n_reductions"], n_alt_tested=ev["n_alt_tested"],
                           independent_confirm=ev["independent_confirm"], n_other_stars_same_time=cm,
                           filter=sysr or "passes"))
    cdf = pd.DataFrame(c_rows)
    cdf.to_csv(HERE / "data" / "control_events_all.csv", index=False)

    # ------------------------------------------------ per-star control false alarms
    fa = []
    for j in stars:
        cs = cdf[cdf.TIC == j["TIC"]] if len(cdf) else cdf
        es = ex[ex.TIC == j["TIC"]] if len(ex) else ex
        fa.append(dict(TIC=j["TIC"], note_cls=j["cls"], n_star_sectors=j["n_star_sectors"],
                       control_raw=len(cs), control_2red=int((cs.n_reductions >= 2).sum()) if len(cs) else 0,
                       control_passing_filters=int((cs["filter"] == "passes").sum()) if len(cs) else 0,
                       extra_raw=len(es), extra_2red=int(es.in_two_reductions.sum()) if len(es) else 0,
                       extra_passing_filters=int((es.category != "d").sum()) if len(es) else 0))
    fa = pd.DataFrame(fa)

    # ------------------------------------------------ recovery of catalogued event and known repeats
    rec = []
    for j in stars:
        Dc = j["dur_h"] / 24
        hit = [ev for ev in j["events"] if abs(ev["T"] - j["event_T"]) < (ev["D"] + Dc) / 2 + 0.1]
        nrep = sum(any(abs(ev["T"] - Tk) < (ev["D"] + Dc) / 2 + 0.1 for ev in j["events"]) for Tk in j["known_repeats"])
        c = j.get("catalogued_event_free") or {}
        rec.append(dict(TIC=j["TIC"], note_cls=j["cls"], event_recovered=bool(hit),
                        event_depth_ppm=round(hit[0]["depth"] * 1e6, 1) if hit else None,
                        event_dur_h=round(hit[0]["D"] * 24, 2) if hit else None, event_snr=round(hit[0]["snr"], 1) if hit else None,
                        ctoi_depth_ppm=j["depth_ppm"], ctoi_dur_h=j["dur_h"],
                        free_fit_snr_at_Dcat=round(c.get("snr", np.nan), 1), free_fit_fail_reason=c.get("reason"),
                        n_known_repeats=len(j["known_repeats"]), n_known_repeats_recovered=nrep))
    rec = pd.DataFrame(rec)

    # ------------------------------------------------ per-star summary + suggestion
    ps = []
    for j in stars:
        es = ex[ex.TIC == j["TIC"]] if len(ex) else ex
        n = {c: int((es.category == c).sum()) if len(es) else 0 for c in "abcd"}
        cls = j["cls"]
        sugg = cls
        if cls == "event-not-recovered":
            sugg = cls + (" (extra events listed separately)" if n["a"] + n["b"] + n["c"] else "")
        elif n["b"] and cls == "single":
            sugg = "repeat found -> duotransit candidate (P = dT/n), needs vetting"
        elif n["a"] and cls == "single":
            sugg = "EB candidate (secondary-like eclipse), needs vetting"
        elif n["c"] and cls == "single":
            sugg = "additional dip of different duration (other body/blend/variability), needs vetting"
        elif cls != "single" and (n["a"] or n["b"] or n["c"]):
            sugg = cls + " (+ extra events: " + ", ".join(f"{k}={v}" for k, v in n.items() if v and k != "d") + ")"
        r = rec[rec.TIC == j["TIC"]].iloc[0]
        ps.append(dict(TIC=j["TIC"], CTOI=j["CTOI"], note_cls=cls, n_star_sectors=j["n_star_sectors"],
                       event_recovered_by_free_search=bool(r.event_recovered),
                       n_known_repeats=int(r.n_known_repeats), n_known_repeats_recovered=int(r.n_known_repeats_recovered),
                       n_extra=len(es), n_a_secondary_like=n["a"], n_b_primary_like=n["b"], n_c_other=n["c"],
                       n_d_systematic=n["d"], n_control=int(fa[fa.TIC == j["TIC"]].control_raw.iloc[0]),
                       n_control_passing=int(fa[fa.TIC == j["TIC"]].control_passing_filters.iloc[0]),
                       suggested_class=sugg))
    ps = pd.DataFrame(ps)

    ex.to_csv(HERE / "extra_events_auto.csv", index=False)
    ps.to_csv(HERE / "data" / "per_star_summary_auto.csv", index=False)
    fa.to_csv(HERE / "control_false_alarms.csv", index=False)
    rec.to_csv(HERE / "recovery_check.csv", index=False)

    # ------------------------------------------------ numbers
    sample = ps[~ps.note_cls.isin(["event-not-recovered", "error"])]
    nr = ps[ps.note_cls == "event-not-recovered"]
    ss = fa.n_star_sectors.sum()
    s272 = fa[fa.TIC.isin(sample.TIC)].n_star_sectors.sum()
    out = dict(
        n_stars=len(ps), n_sample=len(sample), n_not_recovered=len(nr), star_sectors_all=int(ss), star_sectors_272=int(s272),
        control_raw=int(fa.control_raw.sum()), control_2red=int(fa.control_2red.sum()),
        control_passing=int(fa.control_passing_filters.sum()),
        control_raw_272=int(fa[fa.TIC.isin(sample.TIC)].control_raw.sum()),
        control_passing_272=int(fa[fa.TIC.isin(sample.TIC)].control_passing_filters.sum()),
        extra_raw_272=int(ex[ex.TIC.isin(sample.TIC)].shape[0]),
        extra_2red_272=int(ex[ex.TIC.isin(sample.TIC)].in_two_reductions.sum()),
        extra_by_cat_272=ex[ex.TIC.isin(sample.TIC)].category.value_counts().to_dict(),
        extra_d_reasons_272=ex[ex.TIC.isin(sample.TIC) & (ex.category == "d")].category_reason.value_counts().to_dict(),
        stars_with_extra_272=int((sample.n_extra > 0).sum()),
        stars_with_abc_272=int(((sample.n_a_secondary_like + sample.n_b_primary_like + sample.n_c_other) > 0).sum()),
        stars_by_cat_272={c: int((sample[f"n_{c}_" + {"a": "secondary_like", "b": "primary_like", "c": "other", "d": "systematic"}[c]] > 0).sum())
                          for c in "abcd"},
        single_with_a=int(((sample.note_cls == "single") & (sample.n_a_secondary_like > 0)).sum()),
        single_with_b=int(((sample.note_cls == "single") & (sample.n_b_primary_like > 0)).sum()),
        single_with_c=int(((sample.note_cls == "single") & (sample.n_c_other > 0)).sum()),
        single_suggestions=sample[sample.note_cls == "single"].suggested_class.value_counts().to_dict(),
        event_recovered_272=int(rec[rec.TIC.isin(sample.TIC)].event_recovered.sum()),
        known_repeats_272=int(rec[rec.TIC.isin(sample.TIC)].n_known_repeats.sum()),
        known_repeats_recovered_272=int(rec[rec.TIC.isin(sample.TIC)].n_known_repeats_recovered.sum()),
        event_recovered_34=int(rec[rec.TIC.isin(nr.TIC)].event_recovered.sum()),
        extra_by_cat_34=ex[ex.TIC.isin(nr.TIC)].category.value_counts().to_dict() if len(ex) else {},
    )
    (HERE / "data" / "report_numbers_auto.json").write_text(json.dumps(out, indent=1, default=int))
    print(json.dumps(out, indent=1, default=int))


if __name__ == "__main__":
    sys.exit(main())
