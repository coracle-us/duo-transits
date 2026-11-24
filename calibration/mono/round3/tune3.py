"""
Round 3 tuning: score rule variants (rule3.Rule3) on the stored already-seen trials (dump3.py).

Per variant and trial set it reports the pre-registered metrics: wrong period (all called trials,
and those whose alias event is a true transit), false events, repeat recovery, and alias power
(median fraction of aliases excluded in called true-event trials, beside round 2's table through
the catalogued time and the same event).

Inputs  : calibration/mono/work/round3_tuning[_noisy]/<tic>.pkl
Outputs : calibration/mono/round3/tuning/<name>.csv (one row per variant and set),
          calibration/mono/work/round3_tune_rows/<name>.parquet (one row per variant and trial)
Usage   : python calibration/mono/round3/tune3.py <family> [--noisy]
"""
import argparse
import multiprocessing as mp
import pickle
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "round2"))
import rule3  # noqa: E402
from rule3 import Rule3, ROUND2  # noqa: E402

WORK = HERE.parent / "work"
TEST_A = set(pd.read_csv(HERE.parent / "round2/test/arm_a_defs.csv").tic)
TEST_B = set(pd.read_csv(HERE.parent / "round2/test/arm_b_defs.csv").tic)


def families():
    f = {}
    f["base"] = {"round2": ROUND2}
    t1 = {"round2": ROUND2}
    for c in (1, 2, 3):
        for cap in (1.0, 2.0, 4.0, np.inf):
            for sds in (0.0, 0.5):
                t1[f"c{c}_cap{cap}_sds{sds}"] = Rule3(c_drift=c, cap_dur=cap, sig_dur_snr=sds)
    f["t1"] = t1
    f["t1ref"] = {k + "_ref": replace(v, use_ref_time=True) for k, v in t1.items()}
    b = Rule3(c_drift=2, cap_dur=2.0)
    t5 = {"round2": ROUND2, "c2cap2": b}
    for sec in (0, 8, 9, 10):
        for cp, frac in ((False, 0), (True, 0.0), (True, 0.5), (True, 1.0)):
            for lo in (0, 10, 12):
                t5[f"sec{sec}_cp{int(cp)}_{frac}_low{lo}"] = replace(
                    b, sec_box_snr_min=sec, common_period=cp, cp_pmin_frac=frac,
                    lowref_snr=7.0 if lo else 0.0, lowref_event_snr=lo)
    f["t5"] = t5
    f["t3"] = {"c2cap2": b, "pmin_off": replace(b, pmin_mode="off"), "pmin_cat": replace(b, pmin_mode="cat"),
               "pmin_0.5": replace(b, pmin_factor=0.5), "pmin_cat_0.5": replace(b, pmin_mode="cat", pmin_factor=0.5)}
    f["final"] = {"round2": ROUND2, "c2cap2": b, "sec8": replace(b, sec_box_snr_min=8),
                  "sec9": replace(b, sec_box_snr_min=9), "sec10": replace(b, sec_box_snr_min=10)}
    return f


class Boxes:
    def __init__(self, t, f):
        import harness as h2
        self.b = h2.Boxes(t, f)

    def __call__(self, c, dur):
        return self.b(c, dur)


def own_removed(tr):
    ev, null = tr["events"], tr["null_events"]
    if null is None:
        return None
    nt = np.array([e["tc"] for e in null])
    return [e for e in ev if e.get("at_truth", False) or not (nt.size and np.min(np.abs(nt - e["tc"])) < tr["dur"])]


def star(args):
    path, variants = args
    rows = []
    for tr in pickle.load(open(path, "rb")):
        d = tr["d"]
        box = Boxes(tr["T"], tr["F"]) if tr["T"] is not None else None
        r2cache = {}

        def r2frac(tc):
            if tc not in r2cache:
                t1, t2 = sorted([d["t_cat"], tc])
                ev = next(e for e in tr["events"] if e["tc"] == tc)
                rr = rule3.alias_table3(t1, t2, 0.0, box, tr["tmin"], tr["tmax"], tr["dur"], tr["depth0"],
                                        ev["depth_err"], tr["sig"], tr["ref_measured"], bool(d["multi"]),
                                        tr["pmin"], ROUND2)
                r2cache[tc] = float(np.mean([x[2] == "excluded" for x in rr])) if rr else np.nan
            return r2cache[tc]

        for name, r in variants.items():
            base = dict(variant=name, arm=d["arm"], tic=d["tic"], id=d.get("toi", d.get("draw")), src=d.get("src"),
                        multi=bool(d["multi"]), ref_snr=tr["ref_snr"], n_clean_repeats=tr["n_clean_repeats"],
                        truth_P=d["truth_P"])
            for tag, ev in (("", tr["events"]), ("_own", own_removed(tr))):
                if ev is None:
                    continue
                o = rule3.call3(ev, tr, r, box, d["truth_P"])
                if o["call"] in ("resolved", "narrowed") and o["event_true"]:
                    o["frac_alias_excluded_r2"] = r2frac(o["t2"])
                base.update({k + tag: v for k, v in o.items()})
            rows.append(base)
    return rows


def rate(k, n):
    return f"{int(k)}/{int(n)}"


def metrics(df):
    out = []
    for name, g in df.groupby("variant", sort=False):
        for label, s, tag in (
                ("A tune", g[g.arm.isin(["A1", "A2"]) & ~g.multi & ~g.tic.isin(TEST_A)], ""),
                ("A test2", g[g.arm.isin(["A1", "A2"]) & ~g.multi & g.tic.isin(TEST_A)], ""),
                ("A all", g[g.arm.isin(["A1", "A2"]) & ~g.multi], ""),
                ("B tune", g[(g.arm == "B") & ~g.tic.isin(TEST_B)], "_own"),
                ("B test2", g[(g.arm == "B") & g.tic.isin(TEST_B)], "_own"),
                ("B all", g[g.arm == "B"], "_own")):
            if not len(s):
                continue
            fe = s[s.arm == "A1"] if label[0] == "A" else s
            x = s[s.arm == "A2"] if label[0] == "A" else s
            called = x[x["call" + tag].isin(["resolved", "narrowed"])]
            wrong = called[~called["truth_in_alias_set" + tag].astype(bool)]
            ct = called[called["event_true" + tag].astype(bool)]
            wt = ct[~ct["truth_in_alias_set" + tag].astype(bool)]
            rep = x[x.n_clean_repeats > 0]
            out.append(dict(
                variant=name, set=label,
                false_event=rate(fe["false_event" + tag].sum(), len(fe)),
                wrong=rate(len(wrong), len(called)), wrong_true=rate(len(wt), len(ct)),
                recovery=rate(rep["repeat_detected" + tag].sum(), len(rep)),
                resolved=int((ct["call" + tag] == "resolved").sum()),
                power=round(float(ct["frac_alias_excluded" + tag].median()), 4) if len(ct) else np.nan,
                power_r2=round(float(ct["frac_alias_excluded_r2" + tag].median()), 4) if len(ct) else np.nan,
                med_n_ok=float(ct["n_alias_ok" + tag].median()) if len(ct) else np.nan))
    return pd.DataFrame(out)


def run(variants, src_dir, workers=12):
    paths = sorted(Path(src_dir).glob("*.pkl"))
    with mp.get_context("spawn").Pool(workers) as pool:
        rows = [r for rs in pool.imap_unordered(star, [(p, variants) for p in paths], chunksize=4) for r in rs]
    return pd.DataFrame(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("family")
    ap.add_argument("--noisy", action="store_true")
    a = ap.parse_args()
    v = families()[a.family]
    tag = a.family + ("_noisy" if a.noisy else "")
    df = run(v, WORK / ("round3_tuning_noisy" if a.noisy else "round3_tuning"))
    (WORK / "round3_tune_rows").mkdir(exist_ok=True)
    df.to_pickle(WORK / "round3_tune_rows" / f"{tag}.pkl")
    m = metrics(df)
    (HERE / "tuning").mkdir(exist_ok=True)
    m.to_csv(HERE / "tuning" / f"{tag}.csv", index=False)
    print(m.to_string(index=False))
