"""
Re-check of the singles project's candidate periods with the robust period finder (2026-10-12).

Why: pilot.periods_from_repeats proposes trial periods only from the repeat nearest the event and
needs every repeat on one linear ephemeris, so one spurious 'robust' repeat can delete the true
period or leave a wrong one (Fable 5.1 review, finding C2; k2singles/review_fable/
test_periods_from_repeats.py). The K2 project replaced it with k2lib.robust_periods: every repeat
can anchor trial periods and each period reports how many repeats it explains (n_support).

What this script does, for every JSON in singles/results/{pilot,batch2,batch3,full,rerun_*}/ with
n_repeats >= 1 (results are read, never modified):
  1. Rebuild the original grid: pilot.grid_windows on the cached light curves with OWN products,
     at the JSON's event_T (re-derived as pilot.run_target does when absent) and the catalogue
     D/depth. Only light-curve files that existed when the JSON was first committed are used
     (file mtime <= commit time), so later downloads and later own photometry do not leak in.
  2. Reproduction checks on that grid: sector count vs n_sectors, the repeat search
     (pilot.robust, same order and D separation as pilot.run_target) vs the stored repeats, and
     the old periods_from_repeats vs the stored candidate_periods.
  3. k2lib.robust_periods on the stored repeats. Periods that explain every repeat are the
     analogue of the old candidate list (n_robust_full). For stars with >= 2 repeats, periods
     that explain all but one repeat (n_robust_allbut1) show whether a single spurious repeat
     could hide a valid period.
  4. Class under the classify_final.py ladder, once with the stored candidate count and once with
     n_robust_full. The first must reproduce the frozen class; any difference in the second is a
     class change.

Rules fixed before looking at the output:
  period_set_changed   robust full-support set differs from the stored list (0.1% tolerance)
  class_changed        class with n_robust_full differs from the class with the stored list
  hidden_period        >= 2 repeats, no full-support period, but >= 1 period explains all but one

Usage:  python singles/scripts/recheck_periods_robust.py [--lc-dir PATH] [--tic TIC ...]
        (--lc-dir defaults to data/raw/lightcurves of this checkout, else of the main checkout)
Output: singles/results/period_recheck/period_recheck.csv (one row per JSON), printed summary.
"""
import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
sys.path.insert(0, str(ROOT / "k2singles" / "scripts"))
import alias_pipeline as ap  # noqa: E402
import classify_final as cf  # noqa: E402
import k2lib  # noqa: E402
import pilot  # noqa: E402

RES = ROOT / "singles" / "results"
OUT = RES / "period_recheck"
RUNS = ["pilot", "batch2", "batch3", "full"] + sorted(p.name for p in RES.glob("rerun_*") if p.is_dir())
FROZEN = RES / "final_classification_2026-10-11_rnaas.csv"
P_TOL = 1e-3
CHECK_TICS = {159614322: "89.965 / 179.930 d (hand analysis incl. S21 partial)",
              147035057: "37.31497 d (hand analysis incl. own S108)",
              229055790: "1838.945/n >= 26.7 d, 38 aliases"}


def default_lc_dir():
    here = ROOT / "data" / "raw" / "lightcurves"
    if here.exists() and any(here.iterdir()):
        return here
    common = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--git-common-dir"],
                            capture_output=True, text=True, check=True).stdout.strip()
    return (ROOT / common).resolve().parent / "data" / "raw" / "lightcurves"


def first_commit_time(f):
    out = subprocess.run(["git", "-C", str(ROOT), "log", "--diff-filter=A", "--format=%ct", "--",
                          str(f.relative_to(ROOT))], capture_output=True, text=True).stdout.split()
    return float(out[-1]) if out else None


def load_asof(lc_dir, tic, cutoff, tmp):
    """ap.load_lightcurves on a copy of the star's cache that holds only files with mtime <= cutoff."""
    src, dst = lc_dir / str(tic), Path(tmp) / str(tic)
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    dropped = []
    for f in src.glob("s*.csv.gz"):
        if cutoff is None or f.stat().st_mtime <= cutoff:
            shutil.copy2(f, dst / f.name)
        else:
            dropped.append(f.name)
    ap.LC = Path(tmp)
    return ap.load_lightcurves(tic, products=pilot.PRODUCTS), sorted(dropped)


def detect_repeats(g, delta, D):
    """The repeat search of pilot.run_target, on a rebuilt grid."""
    reps = []
    for _, r in g[~g.near_event].sort_values("snr", ascending=False).iterrows():
        if pilot.robust(r, delta) and all(abs(r.T_pred - x) > D for x in reps):
            reps.append(r.T_pred)
    return sorted(reps)


def same_periods(a, b):
    a, b = sorted(a), sorted(b)
    return len(a) == len(b) and all(abs(x - y) / y < P_TOL for x, y in zip(a, b))


def same_times(a, b, tol):
    a, b = sorted(a), sorted(b)
    return len(a) == len(b) and all(abs(x - y) < tol for x, y in zip(a, b))


def star_class(j, ncand, gaia):
    """The classify_final.main decision ladder, with the candidate-period count as a parameter."""
    tic = int(j["TIC"])
    imgs = [d for d in (j.get("diff_images") or []) if isinstance(d, dict)]
    sv, _ = cf.source_verdict(imgs)
    g = gaia.loc[tic] if tic in gaia.index else None
    gaia_bin = bool(g is not None and ((g.get("non_single_star", 0) or 0) > 0 or bool(g.get("rv_variable"))))
    nrep = j.get("n_repeats") or 0
    if j.get("status") == "error":
        return "error"
    if j.get("event_verdict") not in ("supported",):
        return "event-not-recovered"
    if sv == "off-target (firm)":
        return "off-target-EB" if nrep else "off-target (single event)"
    if nrep and gaia_bin:
        return "EB-on-target (Gaia)"
    if nrep and ncand == 1:
        return "periodic-on-target" if sv in ("on-target", "inconclusive") else "periodic (source " + sv + ")"
    if nrep and ncand > 1:
        return "duotransit / narrowed"
    if nrep:
        return "repeats-inconsistent"
    return "single"


def fmt(ps, k=8):
    return " ".join(f"{p:.4f}" for p in sorted(ps, reverse=True)[:k]) + (" ..." if len(ps) > k else "")


def recheck(f, run, lc_dir, tmp, gaia):
    j = json.loads(f.read_text())
    tic = int(j["TIC"])
    D, delta, epoch = j["dur_h"] / 24, j["depth_ppm"] * 1e-6, float(j["epoch_btjd"])
    cutoff = first_commit_time(f)
    lc, dropped = load_asof(lc_dir, tic, cutoff, tmp)
    t, fl, sg = lc.time.to_numpy(), lc.flux.to_numpy(), lc.segment.to_numpy()
    ev_T = j.get("event_T")
    if ev_T is None:
        ev = ap.test_window(t, fl, epoch, pilot.EPOCH_ERR, D, delta, sg)
        ev_T = epoch + (ev.get("offset", 0.0) or 0.0)
    g = pilot.grid_windows(lc, D, delta, ev_T)
    reps = j["repeats"]
    nrep = len(reps)
    stored = [p for _, p in (j.get("candidate_periods") or [])]
    old = [p for _, p in pilot.periods_from_repeats(g, ev_T, reps, D)]
    redet = detect_repeats(g, delta, D)
    rp = k2lib.robust_periods(g, ev_T, reps, D)
    full = [d["P"] for d in rp if d["n_support"] == nrep]
    allbut1 = [d for d in rp if nrep >= 2 and d["n_support"] == nrep - 1]
    max_sup = max([d["n_support"] for d in rp], default=0)
    # which repeat each all-but-one period leaves out
    left_out = sorted({next(x["T"] for x in reps if x["T"] not in d["support_T"]) for d in allbut1})
    cls_old = star_class(j, len(stored), gaia)
    cls_new = star_class(j, len(full), gaia)
    return dict(
        TIC=tic, run=run, event_verdict=j.get("event_verdict"), n_repeats=nrep,
        repeats=" ".join(f"{x['T']:.3f}:S{x['sector']}:{x['product']}:snr{x['snr']}" for x in reps),
        n_sectors_json=j.get("n_sectors"), n_sectors_rebuilt=int(lc.sector.nunique()),
        files_after_cutoff=" ".join(dropped),
        repeats_reproduced=same_times(redet, [x["T"] for x in reps], D / 8),
        old_periods_reproduced=same_periods(old, stored),
        n_candidate_periods_stored=len(stored), periods_stored=fmt(stored),
        n_robust_full=len(full), periods_robust_full=fmt(full),
        n_robust_allbut1=len(allbut1), periods_robust_allbut1=fmt([d["P"] for d in allbut1]),
        allbut1_left_out_T=" ".join(f"{T:.3f}" for T in left_out),
        max_support=max_sup, n_robust_any=len(rp),
        period_set_changed=not same_periods(full, stored),
        class_stored=cls_old, class_robust=cls_new, class_changed=cls_old != cls_new,
        hidden_period=nrep >= 2 and not full and bool(allbut1),
        check_note=CHECK_TICS.get(tic, ""))


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--lc-dir", type=Path, default=None)
    ap_.add_argument("--tic", type=int, nargs="*")
    a = ap_.parse_args()
    lc_dir = a.lc_dir or default_lc_dir()
    print(f"light curves: {lc_dir}")
    gaia = pd.read_csv(ROOT / "singles" / "data" / "gaia_binarity_flags.csv").set_index("tic")
    frozen = pd.read_csv(FROZEN).set_index("TIC")
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for run in RUNS:
            for f in sorted((RES / run).glob("*.json")):
                j = json.loads(f.read_text())
                if not (j.get("n_repeats") or 0) >= 1 or (a.tic and int(j["TIC"]) not in a.tic):
                    continue
                r = recheck(f, run, lc_dir, tmp, gaia)
                fz = frozen.loc[r["TIC"]] if r["TIC"] in frozen.index else None
                r["in_frozen_note"] = bool(fz is not None and fz["run"] == run)
                r["class_frozen"] = fz["cls"] if r["in_frozen_note"] else ""
                rows.append(r)
                print(f"TIC {r['TIC']:>10} {run:<17} reps {r['n_repeats']:>2} | stored {r['n_candidate_periods_stored']:>3}"
                      f" robust-full {r['n_robust_full']:>3} all-but-1 {r['n_robust_allbut1']:>3} | "
                      f"repro reps {r['repeats_reproduced']!s:<5} periods {r['old_periods_reproduced']!s:<5} | "
                      f"{r['class_stored']} -> {r['class_robust']}", flush=True)
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "period_recheck.csv", index=False)
    fz = df[df.in_frozen_note]
    print(f"\nJSONs rechecked: {len(df)}; used by the frozen RNAAS classification: {len(fz)}")
    print(f"frozen class reproduced by the ladder: {(fz.class_stored == fz.class_frozen).sum()}/{len(fz)}")
    print(f"old periods reproduced on the rebuilt grid: {df.old_periods_reproduced.sum()}/{len(df)}; "
          f"repeats reproduced: {df.repeats_reproduced.sum()}/{len(df)}")
    print(f"period set changed: {df.period_set_changed.sum()}; class changed: {df.class_changed.sum()}; "
          f"hidden period (all but one repeat): {df.hidden_period.sum()}")
    print("\nfrozen-note class counts, stored -> robust:")
    print(pd.concat([fz.class_stored.value_counts().rename("stored"),
                     fz.class_robust.value_counts().rename("robust")], axis=1).fillna(0).astype(int).to_string())


if __name__ == "__main__":
    main()
