"""
Partial-transit sweep: look for transits cut off by data gaps, at every data-segment edge,
for every analysed single-transit CTOI.

Why: the repeat search (pilot.py) needs baseline on both sides of a window, so a transit whose
egress (or ingress) falls in a data gap is never counted. TIC 159614322's third transit was
exactly that: ingress plus ~9 h of floor in S21, ending at the orbit-1 downlink gap.

Rules (fixed before looking at any results, 2026-10-11):
  Edges     the end (and, mirrored, the start) of every contiguous data segment
            (gap > 0.5 d), in each product separately: SPOC_120, TESS-SPOC, QLP, OWN.
  Model     one-sided: flux = c0 + c1 (t - t_step) - A [t past the step], fitted on a baseline
            stretch before the step, plus the in-transit stretch between the step and the edge.
            The step time is scanned so that 25-90% of the catalogued duration D lies in the data,
            with at least 1.5 h of in-transit data.
            Baseline: max(1.0 d, 2 x in-transit span), at least 10 points and 0.4 d.
  Noise     red-noise aware: residuals of the baseline binned at w = max(in-span/4, 3 cadences);
            sigma_A = sigma_bin * sqrt(1/n_bins_in + 1/n_bins_base).
  Per-product detection
            A in [0.5, 2] x catalogued depth, SNR >= 5, and a flat transit floor: a linear fit to
            the in-transit residuals changes by less than 0.5 A across the in-transit span.
  Candidate ('2-product')
            detections in >= 2 distinct products, with step times within max(D/4, 1 h).
  'single-product'
            SNR >= 8 in the only product covering that time (weaker; listed separately).
  Exclusions
            steps within D of the catalogued event or of a known (full) repeat.
  Control   the identical search for brightenings (A < 0, same |A| range), run on the same data.
            Systematics produce steps of both signs, so its count estimates the false-alarm rate.
  Positive control
            TIC 159614322 S21 (step at BTJD ~1883.57) must be recovered.

  REVISED after the positive control failed (2026-10-11, before any other star was looked at):
    - in-transit span allowed up to 1.2 D (was 0.9 D): catalogued durations are rough
      (TIC 159614322: catalogue 8.8 h, fitted 9.7 h);
    - the best step is chosen among steps that pass the flat-floor test (it was chosen first and
      tested after, so a slightly mis-placed step failed);
    - the segment end may be trimmed by 0, 0.5 or 1 h, because cadences right at a downlink gap
      are often degraded.
  The brightening control uses the identical rules, so it still measures the false-alarm rate.

Output: singles/results/partial_sweep/<TIC>.json, summary.csv, stdout log.
Usage:  python singles/scripts/partial_sweep.py [n_workers]
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(HERE))
OUT = ROOT / "singles" / "results" / "partial_sweep"
RES = ROOT / "singles" / "results"
RUNS = ["pilot", "batch2", "batch3", "full"]
PRODS = ["SPOC_120", "TESS-SPOC", "QLP", "OWN"]
SNR_MIN, SNR_SINGLE, DMIN, DMAX = 5.0, 8.0, 0.5, 2.0


def edge_scan(t, f, te, D, sign):
    """t, f: one segment (sorted). Step before the END te of the segment (optionally trimmed by
    up to 1 h). sign=+1 dip, -1 bump. Returns the best (by SNR) flat-floor step, else the best
    step overall (marked flat=False), or None."""
    best, best_any = None, None
    for trim in (0.0, 0.5 / 24, 1.0 / 24):
        keep = t <= te - trim + 1e-9
        r = _edge_scan(t[keep], f[keep], te - trim, D, sign)
        for x in r:
            x["trim_h"] = trim * 24
            if best_any is None or x["snr"] > best_any["snr"]:
                best_any = x
            if x["flat"] and (best is None or x["snr"] > best["snr"]):
                best = x
    return best if best is not None else best_any


def _edge_scan(t, f, te, D, sign):
    out = []
    if len(t) < 30:
        return out
    cad = np.median(np.diff(t)) if len(t) > 2 else 0.02
    hi_span = 1.2 * D
    lo_span = max(0.25 * D, 1.5 / 24)
    if lo_span > hi_span:
        return out
    for span in np.arange(lo_span, hi_span + 1e-9, max(min(D / 12, 0.5 / 24), cad)):
        ts = te - span
        Lb = max(1.0, 2 * span)
        base = (t >= ts - Lb) & (t < ts - 0.02)
        inn = t >= ts
        if base.sum() < 10 or inn.sum() < 5 or (t[base].max() - t[base].min()) < 0.4:
            continue
        m = base | inn
        A_mat = np.column_stack([np.ones(m.sum()), t[m] - ts, -(t[m] >= ts).astype(float)])
        coef, *_ = np.linalg.lstsq(A_mat, f[m], rcond=None)
        A = coef[2] * sign
        resid = f[m] - A_mat @ coef
        rb = resid[(t[m] < ts)]
        tb = t[m][t[m] < ts]
        w = max(span / 4, 3 * cad)
        nb_base = int((tb.max() - tb.min()) / w)
        if nb_base < 4:
            continue
        idx = np.floor((tb - tb.min()) / w).astype(int)
        binned = np.array([rb[idx == i].mean() for i in np.unique(idx) if (idx == i).sum() >= 2])
        if len(binned) < 4:
            continue
        sig_bin = np.std(binned, ddof=1)
        nb_in = max(span / w, 1)
        sA = sig_bin * np.sqrt(1 / nb_in + 1 / len(binned))
        snr = A / sA if sA > 0 else 0
        ri = resid[t[m] >= ts]
        ti = t[m][t[m] >= ts]
        slope = np.polyfit(ti - ti.mean(), ri, 1)[0] if len(ti) > 3 else 0.0
        flat = abs(slope * span) < 0.5 * abs(A) if A != 0 else False
        out.append(dict(t_step=float(ts), span_h=float(span * 24), A=float(A), snr=float(snr), flat=bool(flat),
                        n_in=int(inn.sum()), n_base=int(base.sum())))
    return out


def scan_product(lc, D, sign, avoid):
    out = []
    for sgid, s in lc.groupby("segment"):
        t, f = s.time.to_numpy(), s.flux.to_numpy()
        if len(t) < 30:
            continue
        for direction in ("end", "start"):
            if direction == "end":
                tt, ff, te = t, f, t[-1]
            else:          # mirror: the egress falls inside the data after a gap
                tt, ff, te = -t[::-1], f[::-1], -t[0]
            r = edge_scan(tt, ff, te, D, sign)
            if r is None:
                continue
            t_step = r["t_step"] if direction == "end" else -r["t_step"]
            t_mid = t_step + D / 2 if direction == "end" else t_step - D / 2
            if any(abs(t_mid - a) < D for a in avoid):
                continue
            r.update(direction=direction, t_step=t_step, t_mid_est=t_mid, sector=int(s.sector.iloc[0]))
            out.append(r)
    return out


def work(args):
    tic, D_h, depth_ppm, ev_T, rep_T, cand_P = args
    import alias_pipeline as ap
    D, delta = D_h / 24, depth_ppm * 1e-6
    avoid = [ev_T] + list(rep_T)
    res = dict(TIC=tic, D_h=D_h, depth_ppm=depth_ppm, event_T=ev_T, dips=[], bumps=[], edges_tested=0)
    per = {}
    for p in PRODS:
        lc = ap.load_lightcurves(tic, products=[p])
        if len(lc) == 0:
            continue
        per[p] = {}
        for sign, key in ((+1, "dips"), (-1, "bumps")):
            rows = scan_product(lc, D, sign, avoid)
            res["edges_tested"] += len(rows) if sign == 1 else 0
            det = [dict(r, product=p) for r in rows
                   if DMIN * delta <= r["A"] <= DMAX * delta and r["snr"] >= SNR_MIN and r["flat"]]
            per[p][key] = det
        per[p]["coverage"] = (float(lc.time.min()), float(lc.time.max()))
        per[p]["segments"] = [(float(s.time.min()), float(s.time.max())) for _, s in lc.groupby("segment")]
    tol = max(D / 4, 1 / 24)
    for key in ("dips", "bumps"):
        allx = [x for p in per for x in per[p].get(key, [])]
        allx.sort(key=lambda x: x["t_step"])
        used = set()
        for i, x in enumerate(allx):
            if i in used:
                continue
            grp = [x] + [y for j, y in enumerate(allx) if j > i and j not in used and abs(y["t_step"] - x["t_step"]) < tol]
            for j, y in enumerate(allx):
                if y in grp[1:]:
                    used.add(j)
            prods = sorted({g["product"] for g in grp})
            # which other products had data at that time at all?
            covering = [p for p in per if any(a - 0.05 <= x["t_mid_est"] <= b + 0.05 or a <= x["t_step"] <= b
                                                for a, b in per[p]["segments"])]
            cls = ("2-product" if len(prods) >= 2 else
                   "single-product" if (len(covering) <= 1 and max(g["snr"] for g in grp) >= SNR_SINGLE) else "1-of-many")
            ent = dict(t_step=float(np.mean([g["t_step"] for g in grp])), t_mid_est=float(np.mean([g["t_mid_est"] for g in grp])),
                       products=prods, covering=covering, cls=cls,
                       max_snr=float(max(g["snr"] for g in grp)), depth_ppm=[round(g["A"] * 1e6) for g in grp],
                       direction=grp[0]["direction"], sector=grp[0]["sector"])
            if cand_P and key == "dips":
                ent["matches_period"] = [P for P in cand_P
                                         if min(abs(((ent["t_mid_est"] - ev_T) / P) % 1), 1 - abs(((ent["t_mid_est"] - ev_T) / P) % 1)) * P < D / 2 + 0.05]
            res[key].append(ent)
    (OUT / f"{tic}.json").write_text(json.dumps(res, indent=1))
    n2 = sum(1 for d in res["dips"] if d["cls"] == "2-product")
    b2 = sum(1 for d in res["bumps"] if d["cls"] == "2-product")
    return tic, n2, b2, res


def targets():
    cls = pd.read_csv(RES / "final_classification.csv")
    keep = cls[~cls.cls.isin(["error", "event-not-recovered"])]
    rows = []
    for _, r in keep.iterrows():
        f = RES / r.run / f"{int(r.TIC)}.json"
        j = json.loads(f.read_text())
        rows.append((int(r.TIC), float(j["dur_h"]), float(j["depth_ppm"]), float(j["event_T"]),
                     [x["T"] for x in j.get("repeats") or []],
                     [p for _, p in (j.get("candidate_periods") or [])]))
    return rows


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    OUT.mkdir(parents=True, exist_ok=True)
    rows = targets()
    print(f"stars: {len(rows)}", flush=True)
    summ = []
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(work, a): a[0] for a in rows}
        for i, fu in enumerate(as_completed(futs), 1):
            try:
                tic, n2, b2, res = fu.result()
            except Exception as e:  # noqa: BLE001
                print(f"[{i}/{len(rows)}] TIC {futs[fu]}: ERROR {str(e)[:150]}", flush=True)
                continue
            for key in ("dips", "bumps"):
                for d in res[key]:
                    summ.append(dict(TIC=tic, kind=key, **{k: (";".join(v) if isinstance(v, list) and v and isinstance(v[0], str) else v)
                                                          for k, v in d.items()}))
            if n2 or b2:
                print(f"[{i}/{len(rows)}] TIC {tic}: 2-product dips {n2}, bumps {b2}", flush=True)
            elif i % 25 == 0:
                print(f"[{i}/{len(rows)}] ...", flush=True)
    s = pd.DataFrame(summ)
    s.to_csv(OUT / "summary.csv", index=False)
    if len(s):
        print(s.groupby(["kind", "cls"]).size().to_string())
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
