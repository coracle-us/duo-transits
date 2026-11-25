"""
tess_wide step 7 (PLAN.md section 9): the pilot report, counts only (no star is named), and the all-sector projection.

Reads the s<NNN>_*.json/csv outputs of steps 1-6 and writes tess_wide/PILOT_S<NNN>.md and
tess_wide/data/s<NNN>_report.json. Completeness: injections binned by true depth (5 log bins, 300-15,000 ppm) and
duration (4 log bins, 2-13 h), at detection and end to end (5.1-5.3 + CNN pass; true Rp > 2 R_J left out of the
end-to-end denominator), with binomial (Wilson) errors.

Projection: QLP files per sector from the size of each sector's MAST bulk script (HEAD requests; bytes per line
calibrated on the pilot sector's script), times the pilot's star-cut pass fraction; file size by cadence era measured
on MAST 2026-11-25 (S10 77,760 B at 30 min, S40 230,400 B at 10 min, S70 653,760 B and S104 682,560 B at 200 s);
per-star times from the pilot.
Usage:  python tess_wide/s7_report.py --sector 104
"""
import argparse
import json
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
DEPTH_EDGES = np.geomspace(300, 15000, 6)
DUR_EDGES = np.geomspace(2, 13, 5)
LAST_QLP = 104


def era_bytes(s):
    return 77760 if s <= 26 else 230400 if s <= 55 else 682560


def wilson(k, n, z=1.0):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, h


def grid(tr, ok, mask):
    rows = []
    for i in range(len(DEPTH_EDGES) - 1):
        for j in range(len(DUR_EDGES) - 1):
            m = mask & tr.depth_true.ge(DEPTH_EDGES[i] * 1e-6) & tr.depth_true.lt(DEPTH_EDGES[i + 1] * 1e-6) \
                & tr.dur_h.ge(DUR_EDGES[j]) & tr.dur_h.lt(DUR_EDGES[j + 1])
            p, h = wilson(int((ok & m).sum()), int(m.sum()))
            rows.append(dict(depth_lo=DEPTH_EDGES[i], depth_hi=DEPTH_EDGES[i + 1], dur_lo=DUR_EDGES[j],
                             dur_hi=DUR_EDGES[j + 1], n=int(m.sum()), k=int((ok & m).sum()), frac=p, err=h))
    return pd.DataFrame(rows)


def table(g):
    lines = ["| depth (ppm) \\ duration (h) | " + " | ".join(f"{a:.1f}–{b:.1f}" for a, b in
                                                          zip(DUR_EDGES[:-1], DUR_EDGES[1:])) + " |",
             "|---|" + "---|" * (len(DUR_EDGES) - 1)]
    for i in range(len(DEPTH_EDGES) - 1):
        cells = []
        for j in range(len(DUR_EDGES) - 1):
            r = g.iloc[i * (len(DUR_EDGES) - 1) + j]
            cells.append("–" if r.n == 0 else f"{100 * r.frac:.0f} ± {100 * r.err:.0f} % ({int(r.k)}/{int(r.n)})")
        lines.append(f"| {DEPTH_EDGES[i]:,.0f}–{DEPTH_EDGES[i + 1]:,.0f} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def script_sizes(s_pilot):
    f = DATA / "qlp_script_sizes.json"
    sizes = json.loads(f.read_text()) if f.exists() else {}
    for s in range(1, LAST_QLP + 1):
        if str(s) in sizes:
            continue
        u = f"https://archive.stsci.edu/hlsps/qlp/download_scripts/hlsp_qlp_tess_ffi_s{s:04d}_tess_v01_llc-fits.sh"
        # curl, not urllib: urllib HEAD requests to MAST failed for every sector on 2026-11-25
        r = subprocess.run(["curl", "-s", "-I", "--max-time", "60", u], capture_output=True, text=True)
        m = re.search(r"(?im)^content-length:\s*(\d+)", r.stdout) if " 200" in r.stdout.split("\n")[0] else None
        sizes[str(s)] = int(m.group(1)) if m else None
    f.write_text(json.dumps(sizes, indent=1))
    return sizes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sector", type=int, default=104)
    s = ap.parse_args().sector
    J = lambda n: json.loads((DATA / f"s{s:04d}_{n}.json").read_text())
    cut, dl, se, stg, cnn, xm = J("starcut"), J("download"), J("search"), J("stages"), J("cnn"), J("crossmatch")
    # completeness
    tr = pd.read_csv(DATA / f"s{s:04d}_inj_truth.csv")
    ci = pd.read_csv(DATA / f"s{s:04d}_cnn_inj.csv")
    v = tr.tic.map(ci.set_index("tic").cnn)
    tr["rec_e2e"] = tr.rec_stage & v.eq("pass")
    tr["rec_e2e_or"] = tr.rec_stage & v.isin(["pass", "out_of_range"])
    allm = pd.Series(True, index=tr.index)
    e2e_mask = tr.rp_true_rj <= 2.0
    g_det, g_e2e = grid(tr, tr.rec_detect, allm), grid(tr, tr.rec_e2e, e2e_mask)
    g_det.to_csv(DATA / f"s{s:04d}_completeness_detect.csv", index=False)
    g_e2e.to_csv(DATA / f"s{s:04d}_completeness_e2e.csv", index=False)
    tot = dict(detect=wilson(int(tr.rec_detect.sum()), len(tr)),
               e2e=wilson(int(tr.rec_e2e[e2e_mask].sum()), int(e2e_mask.sum())),
               e2e_or=wilson(int(tr.rec_e2e_or[e2e_mask].sum()), int(e2e_mask.sum())))
    # stage counts
    R, C = stg["real"], stg["ctrl"]
    cr, cc = cnn["real"]["verdicts"], cnn["ctrl"]["verdicts"]
    final_real, final_ctrl = cr.get("pass", 0), cc.get("pass", 0)
    control_level = final_real == 0 or final_ctrl >= 0.5 * final_real
    # PLAN 5.4 flags (not cuts) among CNN passes, real vs control
    flags = {}
    for kind in ("real", "ctrl"):
        st = pd.read_csv(DATA / f"s{s:04d}_stage_{kind}.csv")
        cn = pd.read_csv(DATA / f"s{s:04d}_cnn_{kind}.csv")
        m = st[st.survivor].reset_index(drop=True).assign(cnn=cn.cnn.values)
        p_ = m[m.cnn == "pass"]
        flags[kind] = dict(n=len(p_), bkg=int((p_.bkg_sig.abs() > 3).sum()),
                           centroid=int(((p_.x_sig.abs() > 3) | (p_.y_sig.abs() > 3)).sum()),
                           either=int(((p_.bkg_sig.abs() > 3) | (p_.x_sig.abs() > 3) | (p_.y_sig.abs() > 3)).sum()),
                           camera=p_.camera.value_counts().sort_index().to_dict())
    # cost and projection
    n_star = cut["passed"]
    prev = dl.get("previous", {})
    n_dl = dl.get("downloaded_now", dl["stars"]) + prev.get("downloaded_now", prev.get("stars", 0))
    s_dl = dl["seconds"] + prev.get("seconds", 0)
    b_dl = dl["downloaded_bytes"] + prev.get("downloaded_bytes", 0)
    sec_per_star = dict(download=s_dl / max(n_dl, 1),
                        search=se["real"]["seconds"] * se["real"]["workers"] / max(se["real"]["stars"], 1))
    kept_per_star = dl["kept_total_bytes"] / max(n_star, 1)
    sizes = script_sizes(s)
    bpl = sizes[str(s)] / cut["qlp_files"]
    frac = n_star / cut["qlp_files"]
    proj = []
    for k, b in sizes.items():
        if b is None:
            continue
        sec = int(k)
        nf = b / bpl
        ns = nf * frac
        proj.append(dict(sector=sec, qlp_files_est=nf, stars_est=ns, dl_bytes=ns * era_bytes(sec),
                         kept_bytes=ns * kept_per_star * era_bytes(sec) / 682560))
    P = pd.DataFrame(proj)
    P.to_csv(DATA / "projection_all_sectors.csv", index=False)
    rate = n_dl / max(s_dl, 1)
    proj_tot = dict(sectors=len(P), stars=P.stars_est.sum(), dl_TB=P.dl_bytes.sum() / 1e12, kept_GB=P.kept_bytes.sum() / 1e9,
                    dl_days=P.stars_est.sum() / rate / 86400,
                    search_cpu_h=P.stars_est.sum() * sec_per_star["search"] / 3600)
    rep = dict(starcut=cut, download=dl, search=se, stages=stg, cnn=cnn, crossmatch=xm,
               completeness_total={k: dict(frac=v[0], err=v[1]) for k, v in tot.items()},
               control_level=control_level, flags_cnn_pass=flags, sec_per_star=sec_per_star, projection=proj_tot)
    (DATA / f"s{s:04d}_report.json").write_text(json.dumps(rep, indent=1, default=float))
    xp, xo = xm.get("pass", {}), xm.get("out_of_range", {})
    md = f"""# Pilot: wide single-transit search, Sector {s} (QLP)

Run under `PLAN.md` v1.0 (committed before any download). Counts only: files naming candidate stars stay in
`tess_wide/data/` (outside git) and are not public. Scripts: `s1_starcut.py` … `s7_report.py`.

## Stars

| | Stars |
|---|---|
| QLP light curves in S{s} | {cut['qlp_files']:,} |
| fail T ≤ 12 | {cut['failed_first'].get('tmag', 0):,} |
| fail R★ ≤ 2 R☉ (or no radius) | {cut['failed_first'].get('rad', 0):,} |
| fail TIC disposition (split/artefact/duplicate) | {cut['failed_first'].get('disposition', 0):,} |
| fail contratio ≤ 1 | {cut['failed_first'].get('contratio', 0):,} |
| fail 3000–7500 K | {cut['failed_first'].get('teff', 0):,} |
| fail log g ≥ 3.5 | {cut['failed_first'].get('logg', 0):,} |
| **pass the cut** | **{n_star:,}** |
| downloaded (failed) | {dl['stars']:,} ({dl['failed']}) |
| searched | {R['stars_searched']:,} |

## Events at each stage (real search, and the identical search on 2 − flux)

| Stage | Real | Control |
|---|---|---|
| detections (SNR ≥ 8 and the other `search_any` rules) | {R['detections']:,} on {R['stars_with_detection']:,} stars | {C['detections']:,} on {C['stars_with_detection']:,} stars |
| common-mode bins flagged | {stg['cm_bins']['real']} | {stg['cm_bins']['ctrl']} |
| after common-mode veto | {R['after_common_mode']:,} | {C['after_common_mode']:,} |
| in stars with exactly one event (single) | {R['single_events']:,} | {C['single_events']:,} |
| (events in "repeating" / "multi-event" stars) | {R['star_class'].get('repeating', 0):,} / {R['star_class'].get('multi', 0):,} | {C['star_class'].get('repeating', 0):,} / {C['star_class'].get('multi', 0):,} |
| after size test (Rp ≤ 2 R_J) | {R['after_size']:,} | {C['after_size']:,} |
| CNN pass | {final_real} | {final_ctrl} |
| CNN fail / out of range / error | {cr.get('fail', 0)} / {cr.get('out_of_range', 0)} / {cr.get('error', 0)} | {cc.get('fail', 0)} / {cc.get('out_of_range', 0)} / {cc.get('error', 0)} |
| CNN pass, already known (cross-match) | {xp.get('known', 0)} | – |
| **CNN pass, not known: new single-event candidates** | **{xp.get('new', 0)}** | – |
| out of CNN range, not known (counted separately) | {xo.get('new', 0)} | – |

Cross-match hits among CNN passes, by catalogue: {json.dumps(xp.get('by_catalogue', {}))}. VizieR queries that failed:
{xp.get('query_failed', 0)}. Catalogue snapshots: {xm['run_day']}.

**Control rule (PLAN 4):** final-stage control count {final_ctrl} against real {final_real}:
{"**control-level** — by the rule fixed in PLAN.md, the search is not yet separating signal from systematics in this sector" if control_level else "below half the real count"}.

**Flags among CNN passes (PLAN 5.4; recorded, not used to cut):**

| | Real | Control |
|---|---|---|
| CNN passes | {flags['real']['n']} | {flags['ctrl']['n']} |
| in-event SAP_BKG change > 3σ | {flags['real']['bkg']} | {flags['ctrl']['bkg']} |
| SAP_X or SAP_Y shift > 3σ | {flags['real']['centroid']} | {flags['ctrl']['centroid']} |
| either flag | {flags['real']['either']} | {flags['ctrl']['either']} |
| by camera (1/2/3/4) | {'/'.join(str(flags['real']['camera'].get(c, 0)) for c in (1, 2, 3, 4))} | {'/'.join(str(flags['ctrl']['camera'].get(c, 0)) for c in (1, 2, 3, 4))} |

## Completeness (injection–recovery, {len(tr):,} injections)

Overall: detection {100 * tot['detect'][0]:.0f} ± {100 * tot['detect'][1]:.0f} %; end to end (through the CNN pass,
true Rp ≤ 2 R_J, {int(e2e_mask.sum()):,} injections) {100 * tot['e2e'][0]:.0f} ± {100 * tot['e2e'][1]:.0f} %
({100 * tot['e2e_or'][0]:.0f} % counting CNN out-of-range as kept). Injections falling in data gaps count as missed.

**At detection** (bins of true depth and injected duration; ± is the 1σ Wilson interval):

{table(g_det)}

**End to end:**

{table(g_e2e)}

## Runtime and disk (this MacBook: 18 cores, 128 GB)

| Step | Time | Notes |
|---|---|---|
| TIC rows for {cut['qlp_files']:,} stars (MAST TAP, 300 ids per query) | {cut['tap_seconds'] / 60:.0f} min | the TAP server slowed for ~10 min mid-run |
| download + extract {n_dl:,} stars (two runs, started on a partial star list; 200 more in a timing test) | {s_dl / 60:.0f} min | {b_dl / 1e9:.1f} GB in; FITS deleted after extraction; 16 parallel requests |
| search (real + control) | {se['real']['seconds'] / 60:.1f} min on {se['real']['workers']} cores | {sec_per_star['search'] * 1000:.0f} ms CPU per star |
| injections ({se['inj']['stars']:,}) | {se['inj']['seconds'] / 60:.1f} min | |
| CNN (real / control / injections) | {cnn['real']['seconds']:.0f} / {cnn['ctrl']['seconds']:.0f} / {cnn['inj']['seconds']:.0f} s | CPU, weights loaded once |
| disk kept (compressed arrays) | {dl['kept_total_bytes'] / 1e9:.1f} GB | {kept_per_star / 1e3:.0f} kB per star |

## Projection for all QLP sectors (S1–S{LAST_QLP}) on this MacBook

Assumes every sector passes the same fraction of its QLP stars as S{s} ({100 * frac:.1f} %), file sizes by cadence era
as measured on MAST, and the pilot's download rate ({rate:.1f} stars/s) and per-star search cost.
- stars (star-sectors): {proj_tot['stars'] / 1e6:.2f} million over {proj_tot['sectors']} sectors (MAST has no QLP bulk
  script for {', '.join('S' + k for k, v in sizes.items() if v is None) or 'none'}; they are left out);
- download: {proj_tot['dl_TB']:.2f} TB in, {proj_tot['kept_GB']:.0f} GB kept, about {proj_tot['dl_days']:.1f} days at the pilot's rate;
- search: {proj_tot['search_cpu_h']:.0f} CPU-hours ({proj_tot['search_cpu_h'] / 16:.1f} h on 16 cores), plus injections;
- the CNN cost is per surviving event and negligible next to downloads.

The download dominates. The pass fraction is the weakest assumption: it depends on galactic latitude and on how
deep QLP went in each sector.
"""
    (HERE / f"PILOT_S{s}.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
