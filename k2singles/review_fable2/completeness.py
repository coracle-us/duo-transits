"""Completeness spot-check of run_full2: for a sample of max_support <= 1 events (deepest / highest TESS SNR per transit),
run an independent sliding-box dip search over the SAME combined TESS+K2 light curves (k2lib.load_combined) with
relaxed criteria (depth 0.4-2.5 x catalogue depth, SNR >= 7, duration-matched box) and compare with the pipeline's
repeats. Also re-examine the catalogued K2 event itself. Writes completeness.txt."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts")); sys.path.insert(0, str(ROOT / "singles" / "scripts")); sys.path.insert(0, str(ROOT / "k2singles" / "scripts"))
import k2lib  # noqa: E402

summ = pd.read_csv(ROOT / "k2singles" / "results" / "full2" / "summary.csv")
cat = pd.read_csv(ROOT / "k2singles" / "data" / "k2_singles.csv")
cat["event_idx"] = cat.groupby("EPIC").cumcount()
m = summ.merge(cat[["EPIC", "event_idx", "tess_snr_per_transit", "epoch_btjd", "campaign"]], on=["EPIC", "event_idx"])
low = m[m.max_support <= 1].copy()
pick = pd.concat([low.sort_values("depth_ppm", ascending=False).head(3),
                  low.sort_values("tess_snr_per_transit", ascending=False).head(3),
                  low[low.max_support == 1].sort_values("tess_snr_per_transit", ascending=False).head(3),
                  low[low.n_repeats_k2 == 0].sample(2, random_state=3)]).drop_duplicates(["EPIC", "event_idx"])
if len(sys.argv) > 1:
    pick = m[m.EPIC.isin([int(x) for x in sys.argv[1:]])]
out = open(HERE / "completeness.txt", "w")


def say(*a):
    print(*a, flush=True); print(*a, file=out)


def find_dips(lc, D, delta, snr_min=7.0, lo=0.4, hi=2.5):
    t, f, sec = lc.time.to_numpy(), lc.flux.to_numpy(), lc.sector.to_numpy()
    hits = []
    for T in np.arange(t.min(), t.max(), D / 6):
        inn = np.abs(t - T) < 0.35 * D
        oot = (np.abs(t - T) > 0.5 * D + 0.05) & (np.abs(t - T) < 0.5 * D + max(0.8, D))
        if inn.sum() < 4 or (oot & (t < T)).sum() < 8 or (oot & (t > T)).sum() < 8:
            continue
        c = np.polyfit(t[oot] - T, f[oot], 1)
        r = f / np.polyval(c, t - T) - 1
        sig = 1.4826 * np.median(np.abs(r[oot] - np.median(r[oot]))) + 1e-7
        d = -np.mean(r[inn]); s = sig / np.sqrt(inn.sum())
        # red-noise guard: scatter of oot means over boxes of the same size
        nb = max(1, int(round(0.7 * D / np.median(np.diff(t[oot])))))
        ro = r[oot]; nbl = len(ro) // nb
        s_red = np.std(ro[:nbl * nb].reshape(nbl, nb).mean(axis=1)) if nbl >= 3 else s
        s = max(s, s_red)
        if lo * delta <= d <= hi * delta and d / s >= snr_min:
            hits.append((T, d, d / s, int(sec[min(np.searchsorted(t, T), len(t) - 1)])))
    hits.sort()
    cl = []
    for h in hits:
        if cl and h[0] - cl[-1][-1][0] < D:
            cl[-1].append(h)
        else:
            cl.append([h])
    return [max(c, key=lambda x: x[2]) for c in cl]


say(f"Sample: {len(pick)} events (of {len(low)} with max_support <= 1)")
for _, r in pick.iterrows():
    key = f"{int(r.EPIC)}_{int(r.event_idx)}"
    j = json.loads((ROOT / "k2singles" / "results" / "full2" / f"{key}.json").read_text())
    D, delta, epoch = r.dur_h / 24, r.depth_ppm * 1e-6, float(r.epoch_btjd)
    say(f"\n=== EPIC {key} TIC {int(r.TIC)}: depth {delta*1e3:.1f} ppt, D {r.dur_h:.1f} h, K2 C{r.campaign} epoch {epoch:.3f}, "
        f"TESS SNR/transit {r.tess_snr_per_transit:.0f}, max_support {int(r.max_support)}, pipeline repeats {len(j['repeats'])}, "
        f"comment: {r.comment}")
    say(f"    sectors: K2 {j['k2_campaigns_loaded']}  TESS {j['tess_sectors_loaded']}; excluded windows {j['n_excluded_windows']}/{j['n_windows']}, frac_logP_excl {j['frac_logP_excluded']:.2f}")
    try:
        lc, k2 = k2lib.load_combined(int(r.TIC), int(r.EPIC))
    except Exception as e:  # noqa: BLE001
        say("    load failed:", e); continue
    dips = find_dips(lc, D, delta)
    pipe = [x["T"] for x in j["repeats"]]
    say(f"    independent search: {len(dips)} dips (depth 0.4-2.5x, SNR>=7):")
    for T, d, snr, sec in dips:
        tag = "catalogue event" if abs(T - epoch) < D else ("= pipeline repeat" if any(abs(T - p) < D for p in pipe) else "** NOT in pipeline **")
        say(f"      T={T:9.3f} sector {sec:5d} depth {d*1e3:7.2f} ppt SNR {snr:5.1f}  {tag}")
    for p in pipe:
        if not any(abs(T - p) < D for T, *_ in dips):
            say(f"      pipeline repeat at {p} not recovered by the independent search (depth {[x['depth_ppm'] for x in j['repeats'] if x['T']==p]})")
    # catalogue event check
    t, f = lc.time.to_numpy(), lc.flux.to_numpy()
    inn = np.abs(t - epoch) < 0.35 * D; oot = (np.abs(t - epoch) > 0.5 * D + 0.05) & (np.abs(t - epoch) < 0.5 * D + max(0.8, D))
    if inn.sum() > 3 and oot.sum() > 10:
        c = np.polyfit(t[oot] - epoch, f[oot], 1); rr = f / np.polyval(c, t - epoch) - 1
        say(f"    catalogue event at {epoch:.3f}: measured box depth {-np.mean(rr[inn])*1e3:.2f} ppt (catalogue {delta*1e3:.1f}) in K2 data")
    else:
        say(f"    catalogue event at {epoch:.3f}: K2 data not loaded around the event (in={inn.sum()}, oot={oot.sum()})")
out.close()
