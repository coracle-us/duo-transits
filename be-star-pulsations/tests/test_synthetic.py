"""Synthetic checks for the pilot pipeline. No network. Run: python tests/test_synthetic.py"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from metrics import prepare, sector_metrics          # noqa: E402
from wise import DiskCurve                           # noqa: E402
import stats4                                        # noqa: E402

FAIL = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    if not ok:
        FAIL.append(name)


def fake_sector(rng, amps=(5.0, 3.0), freqs=(1.1, 2.2), white=1.0, red=2.0, tau=0.3):
    dt = 30 / 1440
    t = np.arange(0, 27.0, dt)
    t = t[(t < 13.0) | (t > 14.0)]                    # mid-sector downlink gap
    y = rng.normal(0, white, len(t))
    phi = np.exp(-dt / tau)
    r = np.zeros(len(t))
    e = rng.normal(0, red * np.sqrt(1 - phi ** 2), len(t))
    for i in range(1, len(t)):
        r[i] = phi * r[i - 1] + e[i]
    y += r
    for a, f in zip(amps, freqs):
        y += a * np.sin(2 * np.pi * f * t + rng.uniform(0, 2 * np.pi))
    return t + 1500.0, 1.0 + y / 1e3


def test_metrics():
    """Brief: A1 back within 5%, f1 within 0.01 c/d. With 2 ppt AR(1) red noise the
    amplitude of a 5 ppt sinusoid itself scatters by ~4% per 27-d sector (a direct sine
    fit to the raw fake data scatters as much), so a single-sector 5% check fails ~20% of
    the time by noise alone. We test the pipeline's own error instead: every sector's A1
    within 5% of a sine fit to the raw fake light curve, and the mean over 20 sectors
    within 5% of the true 5 ppt. The per-sector rate within 5% of truth is reported."""
    rng = np.random.default_rng(1)
    worst_pipe, worst_f, rel = 0, 0, []
    from metrics import sinefit
    for k in range(20):
        t, f = fake_sector(rng)
        tb, yb = prepare(t, f)
        m, peaks, spec = sector_metrics(tb, yb)
        _, a_raw = sinefit(t, (f - 1) * 1e3, [1.1])
        worst_pipe = max(worst_pipe, abs(m["A1"] / a_raw[0] - 1))
        worst_f = max(worst_f, abs(m["f1"] - 1.1))
        rel.append(m["A1"] / 5.0)
    rel = np.array(rel)
    check("A1 pipeline error < 5% vs raw sine fit (20 sectors)", worst_pipe < 0.05,
          f"worst {worst_pipe:.3f}")
    check("mean A1 within 5% of truth", abs(rel.mean() - 1) < 0.05,
          f"mean {rel.mean():.3f}, scatter {rel.std():.3f}, "
          f"{np.mean(np.abs(rel - 1) < 0.05):.0%} of single sectors within 5%")
    check("f1 within 0.01 c/d (20 sectors)", worst_f < 0.01, f"worst {worst_f:.4f}")
    check("second peak found near 2.2", any(abs(p["freq"] - 2.2) < 0.01 for p in peaks),
          f"{[round(p['freq'], 3) for p in peaks][:4]}")
    check("log spectrum shape", spec.shape == (600,))
    # band rms of a lone 5 ppt sinusoid in B2 should be ~5/sqrt(2)
    t, f = fake_sector(rng, amps=(5.0,), freqs=(1.1,), white=0.3, red=0.0)
    m, _, _ = sector_metrics(*prepare(t, f))
    check("B2 rms of a pure sinusoid ~ A/sqrt2", abs(m["B2"] / (5 / np.sqrt(2)) - 1) < 0.1,
          f"{m['B2']:.3f} vs {5/np.sqrt(2):.3f}")


VISITS = np.r_[55210.0, 55390.0, np.arange(56660.0, 60520.0, 182.6)]
SECTOR_MID = 58324.0 + 13.7 + 27.3 * np.arange(80)


def fake_disk(rng, trend=0.0):
    w1 = np.full(len(VISITS), rng.uniform(9, 12))
    for _ in range(rng.integers(1, 4)):
        c, w, a = rng.uniform(55500, 60500), rng.uniform(200, 900), rng.uniform(0.1, 0.6)
        w1 -= a * np.exp(-0.5 * ((VISITS - c) / w) ** 2)
    w1 -= trend * (VISITS - 56660.0) / 1000.0
    w1 += rng.normal(0, 0.02, len(VISITS))
    return DiskCurve(VISITS, w1)


def evaluate(d, t):
    o = d.at(t)
    return o["disk_excess"], o["phase_label"], o["wise_usable"]


def fake_sample(rng, n_star=200, beta=0.0, drift=0.0, trend=0.0, scatter=0.15):
    disks, sid, x, lab, use, t, sec = [], [], [], [], [], [], []
    for i in range(n_star):
        d = fake_disk(rng, trend)
        disks.append(d)
        s = np.sort(rng.choice(80, size=rng.integers(2, 7), replace=False))
        xx, ll, uu = evaluate(d, SECTOR_MID[s])
        sid += [i] * len(s); x += list(xx); lab += list(ll); use += list(uu)
        t += list(SECTOR_MID[s]); sec += list(s)
    sid, x, use, t, sec = map(np.array, (sid, x, use, t, sec))
    lab = np.array(lab, dtype=object)
    a0 = rng.normal(0, 1, n_star)[sid]
    y = a0 + beta * x + drift * sec + rng.normal(0, scatter, len(x))
    return sid, x, y, lab, use, t, disks, {i: i for i in range(n_star)}


def test_stage4():
    rng = np.random.default_rng(2)
    # Case A: coupling
    sid, x, y, lab, use, t, disks, dos = fake_sample(rng, beta=0.8)
    r = stats4.run_test(sid, x, y, lab, use, t, disks, dos, evaluate, n_null=2000, n_boot=500)
    check("case A beta ~ 0.8", abs(r["beta"] - 0.8) < 0.15,
          f"beta {r['beta']:.3f} [{r['beta_lo']:.3f},{r['beta_hi']:.3f}] n_stars {r['n_stars']}")
    check("case A p < 0.01 both nulls", r["p_perm"] < 0.01 and r["p_swap"] < 0.01,
          f"p_perm {r['p_perm']:.4f} p_swap {r['p_swap']:.4f}")
    # Case B: no coupling, 100 repeats
    pp, ps = [], []
    for k in range(100):
        sid, x, y, lab, use, t, disks, dos = fake_sample(rng)
        r = stats4.run_test(sid, x, y, lab, use, t, disks, dos, evaluate, n_null=400,
                            n_boot=10, seed=k)
        pp.append(r["p_perm"]); ps.append(r["p_swap"])
    rp, rs = np.mean(np.array(pp) < 0.05), np.mean(np.array(ps) < 0.05)
    check("case B false-positive rate ~5% (perm)", 0.0 <= rp <= 0.12, f"rate {rp:.2f}")
    check("case B false-positive rate ~5% (swap)", 0.0 <= rs <= 0.12, f"rate {rs:.2f}")
    # Case C: amplitude drifts with sector number; disks share a common time trend
    pp, ps = [], []
    for k in range(20):
        sid, x, y, lab, use, t, disks, dos = fake_sample(rng, drift=0.01, trend=0.08)
        r = stats4.run_test(sid, x, y, lab, use, t, disks, dos, evaluate, n_null=400,
                            n_boot=10, seed=100 + k)
        pp.append(r["p_perm"]); ps.append(r["p_swap"])
    rp, rs = np.mean(np.array(pp) < 0.05), np.mean(np.array(ps) < 0.05)
    print(f"      case C: within-star perm null p<0.05 in {rp:.2f} of 20 (expected high)")
    check("case C cross-star null not significant", rs <= 0.15, f"rate {rs:.2f}")


if __name__ == "__main__":
    test_metrics()
    test_stage4()
    if FAIL:
        print("FAILED:", FAIL)
        sys.exit(1)
    print("ALL CHECKS PASSED")
