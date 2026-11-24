"""Synthetic injection tests for the eclipse and pulsation clocks.

Run: .venv/bin/python tests/test_synthetic.py   (needs results/template_S58.pkl, prepared_main.pkl)
1. Eclipse times: S58 cadences, template shifted by a known delay d(t) = 40 s sin(2 pi t / 17 d),
   plus the S58 pulsation model and white noise at the S58 residual rms. Recovered O-C must follow
   d(t): rms of (recovered - injected) within 2x the median error, and the bootstrap errors must
   not underestimate that scatter by more than a factor of 1.5.
2. Pulsation delays: the S58 pulsation model evaluated at t - tau(t), tau = -60 s cos(2 pi phi_orb),
   plus noise. The 0.3-d common-tau fit must recover the cos term as -60 +/- 3 sigma, with the sign
   the right way round (a < 0).
"""
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import hip  # noqa: E402
import s2_pulsation as s2  # noqa: E402

rng = np.random.default_rng(42)
TPL = pickle.load(open(hip.OUTDIR / "template_S58.pkl", "rb"))
P = pickle.load(open(hip.OUTDIR / "prepared_main.pkl", "rb"))
d58 = P[58]
t = d58["t"]
FN = list(hip.FREQS.values())
coef = d58["coef"]
noise = 0.0009


def test_eclipses():
    delay = lambda tt: 40.0 / 86400 * np.sin(2 * np.pi * tt / 17.0)
    flux = TPL(hip.phase(t - delay(t))) + hip.sines_model(t, FN, coef) + rng.normal(0, noise, t.size)
    y = flux - hip.sines_model(t, FN, coef)
    e = np.full_like(t, noise)
    rec, inj, err = [], [], []
    for E in range(int((t.min() - hip.T0_LEE) / hip.P_LEE), int((t.max() - hip.T0_LEE) / hip.P_LEE) + 1):
        for which, off in (("p", 0.0), ("s", 0.5)):
            tc = hip.T0_LEE + (E + off + TPL.shift[which]) * hip.P_LEE
            win = np.abs(t - tc) < hip.FIT_HALF * hip.P_LEE
            core = np.abs(t - tc) < hip.ECL_HALF * hip.P_LEE
            if core.sum() < 0.8 * 2 * hip.ECL_HALF * hip.P_LEE * 86400 / 120:
                continue
            prof = TPL.eclipse_profile(which)
            best, _ = hip.fit_eclipse(t[win], y[win], e[win], tc, prof, 120.0)
            sig = hip.block_bootstrap(t[win], y[win], e[win], tc, prof, 120.0, best, nboot=60)
            rec.append(best[0] * 86400); inj.append(delay(tc) * 86400); err.append(sig * 86400)
    rec, inj, err = map(np.array, (rec, inj, err))
    scat = np.std(rec - inj)
    print(f"eclipse test: {len(rec)} eclipses, injected amplitude 40 s, rms(rec-inj) {scat:.2f} s, "
          f"median boot err {np.median(err):.2f} s, mean offset {np.mean(rec-inj):+.2f} s")
    assert len(rec) >= 15
    assert scat < 2 * np.median(err) + 1.0
    assert np.median(err) > scat / 1.5
    assert abs(np.mean(rec - inj)) < 3 * scat / np.sqrt(len(rec)) + 0.5


def test_pulsation_sign():
    A = 60.0
    tau = -A / 86400 * np.cos(2 * np.pi * (hip.phase(t) % 1.0))
    r = hip.sines_model(t - tau, FN, coef) + rng.normal(0, noise, t.size)
    oot = ~hip.in_eclipse(t)
    rows = []
    for idx in s2.segments(t, length=0.3):
        idx = idx[oot[idx]]
        if len(idx) < 100:
            continue
        tt, st, _ = s2.common_tau(t[idx], r[idx], coef)
        rows.append([hip.phase(np.mean(t[idx])) % 1.0, tt * 86400, st * 86400])
    ph, y, e = np.array(rows).T
    X = np.column_stack([np.ones_like(ph), np.cos(2 * np.pi * ph), np.sin(2 * np.pi * ph)])
    W = 1 / e ** 2
    cov = np.linalg.inv((X * W[:, None]).T @ X)
    c = cov @ (X * W[:, None]).T @ y
    print(f"pulsation test: injected a = -{A:.0f} s, recovered a = {c[1]:+.1f} +/- {np.sqrt(cov[1,1]):.1f} s, "
          f"b = {c[2]:+.1f} s")
    assert abs(c[1] + A) < 3 * np.sqrt(cov[1, 1]) + 5
    assert c[1] < 0


if __name__ == "__main__":
    test_eclipses()
    test_pulsation_sign()
    print("synthetic tests passed")
