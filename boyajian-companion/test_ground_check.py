#!/usr/bin/env python3
"""Synthetic check of ground_check.test_window on AAVSO-format data (no network).

Observer A: time series, 90 points a night, 0.008 mag nightly zero-point offsets, no transit.
  The per-point statistic calls many of its nights 'excluded' or 'consistent' (false dips);
  the nightly-mean statistic must leave them 'untested'.
Observer B: precise (0.001 mag nightly offsets) with a 1.1% transit injected in window 3.
  Its transit night must come out 'consistent', its other nights 'excluded'.

  python test_ground_check.py
"""
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import ground_check as gc

rng = np.random.default_rng(3)
W3 = gc.WINDOWS[0]
TC_TRUE = gc.T0 + 2 * 1302.0  # inside window 3


def observer(code, nights, offset_sig, point_sig, transit=None):
    rows = []
    for night in nights:
        t = night + 0.15 + np.arange(90) * 2 / 1440  # 3 h from ~03:36 UTC
        m = 11.9 + rng.normal(0, offset_sig) + rng.normal(0, point_sig, t.size)
        if transit is not None:
            m[np.abs(t - transit) < 0.4] += gc.DEPTH_MAG
        rows += [dict(JD=x, Magnitude=y, Uncertainty=point_sig, Band="V", **{"Observer Code": code, "Fainter Than": 0})
                 for x, y in zip(t, m)]
    return rows


def main():
    nights = np.arange(np.floor(TC_TRUE) - 20, np.floor(TC_TRUE) + 20)
    nights_b = nights[::2]
    tc_b = nights_b[np.argmin(np.abs(nights_b + 0.25 - TC_TRUE))] + 0.25  # B observes the transit core
    rows = observer("AAA", nights, 0.008, 0.01) + observer("BBB", nights_b, 0.001, 0.005, transit=tc_b)
    path = Path(tempfile.mkdtemp()) / "aavso_fake.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    data = gc.load_aavso(str(path))

    a = gc.test_window(data[data.group.str.contains("AAA")], *W3)
    on_a = a[a.groups != ""]
    print("observer A only:", on_a.status.value_counts().to_dict(), "| per point:", on_a.status_per_point.value_counts().to_dict())
    assert (on_a.status_per_point != "untested").mean() > 0.3, "per-point statistic gives verdicts on pure nightly noise"
    assert (on_a.status == "excluded").mean() < 0.05, "nightly-mean statistic must not exclude on 0.008 mag nightly noise"

    b = gc.test_window(data[data.group.str.contains("BBB")], *W3)
    at = b.iloc[(b.tc_bjd - tc_b).abs().argmin()]
    other = b[(b.groups != "") & ((b.tc_bjd - tc_b).abs() > 0.8)]
    print("observer B at the transit:", at.status, round(at.dchi2, 1), "| other covered tc:", other.status.value_counts().to_dict())
    assert at.status == "consistent", at
    assert (other.status == "excluded").mean() > 0.9
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
