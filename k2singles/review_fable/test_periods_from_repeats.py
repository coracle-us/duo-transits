"""periods_from_repeats: one spurious 'robust' repeat wipes out every candidate period, including the
right one.  Also: a period is only ever tried as (T_nearest_repeat - epoch)/n, so if the nearest repeat
is the spurious one, the true period is never even proposed.

Run: .venv/Scripts/python.exe k2singles/review_fable/test_periods_from_repeats.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "singles" / "scripts"))
sys.path.insert(0, str(ROOT / "scripts"))
import pilot  # noqa: E402

D = 3.0 / 24
# a grid with NO flat windows at all (so excluded_by_flat never fires) -> isolates the ephemeris logic
g = pd.DataFrame({"T_pred": np.arange(0, 3300, D / 4), "verdict": "unconstrained"})
epoch = 179.9
P_true = 140.494
real = [dict(T=epoch + 7 * P_true), dict(T=epoch + 17 * P_true), dict(T=epoch + 22 * P_true)]

print("real repeats only       ->", pilot.periods_from_repeats(g, epoch, real, D)[:6])
spur = dict(T=epoch + 11.37 * P_true)  # a systematic somewhere between the real ones
print("real + 1 spurious repeat->", pilot.periods_from_repeats(g, epoch, real + [spur], D)[:6])
spur_near = dict(T=epoch + 3.3 * P_true)  # a systematic CLOSER to the event than any real repeat
print("real + spurious nearest ->", pilot.periods_from_repeats(g, epoch, real + [spur_near], D)[:6])

# and the cheap alternative: drop one repeat at a time
def leave_one_out(g, epoch, reps, D):
    out = {}
    for i in range(len(reps)):
        sub = reps[:i] + reps[i + 1:]
        for n, P in pilot.periods_from_repeats(g, epoch, sub, D):
            out.setdefault(round(P, 2), []).append(i)
    return out

print("\nleave-one-out with the spurious repeat present:", leave_one_out(g, epoch, real + [spur], D))
