"""Independent timing check: (a) Kwee-van Woerden-style mirror symmetry minimisation; (b) season template cross-fit
(stack the season's eclipses of the same kind after aligning on (a), then fit each event with shift + depth scale + linear baseline).
Uses own TESScut r=1.5 photometry and also SPOC PDCSAP (S43-45) / QLP (S71-72)."""
import numpy as np, pandas as pd, glob, sys
from scipy.optimize import minimize_scalar, least_squares
from scipy.interpolate import interp1d
sys.path.insert(0, "scripts")
ev = pd.read_csv("data/eclipse_times.csv"); ev = ev[ev.src.str.startswith("TESS")]
def load(sec):
    d = pd.read_csv(f"data/phot/tess_s{sec}.csv"); d = d[d.q == 0]
    return d.time.values, (d["f_r1.5"]/np.nanmedian(d["f_r1.5"])).values
def norm_local(t, f, tc, half):
    out = (np.abs(t - tc) > 1.3*half) & (np.abs(t - tc) < 1.3*half + 0.5)
    c = np.polyfit(t[out] - tc, f[out], 1)
    return f/np.polyval(c, t - tc)
def kvw(t, f, tc0, half):
    def cost(tc):
        m = np.abs(t - tc) < half
        if m.sum() < 10: return 1e9
        fi = interp1d(t, f, bounds_error=False, fill_value=np.nan)
        refl = fi(2*tc - t[m]); ok = np.isfinite(refl)
        return np.mean((f[m][ok] - refl[ok])**2)
    r = minimize_scalar(cost, bounds=(tc0 - 0.02, tc0 + 0.02), method="bounded", options={"xatol": 1e-6})
    return r.x
out = []
for _, e in ev.iterrows():
    t, f = load(int(e.sector))
    half = 2.2*e.w                     # ~ eclipse half-duration
    fn = norm_local(t, f, e.tc, half)
    tk = kvw(t, fn, e.tc, half*0.9)
    out.append(dict(kind=e.kind, cycle=e.cycle, sector=e.sector, tc_gg=e.tc, tc_kvw=tk, diff_min=(tk - e.tc)*1440, coverage=e.coverage))
O = pd.DataFrame(out)
print(O.to_string(index=False, float_format=lambda x: f"{x:.5f}"))
O.to_csv("data/timing_check_kvw.csv", index=False)
