"""Collect final numbers: event table, ephemerides, secondary phase / rough e-omega, future eclipse times, results.csv."""
import sys, os, numpy as np, pandas as pd
from scipy.optimize import brentq
from astropy.time import Time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib2

H = lib2.HERE
TIC = {'246194159': 301328232, '212012030': 14394716, '212694013': 178941241}
allev = []
for e in TIC:
    ev = pd.read_csv(os.path.join(H, f'events_combined_{e}.csv'))
    ev.insert(0, 'epic', e)
    allev.append(ev)
allev = pd.concat(allev)
allev.to_csv(os.path.join(H, 'events_all.csv'), index=False)
print(allev[['epic', 'event', 'type', 't0', 'err', 'depth', 'T14_h', 'T12_h', 'sources']].to_string())

# ---- EPIC 246194159 ephemeris
c = pd.read_csv(os.path.join(H, 'period_candidates_246194159.csv'))
a = c[c.status == 'ALLOWED'].iloc[0]
P, sP = a.P, a.sigP
T0p, sT0p, T0s, sT0s = a.T0_narrow, a.sT0_narrow, a.T0_wide, a.sT0_wide
# put the secondary reference epoch in the same cycle (after primary)
phi = ((T0s - T0p) / P) % 1
sphi = np.hypot(sT0p, sT0s) / P
print(f'\n246194159: P = {P:.6f} +- {sP:.6f} d (chi2_timing={a.chi2_timing:.2f}); T0_prim(narrow) = {T0p:.5f} +- {sT0p:.5f}; '
      f'T0_sec(wide) = {T0s:.5f} +- {sT0s:.5f}; secondary phase = {phi:.5f} +- {sphi:.5f}')
ev = pd.read_csv(os.path.join(H, 'events_combined_246194159.csv'))
Dp = ev[ev.type == 'narrow'].T14_h.mean(); Ds = ev[ev.type == 'wide'].T14_h.mean()
esinw = (Ds - Dp) / (Ds + Dp)


def sterne(ecw, esw):
    e = np.hypot(ecw, esw)
    return 0.5 + (np.arctan(ecw / np.sqrt(1 - e ** 2)) + ecw * np.sqrt(1 - e ** 2) / (1 - esw ** 2)) / np.pi


ecosw = brentq(lambda x: sterne(x, esinw) - phi, 0.0, np.sqrt(1 - esinw ** 2) - 1e-6)
print(f'  rough (i=90 deg, non-grazing approx.): e sin w ~ {esinw:.3f} (from durations {Dp:.2f}h / {Ds:.2f}h), '
      f'e cos w ~ {ecosw:.3f}, e ~ {np.hypot(ecosw, esinw):.2f}')
now = 4325.5  # BTJD of 2026-10-12
k = np.ceil((now - T0p) / P)
print('  next eclipses (BTJD, UTC):')
for j in range(3):
    for lab, T in (('primary(narrow)', T0p), ('secondary(wide)', T0p + phi * P)):
        kk = np.ceil((now - T) / P) + j
        t = T + kk * P
        st = np.hypot(sT0p, kk * sP)
        print(f'    {lab:16s} {t:.4f} +- {st*1440:.1f} min  {Time(t + 2457000, format="jd", scale="tdb").utc.iso[:16]}')

# ---- results table
rows = []
rows.append(dict(epic='246194159', tic=301328232,
                 events='747.3850(narrow,K2C12);810.9061(wide,K2C12);3227.9968(wide,TESS S70);3813.9579(wide,TESS S92);3823.6824(narrow,TESS S92)',
                 n_events=5, allowed_periods_d=f'{P:.5f}+-{max(sP,1e-5):.5f}', excluded='all P=73.245/j for j>=2 (K2 C12, C19, TESS S42/S70/S92)',
                 secondary_phase=f'{phi:.4f}', depth_primary=0.375, dur_primary_h=round(Dp, 2), depth_secondary=0.299, dur_secondary_h=round(Ds, 2),
                 twin_2x='no (two distinct eclipse types, both seen)', nature='eccentric eclipsing binary (grazing/V-shaped eclipses, e>~0.6)',
                 confidence='high'))
c2 = pd.read_csv(os.path.join(H, 'period_candidates_212012030.csv'))
al2 = c2[c2.status == 'ALLOWED']
rows.append(dict(epic='212012030', tic=14394716,
                 events='174.4713(K2C5);1276.0801(K2C18);3258.9766(TESS S71)', n_events=3,
                 allowed_periods_d=';'.join(f'{p:.4f}' for p in al2.P), excluded='P=220.32/k for k>=6 (K2 C5/C18 coverage)',
                 secondary_phase='none detected (100% phase coverage for P<=73.4 d, no dip >~0.1%)',
                 depth_primary=0.202, dur_primary_h=3.89, depth_secondary=np.nan, dur_secondary_h=np.nan,
                 twin_2x='possible (all 3 events identical depth/shape); 2P = ' + ';'.join(f'{2*p:.3f}' for p in al2.P),
                 nature='deep (20%) V/U-shaped eclipse on K dwarf; likely EB or grazing large companion', confidence='high (alias set); unresolved'))
c3 = pd.read_csv(os.path.join(H, 'period_candidates_212694013.csv'))
rows.append(dict(epic='212694013', tic=178941241,
                 events='242.158(K2C6);1215.947(K2C17)', n_events=2,
                 allowed_periods_d='973.789/N for N in {1,2,3,4,5,6,7,9,10,11,12,14,15,17}; N=14,17 disfavoured',
                 excluded='N>=18 (K2 C6/C17); N=13 (TESS S46, 3 reductions); N=8,16 likely excluded (S46 orbit start)',
                 secondary_phase='none detected', depth_primary=0.0056, dur_primary_h=40.5, depth_secondary=np.nan, dur_secondary_h=np.nan,
                 twin_2x='implausible (flat-bottomed 0.55% transit of a giant; equal-depth occultation impossible)',
                 nature='transit/eclipse of a ~0.3 Rsun body across a 4.4 Rsun giant',
                 confidence='high (2 K2 events); period unresolved'))
pd.DataFrame(rows).to_csv(os.path.join(H, 'results.csv'), index=False)
print(c3[c3.status == 'ALLOWED'][['hypothesis', 'P', 'sigP']].to_string())
print(al2[['hypothesis', 'P', 'sigP', 'T0_A', 'sT0_A']].to_string())
