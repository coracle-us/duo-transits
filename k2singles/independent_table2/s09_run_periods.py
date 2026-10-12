"""Driver: combine per-source event fits, build candidate periods, run coverage tests. Usage: s09_run_periods.py EPIC"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib2
from s08_periods import analyse

EPIC = sys.argv[1]
fits_df = pd.read_csv(os.path.join(lib2.HERE, f'event_fits_{EPIC}.csv'))
# combine sources per event: unweighted mean; error = max(median per-source error, scatter between sources)
ev = []
for lab, g in fits_df.groupby('event', sort=False):
    t0 = g.t0.mean()
    err = max(g.t0_err.median(), g.t0.std(ddof=1) if len(g) > 1 else 0)
    ev.append(dict(event=lab, type=g.type.iloc[0], t0=t0, err=err, depth=g.depth.median(), T14_h=g.T14_h.median(),
                   T12_h=g.T12_h.median(), sources=';'.join(g.source)))
ev = pd.DataFrame(ev).sort_values('t0')
ev.to_csv(os.path.join(lib2.HERE, f'events_combined_{EPIC}.csv'), index=False)
print(ev.to_string())


def linfit(times, errs, cycles):
    A = np.vstack([np.ones_like(cycles), cycles]).T.astype(float)
    W = np.diag(1 / np.asarray(errs) ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    p = cov @ A.T @ W @ np.asarray(times)
    res = np.asarray(times) - A @ p
    chi2 = np.sum(res ** 2 / np.asarray(errs) ** 2)
    return p, cov, res, chi2


def cls_from(evsub, name):
    return dict(name=name, T14=evsub.T14_h.median() / 24, T12=evsub.T12_h.median() / 24, depth=evsub.depth.median())


TOL = 0.005  # d, extra timing tolerance allowed for consistency (ETVs / systematics)
PMIN = 1.0
out_rows, summ = [], []

if EPIC == '246194159':
    sources = ['K2-c12-PDCSAP-LC', 'K2-c12-K2SFF', 'K2-c19-SAP-LC', 'K2-c19-K2SFF',
               'TESS-S42-SPOC120', 'TESS-S42-QLP', 'TESS-S42-CUT', 'TESS-S70-TESSSPOC', 'TESS-S70-QLP', 'TESS-S70-CUT',
               'TESS-S92-QLP', 'TESS-S92-CUT']
    nar = ev[ev.type == 'narrow']
    wid = ev[ev.type == 'wide']
    classes = [cls_from(nar, 'narrow'), cls_from(wid, 'wide')]
    D = nar.t0.iloc[-1] - nar.t0.iloc[0]
    cands = []
    for N in range(1, int(D / PMIN) + 1):
        P = D / N
        # wide events must share one phase
        ph = ((wid.t0.values - wid.t0.values[0]) / P)
        dev = np.abs(ph - np.round(ph)) * P
        if np.any(dev > 5 * np.hypot(wid.err.values, wid.err.values[0]) + TOL):
            continue
        # joint fit: common P, separate T0 per class
        tt = np.r_[nar.t0.values, wid.t0.values]
        ee = np.r_[nar.err.values, wid.err.values]
        cn = np.round((nar.t0.values - nar.t0.values[0]) / P)
        cw = np.round((wid.t0.values - wid.t0.values[0]) / P)
        A = np.zeros((len(tt), 3))
        A[:len(cn), 0] = 1; A[len(cn):, 1] = 1
        A[:, 2] = np.r_[cn, cw]
        W = np.diag(1 / ee ** 2)
        cov = np.linalg.inv(A.T @ W @ A)
        p = cov @ A.T @ W @ tt
        res = tt - A @ p
        chi2 = np.sum(res ** 2 / ee ** 2)
        cands.append(dict(P=p[2], sigP=np.sqrt(cov[2, 2]), tag=f'N={N}', T0_narrow=p[0], sT0_narrow=np.sqrt(cov[0, 0]),
                          T0_wide=p[1], sT0_wide=np.sqrt(cov[1, 1]), chi2=chi2, res=res))
    print(f'{len(cands)} candidate periods consistent with all event times')
    for c in cands[:5]:
        print('  ', c['tag'], c['P'], c['sigP'], 'chi2', c['chi2'], 'res(min)', np.round(c['res'] * 1440, 2))
    Bw = 0.5
    known = ev.t0.values
elif EPIC in ('212012030', '212694013'):
    if EPIC == '212012030':
        sources = ['K2-c05-PDCSAP-LC', 'K2-c05-K2SFF', 'K2-c18-PDCSAP-LC', 'K2-c18-K2SFF',
                   'TESS-S71-SPOC120', 'TESS-S71-QLP', 'TESS-S71-CUT', 'TESS-S72-SPOC120', 'TESS-S72-QLP', 'TESS-S72-CUT',
                   'TESS-S1751-TESSSPOC', 'TESS-S1751-CUT']
        Bw = 0.3
    else:
        sources = ['K2-c06-PDCSAP-LC', 'K2-c06-K2SFF', 'K2-c06-EVEREST', 'K2-c17-PDCSAP-LC', 'K2-c17-K2SFF', 'K2-c17-EVEREST',
                   'TESS-S46-QLP', 'TESS-S46-QLPSAP', 'TESS-S46-CUT', 'TESS-S91-QLP', 'TESS-S91-QLPSAP', 'TESS-S91-CUT']
        Bw = 1.0
    classes = [cls_from(ev, 'A')]
    D = ev.t0.iloc[-1] - ev.t0.iloc[0]
    cands = []
    for N in range(1, int(D / PMIN) + 1):
        P = D / N
        cyc = np.round((ev.t0.values - ev.t0.values[0]) / P)
        dev = np.abs((ev.t0.values - ev.t0.values[0]) / P - cyc) * P
        if np.any(dev > 5 * np.hypot(ev.err.values, ev.err.values[0]) + TOL):
            continue
        p, cov, res, chi2 = linfit(ev.t0.values, ev.err.values, cyc)
        cands.append(dict(P=p[1], sigP=np.sqrt(cov[1, 1]), tag=f'N={N}', T0_A=p[0], sT0_A=np.sqrt(cov[0, 0]), chi2=chi2,
                          res=res))
    print(f'{len(cands)} candidate periods consistent with all event times')
    known = ev.t0.values

summ = analyse(EPIC, classes, sources, cands, Bw, known, out_rows)
sm = pd.DataFrame(summ)
for c, s in zip(cands, summ):
    s['chi2_timing'] = c['chi2']
    for k in c:
        if k.startswith('T0_') or k.startswith('sT0_'):
            s[k] = c[k]
sm = pd.DataFrame(summ)
sm.to_csv(os.path.join(lib2.HERE, f'period_candidates_{EPIC}.csv'), index=False)
pd.concat(out_rows).to_csv(os.path.join(lib2.HERE, f'epoch_tests_{EPIC}.csv'), index=False)
print('\nALLOWED:')
print(sm[sm.status == 'ALLOWED'][['hypothesis', 'P', 'sigP', 'chi2_timing']].to_string())
