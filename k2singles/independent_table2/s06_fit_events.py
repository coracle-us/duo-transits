"""Fit every detected event in every independent reduction (generic). Usage: s06_fit_events.py EPIC"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib2, fitlib

K2_12 = ['K2-c12-PDCSAP-LC', 'K2-c12-K2SFF', 'K2-c12-EVERESTallq']
CFG = {
    '246194159': dict(win=1.2, npoly=2, events=[
        ('C12a', 'narrow', 747.38, 0.36, 0.38, K2_12),
        ('C12b', 'wide', 810.88, 0.63, 0.31, K2_12),
        ('S70', 'wide', 3227.99, 0.63, 0.31, ['TESS-S70-TESSSPOC', 'TESS-S70-QLP', 'TESS-S70-CUT']),
        ('S92a', 'wide', 3813.95, 0.63, 0.31, ['TESS-S92-QLP', 'TESS-S92-CUT']),
        ('S92b', 'narrow', 3823.68, 0.36, 0.38, ['TESS-S92-QLP', 'TESS-S92-CUT'])]),
    '212012030': dict(win=0.5, npoly=1, events=[
        ('C05', 'A', 174.47, 0.17, 0.2, ['K2-c05-PDCSAP-LC', 'K2-c05-K2SFF', 'K2-c05-EVERESTallq']),
        ('C18', 'A', 1276.07, 0.17, 0.2, ['K2-c18-PDCSAP-LC', 'K2-c18-K2SFF', 'K2-c18-EVERESTallq']),
        ('S71', 'A', 3258.977, 0.17, 0.2, ['TESS-S71-SPOC120', 'TESS-S71-TESSSPOC', 'TESS-S71-QLP'])]),
    '212694013': dict(win=3.0, npoly=2, events=[
        ('C06', 'A', 242.15, 1.6, 0.006, ['K2-c06-PDCSAP-LC', 'K2-c06-K2SFF', 'K2-c06-EVEREST']),
        ('C17', 'A', 1215.95, 1.6, 0.006, ['K2-c17-PDCSAP-LC', 'K2-c17-K2SFF', 'K2-c17-EVEREST', 'K2-c17-PDCSAP-SC'])]),
}

EPIC = sys.argv[1]
cfg = CFG[EPIC]
segs = {s['src']: s for s in lib2.load_all(EPIC)}
rows = []
ev = cfg['events']
fig, axs = plt.subplots(len(ev), 1, figsize=(9, 3 * len(ev)))
for ax, (lab, typ, tg, T14, dep, srcs) in zip(np.atleast_1d(axs), ev):
    for src in srcs:
        s = segs[src]
        t, f, e = s['t'], s['f'], s['e']
        if src.endswith('-SC'):  # bin K2 short cadence to 10 min for speed
            b = np.floor((t - t.min()) / (10 / 1440)).astype(int)
            nb = np.bincount(b)
            ok = nb > 3
            t = (np.bincount(b, t)[ok] / nb[ok]); f = (np.bincount(b, f)[ok] / nb[ok])
            e = np.full_like(f, np.std(np.diff(f)) / np.sqrt(2))
        r = fitlib.fit_event(t, f, e, tg, T14, dep, s['exptime'] if not src.endswith('-SC') else 10 / 1440,
                             win=cfg['win'], npoly=cfg['npoly'])
        if r is None:
            print('no data', lab, src)
            continue
        rows.append(dict(epic=EPIC, event=lab, type=typ, source=src, t0=r['t0'], t0_err=r['t0_err'],
                         t0_err_cov=r['t0_err_cov'], t0_err_pb=r['t0_err_pb'], depth=r['depth'], depth_err=r['depth_err'],
                         T14_h=r['T14'] * 24, T14_err_h=r['T14_err'] * 24, T12_h=r['T12'] * 24, chi2r=r['chi2r'], n=r['n']))
        print(f"{lab:5s} {typ:6s} {src:22s} t0={r['t0']:.5f} +- {r['t0_err']:.5f} (cov {r['t0_err_cov']:.5f}, pb {r['t0_err_pb']:.5f}) "
              f"depth={r['depth']:.5f}+-{r['depth_err']:.5f} T14={r['T14']*24:.2f}+-{r['T14_err']*24:.2f}h T12={r['T12']*24:.2f}h chi2r={r['chi2r']:.2f} n={r['n']}")
        ax.plot(r['tt'] - tg, r['ff'] / r['base'], '.', ms=2, label=src)
        ax.plot(r['tt'] - tg, r['mm'] / r['base'], '-', lw=1)
    ax.set_title(f'{EPIC} {lab} ({typ})  t - {tg}', fontsize=9)
    ax.legend(fontsize=6)
fig.tight_layout()
fig.savefig(os.path.join(lib2.HERE, 'plots', f'fits_{EPIC}.png'), dpi=80)
pd.DataFrame(rows).to_csv(os.path.join(lib2.HERE, f'event_fits_{EPIC}.csv'), index=False)
