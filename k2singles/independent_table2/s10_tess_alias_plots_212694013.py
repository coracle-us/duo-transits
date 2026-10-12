"""EPIC 212694013: visual check of every TESS epoch predicted by the K2-allowed aliases (N<=17)."""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib2

EPIC = '212694013'
d = pd.read_csv(os.path.join(lib2.HERE, f'epoch_tests_{EPIC}.csv'))
ev = pd.read_csv(os.path.join(lib2.HERE, f'events_combined_{EPIC}.csv'))
T14 = ev.T14_h.median() / 24; T12 = ev.T12_h.median() / 24; dep = ev.depth.median()
d['N'] = d.hypothesis.str[2:].astype(int)
sel = d[(d.N <= 17) & d.src.str.startswith('TESS') & (d.verdict != 'NOCOV')]
ep = sel.groupby(['N', 'cycle']).agg(t_pred=('t_pred', 'first'), verdicts=('verdict', lambda v: ','.join(v)),
                                      srcs=('src', lambda v: ','.join(s.replace('TESS-', '') for s in v))).reset_index()
print(ep.to_string())
segs = {s['src']: s for s in lib2.load_all(EPIC) if s['src'].startswith('TESS')}


def binlc(t, f, w):
    b = np.floor((t - t.min()) / w).astype(int)
    n = np.bincount(b)
    ok = n >= 3
    return np.bincount(b, t)[ok] / n[ok], np.bincount(b, f)[ok] / n[ok]


fig, axs = plt.subplots(len(ep), 1, figsize=(10, 2.6 * len(ep)), squeeze=False)
for ax, r in zip(axs[:, 0], ep.itertuples()):
    for k, src in enumerate(['S46-QLP', 'S46-QLPSAP', 'S46-CUT', 'S91-QLP', 'S91-QLPSAP', 'S91-CUT']):
        s = segs.get('TESS-' + src)
        if s is None:
            continue
        m = np.abs(s['t'] - r.t_pred) < 3.5
        if m.sum() < 10:
            continue
        tb, fb = binlc(s['t'][m], s['f'][m], 1 / 24)
        med = np.median(fb)
        ax.plot(tb - r.t_pred, fb / med - 0.012 * (k % 3), '.', ms=3, label=src)
        mod = lib2.trapz(tb, r.t_pred, dep, T14, T12)
        ax.plot(tb - r.t_pred, mod - 0.012 * (k % 3), 'k-', lw=0.7)
    ax.set_title(f'N={r.N} P={973.789/r.N:.3f}  cycle {r.cycle}  t_pred={r.t_pred:.2f}  verdicts {r.verdicts} ({r.srcs})', fontsize=8)
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(0.97, 1.01)
    ax.legend(fontsize=6, ncol=3, loc='lower left')
fig.tight_layout()
fig.savefig(os.path.join(lib2.HERE, 'plots', f'tess_alias_epochs_{EPIC}.png'), dpi=70)
