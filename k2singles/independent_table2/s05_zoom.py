"""Zoom plots around given times for all (non-redundant) sources. Usage: s05_zoom.py EPIC halfwidth t1 t2 ..."""
import sys, os, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib2

epic = sys.argv[1]
hw = float(sys.argv[2])
times = [float(x) for x in sys.argv[3:]]
tag = os.environ.get('ZTAG', '')
segs = lib2.load_all(epic)
fig, axs = plt.subplots(len(times), 1, figsize=(10, 3.2 * len(times)), squeeze=False)
for ax, tc in zip(axs[:, 0], times):
    k = 0
    for s in segs:
        m = np.abs(s['t'] - tc) < hw
        if m.sum() < 3:
            continue
        ax.plot(s['t'][m] - tc, s['f'][m] - 0.0 * k, '.', ms=2.5, label=f"{s['src']}")
        k += 1
    ax.set_title(f'EPIC {epic}  t0={tc}', fontsize=9)
    ax.legend(fontsize=6, ncol=2)
    ax.set_xlabel('t - t0 [d]')
fig.tight_layout()
fn = os.path.join(lib2.HERE, 'plots', f'zoom_{epic}{tag}.png')
fig.savefig(fn, dpi=80)
print(fn)
