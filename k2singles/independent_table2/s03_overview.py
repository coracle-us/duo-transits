"""Overview: coverage table + full light-curve plots for every source, each star."""
import sys, os, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib2

STARS = {'246194159': 747.38780, '212012030': 174.46685, '212694013': 242.175}
OUT = os.path.join(lib2.HERE, 'plots')
os.makedirs(OUT, exist_ok=True)

for epic, ep in STARS.items():
    segs = lib2.load_all(epic)
    print('=====', epic, 'epoch', ep)
    for s in segs:
        t = s['t']
        rms = 1.4826 * np.median(np.abs(np.diff(s['f']))) / np.sqrt(2)
        print(f"  {s['src']:28s} n={len(t):6d} {t.min():9.3f}-{t.max():9.3f}  p2p_rms={rms*1e6:8.0f} ppm  min={s['f'].min():.3f}")
    show = [s for s in segs if not any(k in s['src'] for k in ['SAP-', 'SAPallq', 'QLPSAP', 'SPOC120SAP'])]
    fig, axs = plt.subplots(len(show), 1, figsize=(14, 2.0 * len(show)))
    for ax, s in zip(axs, show):
        ax.plot(s['t'], s['f'], 'k.', ms=1)
        ax.set_title(s['src'], fontsize=8, loc='left')
        lo, hi = np.percentile(s['f'], [0.5, 99.5])
        ax.set_ylim(min(lo, s['f'].min()) - 0.01, hi + 0.01)
        if s['t'].min() < ep < s['t'].max():
            ax.axvline(ep, color='r', alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f'overview_{epic}.png'), dpi=80)
    plt.close(fig)
