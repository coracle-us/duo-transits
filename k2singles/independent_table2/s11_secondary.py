"""Search for secondary eclipses at any phase for each allowed period.
Each data set is detrended (robust running median, window Wd) with the known primary events masked,
then per phase bin we compute a box depth (durations 0.5-2 x primary T14) from all data sets combined
(inverse-variance weighted, per-set sigma = empirical point scatter x red-noise factor).
Reports the most significant dip, and the fraction of phase covered.
Usage: s11_secondary.py EPIC P1 [P2 ...]"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib2
from s04_find_events import robust_detrend

EPIC = sys.argv[1]
Ps = [float(x) for x in sys.argv[2:]]
SRC = {
    '212012030': (['K2-c05-K2SFF', 'K2-c18-K2SFF', 'TESS-S71-SPOC120', 'TESS-S72-SPOC120', 'TESS-S1751-TESSSPOC'], 1.0),
    '212694013': (['K2-c06-K2SFF', 'K2-c17-K2SFF', 'TESS-S46-QLP', 'TESS-S46-CUT', 'TESS-S91-QLP'], 5.0),
    '246194159': (['K2-c12-K2SFF', 'K2-c19-K2SFF', 'TESS-S42-SPOC120', 'TESS-S70-TESSSPOC', 'TESS-S92-QLP'], 2.5),
}
srcs, Wd = SRC[EPIC]
ev = pd.read_csv(os.path.join(lib2.HERE, f'events_combined_{EPIC}.csv'))
T14 = ev.T14_h.min() / 24
segs = [s for s in lib2.load_all(EPIC) if s['src'] in srcs]
data = []
for s in segs:
    t, f = s['t'], s['f']
    inev = np.zeros_like(t, bool)
    for r in ev.itertuples():
        inev |= np.abs(t - r.t0) < 0.75 * r.T14_h / 24
    t, f = t[~inev], f[~inev]
    fd, sig = robust_detrend(t, f, Wd)
    data.append((s['src'], t, fd, sig))
T0 = ev.t0.iloc[0]
rows = []
fig, axs = plt.subplots(len(Ps), 1, figsize=(10, 2.6 * len(Ps)), squeeze=False)
for ax, P in zip(axs[:, 0], Ps):
    best = (0, None)
    phgrid = np.arange(0, 1, min(T14 / P / 4, 0.002))
    cover = np.zeros_like(phgrid, bool)
    for durf in (0.5, 1.0, 1.5, 2.0):
        D = durf * T14
        num = np.zeros_like(phgrid); den = np.zeros_like(phgrid)
        for src, t, fd, sig in data:
            ph = ((t - T0) / P) % 1
            # red-noise factor for this duration
            nb = max(int((t.max() - t.min()) / D), 1)
            idx = np.floor((t - t.min()) / D).astype(int)
            cnt = np.bincount(idx)
            means = np.bincount(idx, fd - 1)[cnt > 3] / cnt[cnt > 3]
            npts = np.median(cnt[cnt > 3]) if (cnt > 3).any() else 1
            beta = max(1, 1.4826 * np.median(np.abs(means - np.median(means))) / (sig / np.sqrt(npts))) if len(means) > 5 else 1
            for i, p0 in enumerate(phgrid):
                dph = (ph - p0 + 0.5) % 1 - 0.5
                m = np.abs(dph * P) < D / 2
                if m.sum() >= 3:
                    dep = 1 - np.mean(fd[m])
                    var = (beta * sig) ** 2 / m.sum()
                    num[i] += dep / var; den[i] += 1 / var
                    if durf == 1.0:
                        cover[i] = True
        with np.errstate(invalid='ignore', divide='ignore'):
            dep = num / den
            snr = dep * np.sqrt(den)
        k = np.nanargmax(np.where(np.isfinite(snr), snr, -1))
        if snr[k] > best[0]:
            best = (snr[k], dict(phase=phgrid[k], depth=dep[k], dur_h=D * 24, snr=snr[k]))
        if durf == 1.0:
            ax.plot(phgrid, dep * 1e6, 'k-', lw=0.6)
            ax.fill_between(phgrid, -3e6 / np.sqrt(den), 3e6 / np.sqrt(den), color='0.85')
    b = best[1]
    print(f'P={P:.5f}  phase coverage {cover.mean()*100:.0f}%  best dip: phase={b["phase"]:.4f} depth={b["depth"]*1e6:.0f}ppm '
          f'dur={b["dur_h"]:.1f}h SNR={b["snr"]:.1f}')
    rows.append(dict(epic=EPIC, P=P, phase_coverage=cover.mean(), best_phase=b['phase'], best_depth=b['depth'],
                     best_dur_h=b['dur_h'], best_snr=b['snr']))
    ax.set_title(f'{EPIC}  P={P:.4f}: box depth (D=T14) vs phase (grey = +-3 sigma); coverage {cover.mean()*100:.0f}%', fontsize=8)
    ax.set_ylabel('ppm')
fig.tight_layout()
fig.savefig(os.path.join(lib2.HERE, 'plots', f'secondary_search_{EPIC}.png'), dpi=70)
pd.DataFrame(rows).to_csv(os.path.join(lib2.HERE, f'secondary_search_{EPIC}.csv'), index=False)
