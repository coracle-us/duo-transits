"""Inspect TESScut cubes: median image, aperture, background and flux time series."""
import sys, os, glob, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib2

epic = sys.argv[1]
fns = lib2._f(epic, '*astrocut.fits')
fig, axs = plt.subplots(len(fns), 3, figsize=(15, 3.2 * len(fns)), squeeze=False)
for row, fn in zip(axs, fns):
    r = lib2.cut_phot(fn, return_info=True)
    h = fits.open(fn)
    hd = h[1].header
    print(os.path.basename(fn), 'napix', r['napix'], 'n_ok', r['ok'].sum(), 'of', len(r['t']),
          'qual!=0', np.sum(r['q'] != 0), '1CRV4P', hd.get('1CRV4P'), '2CRV4P', hd.get('2CRV4P'))
    im = r['img']
    row[0].imshow(np.log10(np.clip(im - np.nanmin(im) + 1, 1, None)), origin='lower')
    yy, xx = np.where(r['ap'])
    row[0].plot(xx, yy, 'r.', ms=3)
    row[0].set_title(os.path.basename(fn)[:30], fontsize=8)
    row[1].plot(r['t'], r['bkg'], 'k.', ms=1)
    row[1].set_title('background/pix')
    ok = r['ok']
    row[2].plot(r['t'][ok], r['f'][ok], 'k.', ms=1)
    row[2].plot(r['t'][~ok], r['f'][~ok], 'r.', ms=1)
    if ok.sum():
        lo, hi = np.percentile(r['f'][ok], [1, 99])
        row[2].set_ylim(lo - 0.3 * (hi - lo), hi + 0.3 * (hi - lo))
    row[2].set_title('aperture flux')
fig.tight_layout()
fig.savefig(os.path.join(lib2.HERE, 'plots', f'cut_inspect_{epic}.png'), dpi=70)
