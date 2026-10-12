"""Download all K2 and TESS light-curve products + TESScut cutouts for the 3 targets."""
import lightkurve as lk, os, warnings
warnings.filterwarnings('ignore')
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
targets = {'246194159':'301328232','212012030':'14394716','212694013':'178941241'}
for epic,tic in targets.items():
    out = os.path.join(D, f'EPIC{epic}'); os.makedirs(out, exist_ok=True)
    for q in [f'EPIC {epic}', f'TIC {tic}']:
        r = lk.search_lightcurve(q)
        for i in range(len(r)):
            try:
                r[i].download(download_dir=out)
                print('ok', q, r.table['mission'][i], r.table['author'][i], r.table['exptime'][i])
            except Exception as e:
                print('FAIL', q, i, e)
    r = lk.search_tesscut(f'TIC {tic}')
    for i in range(len(r)):
        try:
            r[i].download(cutout_size=21, download_dir=out)
            print('ok cut', r.table['mission'][i])
        except Exception as e:
            print('FAIL cut', i, e)
