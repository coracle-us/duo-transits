import lightkurve as lk, sys
targets = {'246194159':'301328232','212012030':'14394716','212694013':'178941241'}
for epic,tic in targets.items():
    print('==== EPIC',epic)
    r = lk.search_lightcurve(f'EPIC {epic}')
    print(r)
    print('==== TIC',tic)
    r = lk.search_lightcurve(f'TIC {tic}')
    print(r)
    r = lk.search_tesscut(f'TIC {tic}')
    print(r)
