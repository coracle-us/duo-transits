"""List all TESS light curve products and TESScut sectors for the four targets."""
import warnings; warnings.filterwarnings("ignore")
import lightkurve as lk
TICS = [275647095, 241550976, 107113048, 230133616]
for t in TICS:
    sr = lk.search_lightcurve(f"TIC {t}", mission="TESS")
    print(f"=== TIC {t}: {len(sr)} light curves")
    if len(sr):
        tab = sr.table
        for r in tab:
            print(f"  {r['mission']:<16} {r['author']:<10} exp={r['exptime']}  {r['productFilename'] if 'productFilename' in tab.colnames else ''}")
    sc = lk.search_tesscut(f"TIC {t}")
    print(f"  TESScut sectors: {[m.split()[-1] for m in sc.table['mission']]}")
