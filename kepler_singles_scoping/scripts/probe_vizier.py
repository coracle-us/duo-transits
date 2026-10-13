from astroquery.vizier import Vizier
import sys
for cat in sys.argv[1:]:
    try:
        r = Vizier.get_catalog_metadata(catalog=cat)
        print("==", cat, r['title'][0] if 'title' in r.colnames else r)
    except Exception as e:
        print("meta fail", cat, e)
    try:
        v = Vizier(row_limit=3, columns=['**'])
        t = v.get_catalogs(cat)
        for tb in t:
            print('  table', tb.meta.get('name'), len(tb), tb.colnames)
    except Exception as e:
        print('  fail', e)
