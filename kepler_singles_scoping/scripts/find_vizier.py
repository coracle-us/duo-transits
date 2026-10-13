from astroquery.vizier import Vizier
import sys
for kw in sys.argv[1:]:
    try:
        c = Vizier.find_catalogs(kw, max_catalogs=40)
        print("##", kw)
        for k, v in c.items():
            print("  ", k, "|", v.description[:110])
    except Exception as e:
        print('fail', kw, e)
