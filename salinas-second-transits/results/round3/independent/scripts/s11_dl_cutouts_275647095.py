"""Download TESScut cutouts (21x21) for TIC 275647095 sectors with full transits that had no
cutout yet (S37, S64, S90), for difference imaging."""
import warnings, os; warnings.filterwarnings("ignore")
import lightkurve as lk
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sc = lk.search_tesscut("TIC 275647095")
want = [i for i, m in enumerate(sc.table["mission"]) if int(m.split()[-1]) in (37, 64, 90)]
print([sc.table["mission"][i] for i in want])
sc[want].download_all(cutout_size=21, download_dir=os.path.join(ROOT, "data", "tesscut_diff"))
