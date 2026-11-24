"""Download 21x21 TESScut cutouts for difference imaging.
Usage: s11_dl_cutouts.py TIC SECTOR [SECTOR ...]   (into data/tesscut_diff/)"""
import sys, os, warnings; warnings.filterwarnings("ignore")
import lightkurve as lk
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tic, secs = int(sys.argv[1]), {int(s) for s in sys.argv[2:]}
sc = lk.search_tesscut(f"TIC {tic}")
want = [i for i, m in enumerate(sc.table["mission"]) if int(m.split()[-1]) in secs]
print([sc.table["mission"][i] for i in want])
sc[want].download_all(cutout_size=21, download_dir=os.path.join(ROOT, "data", "tesscut_diff"))
