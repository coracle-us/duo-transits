"""Download target pixel files: K2 C5 (long cadence) and TESS SPOC 2-min S46, S72."""
import os, warnings
import lightkurve as lk
warnings.filterwarnings("ignore")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
sr = lk.search_targetpixelfile("EPIC 211924561", mission="K2", campaign=5)
print(sr)
for r in sr:
    if r.table["exptime"][0] > 1000:
        r.download(download_dir=OUT, quality_bitmask="none")
sr = lk.search_targetpixelfile("TIC 238606535", mission="TESS", author="SPOC")
print(sr)
for r in sr:
    if int(r.table["sequence_number"][0]) in (44, 45, 46, 72) and r.table["exptime"][0] == 120:
        r.download(download_dir=OUT, quality_bitmask="none")
