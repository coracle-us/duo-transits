"""
Leads triage, step 3a: pixel data for the two leads that a difference image can still decide (Jon's OK, 2026-11-24):
  EPIC 245929407 = TIC 49733422   K2 C12 TPF (event BTJD 760.12), TESScut S69 (event 3196.34)
  EPIC 211894612 = TIC 14678155   K2 C5 and C18 TPFs (events 182.86, 1297.24), TESScut S44 (event 2523.03)
Both TESS events are in FFI-only sectors, so the TESS pixels are 15x15 TESScut cutouts, as in scripts/diff_image.py.

Output: k2singles/leads_triage/data/pixels/{k2,tesscut}/<id>/ (data/ is a symlink to
duo-transits-cache/leads_triage, outside git). Re-running skips files already present.
Usage:  python k2singles/leads_triage/pixels_download.py
"""
from pathlib import Path

import lightkurve as lk
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import Tesscut

HERE = Path(__file__).resolve().parent
PIX = HERE / "data" / "pixels"
K2 = {"245929407": [12], "211894612": [5, 18]}
TESS = {49733422: [69], 14678155: [44]}


def main():
    import pandas as pd
    cat = pd.read_csv(HERE.parent / "data" / "k2_singles.csv", dtype={"EPIC": str}).drop_duplicates("EPIC").set_index("EPIC")
    for epic, camps in K2.items():
        d = PIX / "k2" / epic
        d.mkdir(parents=True, exist_ok=True)
        for c in camps:
            if list(d.rglob(f"ktwo{epic}-c{c:02d}_lpd-targ.fits*")):
                print("have", epic, c)
                continue
            sr = lk.search_targetpixelfile(f"EPIC {epic}", mission="K2", campaign=c, exptime="long")
            print(sr)
            sr[0].download(download_dir=str(d), quality_bitmask="none")
    for tic, secs in TESS.items():
        epic = cat.index[cat.TIC == tic][0]
        r = cat.loc[epic]
        d = PIX / "tesscut" / str(tic)
        d.mkdir(parents=True, exist_ok=True)
        for s in secs:
            if list(d.glob(f"*s{s:04d}*.fits")):
                print("have", tic, s)
                continue
            m = Tesscut.download_cutouts(coordinates=SkyCoord(r.ra * u.deg, r.dec * u.deg), size=15, sector=s, path=str(d))
            print(m)
    for f in sorted(PIX.rglob("*.fits*")):
        print(f.relative_to(PIX), f"{f.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
