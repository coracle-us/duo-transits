"""Stage 2a (full run): fetch IRSA NEOWISE-R single exposures for every non-EB star up front,
independent of stage 1, into the shared cache used by s2_wise.py. Resumable."""
import sys
import pandas as pd
from common import RESULTS
from s2_wise import fetch_irsa

if __name__ == "__main__":
    chunk = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    cat = pd.read_parquet(RESULTS / "catalogue.parquet")
    stars = cat[~cat.eb_type.isin([1, 2, 3])][["tic", "ra", "dec"]].reset_index(drop=True)
    fetch_irsa(stars, chunk=chunk)
    print("done")
