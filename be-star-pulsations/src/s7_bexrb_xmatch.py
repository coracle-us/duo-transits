"""Stage 7 (literature follow-up): are any Be X-ray binaries in the sample?

Niwano et al. 2024 (MNRAS, doi:10.1093/mnras/stae2160) find the main TESS pulsation
amplitudes *falling* as the J-band disk excess grows in 4 of 5 giant-outbursting Be X-ray
binaries. We find A1 rising with W1 excess. This script checks for overlap.

SIMBAD and VizieR (CDS) were unavailable on 2026-11-23, so there is no complete
X-ray-binary list. Instead:
  - Niwano et al.'s 5 sources;
  - an ad hoc list of well-known Galactic Be X-ray binaries (names from memory). This is
    NOT complete; it only checks whether famous systems are in the sample.
Each name is resolved with the MAST name lookup (read-only; it serves cached SIMBAD/NED
positions) and cone-matched (5 arcsec) to the Jia et al. 2026 catalogue.
Output: results/<run>/bexrb_xmatch.csv. Lookups are cached in CACHE/simbad/mast_lookup.json.
"""
import json, time, urllib.parse, urllib.request
import numpy as np, pandas as pd
from common import PILOT, CAT_DIR, CACHE

NIWANO = ["4U 0115+634", "GRO J2058+42", "RX J0440.9+4431", "SAX J2103.5+4545", "V0332+53"]
OTHERS = ["X Per", "A0535+26", "GX 304-1", "EXO 2030+375", "4U 1145-619", "A1118-616",
          "GRO J1008-57", "KS 1947+300", "Swift J0243.6+6124", "GS 0834-430", "2S 1417-624",
          "XTE J1946+274", "RX J0146.9+6121", "MXB 0656-072", "SAX J2239.3+6116",
          "Cep X-4", "SAX J0635.2+0533", "IGR J01363+6610", "RX J0812.4-3114", "3A 0726-260",
          "RX J2030.5+4751", "4U 2206+54", "IGR J06074+2205", "XTE J0658-073",
          "4U 1901+03", "GRO J1750-27", "XTE J1543-568", "1A 1744-361",
          "IGR J21343+4738", "SAX J2103.5+4545", "LS I +61 303", "PSR B1259-63",
          "HESS J0632+057", "MWC 148"]
RADIUS = 5.0  # arcsec
LOOKUP = CACHE / "simbad" / "mast_lookup.json"


def mast_lookup(name):
    req = json.dumps({"service": "Mast.Name.Lookup", "params": {"input": name, "format": "json"}})
    data = urllib.parse.urlencode({"request": req}).encode()
    r = urllib.request.Request("https://mast.stsci.edu/api/v0/invoke", data=data,
                               headers={"User-Agent": "research-script/1.0"})
    with urllib.request.urlopen(r, timeout=60) as f:
        d = json.load(f)
    res = d.get("resolvedCoordinate") or []
    if not res:
        return None
    x = res[0]
    return dict(canonical=x["canonicalName"], ra=x["ra"], dec=x["decl"], resolver=x["resolver"])


def main():
    LOOKUP.parent.mkdir(parents=True, exist_ok=True)
    cache = json.loads(LOOKUP.read_text()) if LOOKUP.exists() else {}
    names = NIWANO + [n for n in OTHERS if n not in NIWANO]
    for n in names:
        if n not in cache:
            cache[n] = mast_lookup(n)
            time.sleep(0.5)
    LOOKUP.write_text(json.dumps(cache, indent=1))
    cat = pd.read_csv(CAT_DIR / "bet_table.csv")
    samp = pd.read_csv(PILOT / "sample.csv")
    status = samp.set_index("tic").status
    cra, cdec = np.radians(cat.ra.values), np.radians(cat.dec.values)
    rows = []
    for n in names:
        c = cache[n]
        row = dict(name=n, niwano=n in NIWANO, resolved=c is not None)
        if c:
            ra, dec = np.radians(c["ra"]), np.radians(c["dec"])
            cosd = (np.sin(dec) * np.sin(cdec) + np.cos(dec) * np.cos(cdec) * np.cos(cra - ra))
            sep = np.degrees(np.arccos(np.clip(cosd, -1, 1))) * 3600
            k = int(np.argmin(sep))
            row.update(canonical=c["canonical"], ra=c["ra"], dec=c["dec"],
                       nearest_sep_arcsec=float(sep[k]))
            if sep[k] <= RADIUS:
                tic = int(cat.TIC.iloc[k])
                row.update(in_catalogue=True, tic=tic, catalogue_name=cat.name.iloc[k],
                           sample_status=status.get(tic, "not in sample.csv"))
            else:
                row.update(in_catalogue=False)
        rows.append(row)
        print(row, flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(PILOT / "bexrb_xmatch.csv", index=False)
    print("resolved", int(out.resolved.sum()), "of", len(out), "; in catalogue",
          int(out.get("in_catalogue", pd.Series(dtype=bool)).fillna(False).sum()))


if __name__ == "__main__":
    main()
