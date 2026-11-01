# guarded: the Sky Patrol client uses multiprocessing, which re-imports this file on macOS
if __name__ == "__main__":
    import pandas as pd
    from astropy.time import Time
    from pyasassn.client import SkyPatrolClient
    lcs = SkyPatrolClient().cone_search(301.564386, 44.456886, radius=1, units="arcsec", download=True)
    d = lcs.data
    print(d.columns.tolist(), len(d))
    print(d.groupby(["phot_filter", "camera", "quality"]).jd.agg(["count", "min", "max"]))
    print(d.sort_values("jd").tail(5)[["jd", "mag", "mag_err", "quality", "phot_filter", "camera"]])
    a, b = Time("2026-10-18 00:00:00", scale="utc").jd, Time("2026-10-25 00:00:00", scale="utc").jd
    print("Oct 18-24 UTC is JD", a, b, "raw points in it:", ((d.jd >= a) & (d.jd < b)).sum())
    print("last raw point:", Time(d.jd.max(), format="jd").iso)
