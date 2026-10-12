"""Period determination + alias exclusion for each star.

For each candidate period we predict every eclipse (each class at its own phase) inside every
data set, and test it with a local measurement:
  * window: points with |t - tp| < T14/2 (in-eclipse) and baseline 0.6*T14 < |t-tp| < 0.6*T14 + Bw
  * a linear baseline is fitted to the baseline points; measured deficit = mean(1 - f/baseline) in eclipse
  * expected deficit = mean(1 - trapezoid model) evaluated at the *actual* in-eclipse timestamps
    (so partial coverage of ingress/egress is handled automatically)
  * sigma = empirical scatter of the same statistic at 300 random epochs in the same data set
    (captures red noise)
Verdict per predicted epoch and data set:
  EVENT  : measured > 0.5*expected and measured/sigma > 5
  FLAT   : (expected - measured) > 5*sigma and measured < 0.5*expected and >=3 in-eclipse points
           ("covered and flat" -> excludes the period)
  NOCOV  : fewer than 3 in-eclipse points
  WEAK   : otherwise (not constraining)
The predicted time is scanned over +-3 sigma_pred and the most event-like value is used, so a
period is only excluded if the epoch is flat for every allowed timing offset.
"""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib2

rng = np.random.default_rng(42)


def local_stat(t, f, tp, T14, T12, depth, exptime, Bw):
    dt = t - tp
    inn = np.abs(dt) < T14 / 2
    bas = (np.abs(dt) > 0.6 * T14) & (np.abs(dt) < 0.6 * T14 + Bw)
    if inn.sum() < 3 or bas.sum() < 5:
        return None
    # need baseline on both sides unless the in-eclipse points are themselves one-sided
    c = np.polyfit(dt[bas], f[bas], 1 if (np.any(dt[bas] < 0) and np.any(dt[bas] > 0)) else 0)
    base = np.polyval(c, dt[inn])
    meas = np.mean(1 - f[inn] / base)
    exp = np.mean(1 - lib2.trapz_binned(t[inn], tp, depth, T14, T12, exptime))
    return meas, exp, inn.sum()


def empirical_sigma(t, f, T14, T12, depth, exptime, Bw, avoid, n=300):
    vals = []
    tries = 0
    while len(vals) < n and tries < 5 * n:
        tries += 1
        tp = rng.uniform(t.min(), t.max())
        if np.any(np.abs(tp - avoid) < T14 + Bw):
            continue
        r = local_stat(t, f, tp, T14, T12, depth, exptime, Bw)
        if r is None or r[2] < 0.5 * T14 / max(np.median(np.diff(t)), 1e-4):
            continue
        vals.append(r[0])
    vals = np.array(vals)
    if len(vals) < 20:
        return np.nan
    return 1.4826 * np.median(np.abs(vals - np.median(vals)))


def test_epoch(seg, tp, sig_pred, cls, Bw):
    T14, T12, depth = cls['T14'], cls['T12'], cls['depth']
    best = None
    for off in np.linspace(-3 * sig_pred, 3 * sig_pred, 7 if sig_pred > 0 else 1):
        r = local_stat(seg['t'], seg['f'], tp + off, T14, T12, depth, seg['exptime'], Bw)
        if r is None:
            continue
        if best is None or r[0] / max(r[1], 1e-9) > best[0] / max(best[1], 1e-9):
            best = r
    if best is None:
        return 'NOCOV', np.nan, np.nan, 0
    meas, exp, nin = best
    s = seg['sig_' + cls['name']]
    if meas > 0.5 * exp and meas / s > 5:
        v = 'EVENT'
    elif (exp - meas) > 5 * s and meas < 0.5 * exp:
        v = 'FLAT'
    else:
        v = 'WEAK'
    return v, meas, exp, nin


def analyse(epic, classes, sources, candidates, Bw, known_times, outrows, label=''):
    segs = [s for s in lib2.load_all(epic) if s['src'] in sources]
    for s in segs:
        for c in classes:
            s['sig_' + c['name']] = empirical_sigma(s['t'], s['f'], c['T14'], c['T12'], c['depth'], s['exptime'], Bw,
                                                    np.array(known_times))
        print(f"   {s['src']:24s} " + ' '.join(f"sigma[{c['name']}]={s['sig_'+c['name']]*1e6:.0f}ppm" for c in classes))
    summary = []
    for cand in candidates:
        P, sigP, tag = cand['P'], cand['sigP'], cand['tag']
        verdicts = []
        for c in classes:
            T0, sT0 = cand['T0_' + c['name']], cand['sT0_' + c['name']]
            for s in segs:
                k0 = int(np.ceil((s['t'].min() - c['T14'] - T0) / P))
                k1 = int(np.floor((s['t'].max() + c['T14'] - T0) / P))
                for k in range(k0, k1 + 1):
                    tp = T0 + k * P
                    sp = np.hypot(sT0, abs(k) * sigP)
                    v, meas, exp, nin = test_epoch(s, tp, sp, c, Bw)
                    known = np.any(np.abs(np.array(known_times) - tp) < c['T14'])
                    verdicts.append(dict(epic=epic, hypothesis=tag, P=P, cls=c['name'], src=s['src'], cycle=k, t_pred=tp,
                                         sig_pred=sp, verdict=v, meas=meas, exp=exp, n_in=nin, known_event=known))
        vd = pd.DataFrame(verdicts)
        outrows.append(vd)
        flats = vd[(vd.verdict == 'FLAT')]
        events_new = vd[(vd.verdict == 'EVENT') & (~vd.known_event)]
        # an epoch is "excluding" if flat in some source and EVENT in none
        excl = []
        for (cl, k), g in vd.groupby(['cls', 'cycle']):
            if (g.verdict == 'FLAT').any() and not (g.verdict == 'EVENT').any():
                excl.append(f"{cl}#{k}@{g.t_pred.iloc[0]:.2f}[{','.join(g[g.verdict=='FLAT'].src.str.replace('TESS-','').str.replace('K2-',''))}]")
        status = 'EXCLUDED' if excl else 'ALLOWED'
        summary.append(dict(epic=epic, hypothesis=tag, P=P, sigP=sigP, status=status, n_excluding_epochs=len(excl),
                            excluding=';'.join(excl[:6]),
                            new_events=';'.join(f"{r.cls}@{r.t_pred:.2f}({r.src})" for r in events_new.itertuples())))
        print(f"  {tag:28s} P={P:12.5f} +- {sigP:.5f}  {status:8s} n_excl={len(excl):3d} {';'.join(excl[:3])}"
              + (f"  NEW EVENTS: {summary[-1]['new_events']}" if len(events_new) else ''))
    return summary
