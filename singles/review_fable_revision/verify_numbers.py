"""Recompute the revised note's numbers from the evidence CSVs. Review only; writes nothing else."""
import json, numpy as np, pandas as pd
from pathlib import Path
R = Path(__file__).resolve().parents[2]
RS = R/'singles/revision_search'; RV = R/'singles/revision'
fc = pd.read_csv(R/'singles/results/final_classification_2026-10-11_rnaas.csv')
print('final_classification rows', len(fc)); print(fc.cls.value_counts().to_string())
dup = [1461133, 3181884, 4689632]
print('duplicates class:', fc[fc.TIC.isin(dup)][['TIC','cls']].to_string(index=False))
# audit
a = pd.read_csv(RV/'audit_34.csv'); print('audit rows', len(a)); print(a.reason.value_counts().to_string())
a31 = a[~a.TIC.isin(dup)]
print('31 by reason', a31.reason.value_counts().to_dict())
d = a31[a31.reason=='d']; print('d notes:'); 
for _,r in d.iterrows(): print('  ',r.TIC, str(r.notes)[:140])
# recovery
rc = pd.read_csv(RS/'recovery_check.csv'); r272 = rc[rc.note_cls!='event-not-recovered']
print('recovery_check rows', len(rc), '272 rows', len(r272), 'recovered by free search', r272.event_recovered.sum())
print(r272[~r272.event_recovered].free_fit_fail_reason.value_counts().to_dict())
# control rates
cf = pd.read_csv(RS/'control_false_alarms.csv'); c272 = cf[cf.note_cls!='event-not-recovered']
ss = c272.n_star_sectors.sum(); print('272 stars', len(c272), 'star-sectors', ss,
  'control passing', c272.control_passing_filters.sum(), round(c272.control_passing_filters.sum()/ss,4),
  'extra passing', c272.extra_passing_filters.sum(), round(c272.extra_passing_filters.sum()/ss,4),
  'control stars', (c272.control_passing_filters>0).sum(), 'extra stars', (c272.extra_passing_filters>0).sum())
# SNR table
ev = pd.read_csv(RS/'extra_events.csv'); ce = pd.read_csv(RS/'data/control_events_all.csv')
e272 = ev[ev.note_cls!='event-not-recovered']; print('extra_events 272 rows', len(e272), 'auto cat', e272.auto_category.value_counts().to_dict())
epass = e272[e272.auto_category!='d']
tics272 = set(c272.TIC); cpass = ce[(ce.TIC.isin(tics272)) & (ce['filter'].isna() | (ce['filter']==''))]
print('control filter values', ce['filter'].value_counts(dropna=False).to_dict())
bins=[8,10,15,20,50,1e9]
print('SNR bins dips vs inverted (passing filters):')
for lo,hi in zip(bins[:-1],bins[1:]):
    print(f'  {lo}-{hi}: dips {((epass.snr>=lo)&(epass.snr<hi)).sum()}  inv {((cpass.snr>=lo)&(cpass.snr<hi)).sum()}')
print('passing extra', len(epass), 'passing control', len(cpass))
print('vetted-real by class:', e272[e272.manual_verdict=='real'].groupby('note_cls').TIC.nunique().to_dict(),
      'events', e272[e272.manual_verdict=='real'].groupby('note_cls').size().to_dict())
print('manual_verdict counts', e272.manual_verdict.value_counts(dropna=False).to_dict())
ps = pd.read_csv(RS/'per_star_summary.csv'); p272 = ps[ps.note_cls!='event-not-recovered']
print('per_star 272', len(p272)); 
ch = p272[~p272.suggested_class.str.contains('unchanged')]
print('suggested_class changes:'); print(ch[['TIC','note_cls','suggested_class']].to_string(index=False))
sing = p272[p272.note_cls=='single']; print('singles', len(sing), 'with >=1 extra', (sing.n_extra_detected>0).sum(), 'with real', ((sing.n_a_real+sing.n_b_real+sing.n_c_real)>0).sum(), 'possible', (sing.n_possible>0).sum())
# 159614322 / 229055790
for t in [159614322,229055790,147035057]:
    print(t, ps[ps.TIC==t][['n_extra_detected','n_extra_two_reductions','n_a_real','n_b_real','n_possible','n_control','n_control_passing']].to_dict('records'))
    print('   events:', ev[ev.TIC==t][['T_btjd','sector','product','depth_ppm','snr','auto_category','manual_verdict','category']].to_string(index=False))
# 147035057 cycles and e cos w
P=37.31497; T0=1515.0661
for T in [1486.6, 3688.4]:  # approx from vetting windows
    print('phase check', T, ((T-T0)/P)%1)
print('cycles S6->S87 approx', (3688.4-1486.6)/P)
ph=0.239; print('e cos w ~', (np.pi/2)*((1-ph)-0.5), 'or', (np.pi/2)*(0.5-ph))
# photometry
pv = pd.read_csv(RV/'photometry_validation.csv')
print(pv[['TIC','Tmag','n_ap_own','n_ap_spoc','ratio_own_sap','ratio_own_pdc','bkg_own_per_px','bkg_spoc_per_px','depth_ratio_own_pdc','midtime_diff_mean_min','midtime_diff_rms_min','seg_level_spread_ppt_1h','seg_level_spread_ppt_05d']].to_string(index=False))
pv['bkg_ratio']=pv.bkg_own_per_px/pv.bkg_spoc_per_px; print('bkg ratio', pv.bkg_ratio.round(3).tolist())
pv['p2p_ratio']=pv.p2p_ppm_own/pv.p2p_ppm_sap; print('p2p own/sap', pv.p2p_ratio.round(2).tolist(), 'median', pv.p2p_ratio.median())
