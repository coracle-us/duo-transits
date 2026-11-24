import pandas as pd, numpy as np
R="/Users/mcbee/Documents/GitHub/coracle/duo-transits-be-star/be-star-pulsations/results/"
sm=pd.read_csv(R+"full/sector_metrics.csv", low_memory=False)
samp=pd.read_csv(R+"full/sample.csv").set_index("tic"); sw=pd.read_csv(R+"full/star_wise.csv").set_index("tic")
sm["preferred"]=sm.tic.map(samp.preferred); sm["w1med"]=sm.tic.map(sw.w1_median)
b=lambda c: sm[c].fillna(False).astype(bool)
d=sm[b("metrics_ok")&b("crowd_ok")&b("wise_usable")&np.isfinite(sm.disk_excess)].copy()
d=d[d.groupby("tic").tic.transform("size")>=2].copy()
print("main",d.tic.nunique(),len(d))
M=["A1","B1","B2","B3","red_noise","flicker_rms","white_noise"]
for m in M: d["l"+m]=np.log(d[m])
rng=np.random.default_rng(1)
def beta(df,ycol,xcol="disk_excess",grp="tic",boot=500):
    df=df[np.isfinite(df[ycol])&np.isfinite(df[xcol])]
    df=df[df.groupby(grp)[ycol].transform("size")>=2]
    xd=df[xcol]-df.groupby(grp)[xcol].transform("mean"); yd=df[ycol]-df.groupby(grp)[ycol].transform("mean")
    g=pd.DataFrame(dict(s=df.tic.values,sxy=(xd*yd).values,sxx=(xd*xd).values)).groupby("s").sum()
    bt=g.sxy.sum()/g.sxx.sum()
    W=rng.multinomial(len(g),np.full(len(g),1/len(g)),size=boot); bs=(W@g.sxy.values)/(W@g.sxx.values)
    return bt,bs.std(),df.tic.nunique(),len(df)
def show(tag,df,ms=M,**k):
    print(tag," ".join(f"{m}:{beta(df,'l'+m,**k)[0]:+.2f}±{beta(df,'l'+m,**k)[1]:.2f}" for m in ms),"N*",beta(df,'l'+ms[0],**k)[2],"Nsec",beta(df,'l'+ms[0],**k)[3])
show("A main      ",d)
# B two-way FE: star + sector, iterate
def twoway(df,ycol):
    x=df.disk_excess.values.copy(); y=df[ycol].values.copy()
    for _ in range(50):
        for g in (df.tic.values,df.sector.values):
            s=pd.Series(x); x=x-s.groupby(g).transform("mean").values
            s=pd.Series(y); y=y-s.groupby(g).transform("mean").values
    return (x*y).sum()/(x*x).sum()
print("B star+sector FE"," ".join(f"{m}:{twoway(d,'l'+m):+.2f}" for m in M))
d["grp_cad"]=d.tic.astype(str)+"_"+d.cadence_s.round().astype(str)
show("H within star x cadence",d,grp="grp_cad")
print("cadence counts",d.cadence_s.round().value_counts().to_dict())
# C leverage
xr=d.groupby("tic").disk_excess.agg(lambda v:v.max()-v.min()); d["xr"]=d.tic.map(xr)
xd=d.disk_excess-d.groupby("tic").disk_excess.transform("mean")
sxx=(xd**2).groupby(d.tic).sum().sort_values(ascending=False)
cs=sxx.cumsum()/sxx.sum()
print("C Sxx share: top1%% stars %.2f top5%% %.2f top10%% %.2f top20%% %.2f; stars for 50%%: %d, 90%%: %d of %d"%(cs.iloc[10],cs.iloc[53],cs.iloc[107],cs.iloc[216],(cs<0.5).sum()+1,(cs<0.9).sum()+1,len(cs)))
print("x-range quantiles",xr.quantile([.1,.25,.5,.75,.9,.99]).round(3).to_dict())
for lo,hi in [(0,0.05),(0.05,0.1),(0.1,0.2),(0.2,0.4),(0.4,9)]:
    s=d[(d.xr>=lo)&(d.xr<hi)]
    # per-star slopes
    out=[]
    for m in ["A1","B1","flicker_rms","white_noise"]:
        sl=s.groupby("tic").apply(lambda g: np.polyfit(g.disk_excess,g["l"+m],1)[0] if np.ptp(g.disk_excess)>0 else np.nan, include_groups=False).dropna()
        bt=beta(s,"l"+m)
        out.append(f"{m}: b={bt[0]:+.2f}±{bt[1]:.2f} pos={np.mean(sl>0):.2f}")
    print(f"  range {lo}-{hi}: N*={s.tic.nunique()} | "+" | ".join(out))
# drop top-leverage stars
for n in (10,50,100):
    show(f"C drop top {n} Sxx stars",d[~d.tic.isin(sxx.index[:n])],ms=["A1","B1","flicker_rms","white_noise"])
# D white noise
d["wn_med"]=d.groupby("tic").white_noise.transform("median")
q=d.drop_duplicates("tic").wn_med.quantile([1/3,2/3]).values
for i,(lo,hi) in enumerate([(0,q[0]),(q[0],q[1]),(q[1],1e9)]):
    s=d[(d.wn_med>=lo)&(d.wn_med<hi)]
    show(f"D white-noise tercile {i} (median wn {s.wn_med.median():.3f} ppt, W1 {s.w1med.median():.1f})",s,ms=["A1","B1","flicker_rms","white_noise"])
for pr in (["QLP"],["SPOC120","TESS-SPOC"]):
    for i,(lo,hi) in enumerate([(0,q[0]),(q[0],q[1]),(q[1],1e9)]):
        s=d[(d.wn_med>=lo)&(d.wn_med<hi)&d.preferred.isin(pr)]
        show(f"D {pr[0]} tercile {i}",s,ms=["A1","B1","white_noise"])
# within-star correlation of white with B1 residuals
def dm(c): return d[c]-d.groupby("tic")[c].transform("mean")
X=dm("disk_excess"); 
for a_,b_ in [("lwhite_noise","lB1"),("lwhite_noise","lA1"),("lwhite_noise","lB3"),("lA1","lB1")]:
    ya=dm(a_); yb=dm(b_); ra=ya-(X*ya).sum()/(X*X).sum()*X; rb=yb-(X*yb).sum()/(X*X).sum()*X
    print(f"D within-star corr (x removed) {a_} vs {b_}: {np.corrcoef(ra,rb)[0,1]:.2f}")
# white beta controlling for lB1 and for lA1: 2-var OLS on demeaned
for ctrl in ("lB1","lA1","lflicker_rms"):
    A=np.c_[X,dm(ctrl)]; co=np.linalg.lstsq(A,dm("lwhite_noise"),rcond=None)[0]; print(f"D white on x + {ctrl}: beta_x={co[0]:+.2f}, coef_{ctrl}={co[1]:+.2f}")
A=np.c_[X,dm("lwhite_noise")]
for m in ["lA1","lB1","lflicker_rms"]:
    co=np.linalg.lstsq(A,dm(m),rcond=None)[0]; print(f"D {m} on x + lwhite: beta_x={co[0]:+.2f}, coef_white={co[1]:+.2f}")
# E f1 stability
d["lowf1"]=(d.f1<0.8).astype(float)
print("E frac sectors f1<0.8:",d.lowf1.mean().round(3),"; LPM slope of I(f1<0.8) on x within star: %.3f±%.3f"%beta(d,"lowf1")[:2])
fr=d.groupby("tic").f1.agg(lambda v:v.max()/v.min()); d["f1ratio"]=d.tic.map(fr)
print("E stars with f1 stable within 10%:",(fr<1.1).sum(),"of",len(fr), "; within 2%:",(fr<1.02).sum())
show("E stable-f1 stars (<10%)",d[d.f1ratio<1.1],ms=["A1","B1","flicker_rms","white_noise"])
show("E stable-f1 & f1>0.8",d[(d.f1ratio<1.1)&(d.groupby('tic').f1.transform('min')>0.8)],ms=["A1","B1","flicker_rms","white_noise"])
show("E unstable-f1 stars",d[d.f1ratio>=1.1],ms=["A1","B1","flicker_rms","white_noise"])
show("E snr1>=8 all sectors",d[d.groupby('tic').snr1.transform('min')>=8],ms=["A1","B1","flicker_rms","white_noise"])
print("E snr1 quantiles",d.snr1.quantile([.1,.25,.5,.75,.9]).round(1).to_dict(), " frac snr1<4:",(d.snr1<4).mean().round(3))
# F QLP vs SPOC difference
for m in ["A1","B1","flicker_rms","white_noise"]:
    bq=beta(d[d.preferred=="QLP"],"l"+m,boot=2000); bs_=beta(d[d.preferred.isin(["SPOC120","TESS-SPOC"])],"l"+m,boot=2000)
    print(f"F {m}: QLP {bq[0]:+.2f}±{bq[1]:.2f}  SPOC {bs_[0]:+.2f}±{bs_[1]:.2f}  diff {bq[0]-bs_[0]:+.2f}±{np.hypot(bq[1],bs_[1]):.2f}")
# sign test numbers
sl=d.groupby("tic").apply(lambda g: np.polyfit(g.disk_excess,g.lA1,1)[0] if np.ptp(g.disk_excess)>0 else np.nan, include_groups=False).dropna()
from math import comb
n=len(sl); k=int((sl>0).sum()); 
from scipy.stats import binomtest, wilcoxon
print("sign A1: n",n,"pos",k,"frac neg",1-k/n,"binom p",binomtest(k,n).pvalue, "wilcoxon p",wilcoxon(sl).pvalue)
for m in ["B1","flicker_rms","red_noise","white_noise"]:
    s2=d.groupby("tic").apply(lambda g: np.polyfit(g.disk_excess,g["l"+m],1)[0] if np.ptp(g.disk_excess)>0 else np.nan, include_groups=False).dropna()
    print("sign",m,"pos frac",round((s2>0).mean(),3),"binom p",binomtest(int((s2>0).sum()),len(s2)).pvalue)
print("sectors per star", d.groupby("tic").size().describe().round(2).to_dict())
print("product in main (sector rows):", d["product"].value_counts().to_dict())
