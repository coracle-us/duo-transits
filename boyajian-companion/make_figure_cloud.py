import sys; sys.path.insert(0,'.')
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
t,f,bl,quiet=np.load('quiet.npy'); T0=3897.164; T14,T23,DEP=0.875,0.50,0.011
deficit=bl-f; vis=deficit<0.008
def status(tc):
    a,b=np.searchsorted(t,[tc-T23/2,tc+T23/2])
    if vis[a:b].sum()>=3: return 2
    return 1 if b-a>=3 else 0
def model(dt):
    a=np.abs(dt); m=np.zeros_like(dt); m[a<=T23/2]=DEP
    s=(a>T23/2)&(a<T14/2); m[s]=DEP*(T14/2-a[s])/((T14-T23)/2); return m
Ps=np.arange(770,812.001,0.01)
S=np.array([[status(T0+n*P) for n in (-4,-3)] for P in Ps])
fig=plt.figure(figsize=(7.2,6.2))
ax=fig.add_axes([0.17,0.62,0.8,0.31])
cols={0:'#bbbbbb',1:'#f2a541',2:'#3b6ea5'}
for k,n in enumerate((-4,-3)):
    for code in (0,1,2):
        m=S[:,k]==code
        ax.scatter(Ps[m],np.full(m.sum(),1-k),c=cols[code],marker='|',s=120,lw=0.6)
surv=np.all(S!=2,axis=1)
ax.scatter(Ps[surv],np.full(surv.sum(),-1),c='#c0392b',marker='|',s=120,lw=0.8)
ax.axvspan(800.9-28.4,800.9+41.1,color='k',alpha=0.06,lw=0)
ax.axvline(781.63,color='k',ls='--',lw=0.8); ax.axvline(776.106,color='k',ls=':',lw=0.8)
ax.text(781.8,1.45,'781.6 d (ML)',fontsize=7); ax.text(770.3,1.45,'776.1 d = (TESS − D800)/4',fontsize=7)
ax.set_yticks([1,0,-1]); ax.set_yticklabels(['n = −4 (Q7–Q9)','n = −3 (Q15–Q17)','allowed'],fontsize=8)
ax.set_ylim(-1.6,1.7); ax.set_xlim(770,812); ax.set_xlabel('Trial period (d)',fontsize=8); ax.tick_params(labelsize=7)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=cols[2],label='excluded (quiet Kepler data)'),Patch(color=cols[1],label='inside a ≥0.8% dip'),Patch(color=cols[0],label='data gap')],fontsize=6.5,loc='upper right',ncol=3,frameon=False,bbox_to_anchor=(1,1.17))
for i,(P,n,title) in enumerate([(781.63,-4,'P = 781.63 d, n = −4'),(1345.648,-2,'Window 4, P = 1345.65 d, n = −2'),(776.30,-4,'P = 776.30 d, n = −4 (D800)'),(776.30,-3,'P = 776.30 d, n = −3 (D1500)')]):
    a=fig.add_axes([0.1+0.47*(i%2),0.07+0.24*(1-i//2),0.4,0.17])
    tc=T0+n*P; m=np.abs(t-tc)<3
    a.plot(t[m]-tc,f[m],'.',ms=1.2,color='0.3')
    if i<2:
        oot=m&(np.abs(t-tc)>0.6); c=np.polyfit(t[oot]-tc,f[oot],1)
        a.plot(t[m]-tc,np.polyval(c,t[m]-tc)-model(t[m]-tc),color='#c0392b',lw=1,label='expected 1.1% transit')
        if i==0: a.legend(fontsize=6,frameon=False,loc='lower left')
    else:
        a.axvspan(-T14/2,T14/2,color='#c0392b',alpha=0.15,lw=0); a.axvspan(-T23/2,T23/2,color='#c0392b',alpha=0.15,lw=0)
    a.set_title(title+f' (BKJD {tc:.1f})',fontsize=7.5); a.tick_params(labelsize=6.5)
    if i%2==0: a.set_ylabel('Rel. flux',fontsize=7)
    if i//2==1: a.set_xlabel('Days from predicted mid-transit',fontsize=7)
fig.savefig('../boyajian_alias_figure.png',dpi=200)
