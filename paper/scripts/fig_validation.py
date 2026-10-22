"""Technical validation of the final model (held out in space, 2000-2025): (a) density of observed vs held-out, (b) station r and KGE by region,
(c) skill by period, (d) 90 % interval coverage by predicted level (cross-calibrated by block)."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
INK='#3d3d3a'; MUTED='#77756e'; GRID='#e6e5df'; BLUE='#2a78d6'; ORANGE='#eb6834'
plt.rcParams.update({'font.size':7.5,'axes.edgecolor':GRID,'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.titlelocation':'left','axes.titlesize':8})
C=pd.read_parquet('work/final/cv_predictions_final.parquet'); C=C[C.salinity.notna()&~C.fresh_site]
ST=pd.read_csv('station_skill_final.csv'); ST=ST[~ST.fresh]; P=pd.read_csv('work/final/gauge/table4_skill_by_period_final.csv'); V=pd.read_csv('interval_coverage_crosscal.csv',index_col=0).iloc[:,0]
fig,ax=plt.subplots(2,2,figsize=(7.09,6.2)); plt.subplots_adjust(hspace=0.5,wspace=0.42)
a=ax[0,0]; h=a.hexbin(C.salinity,C.pred_cv,gridsize=70,extent=(0,60,0,60),bins='log',cmap='Blues',mincnt=1,lw=0)
a.plot([0,60],[0,60],color=INK,lw=0.8,ls='--'); e=C.pred_cv-C.salinity; r2=1-(e**2).sum()/((C.salinity-C.salinity.mean())**2).sum()
a.set_xlabel('gauge salinity (psu)'); a.set_ylabel('held-out model (psu)'); a.set_title(f'a  All held-out station-days (n = {len(C):,})')
a.text(2,55,f'R² = {r2:.2f}\nRMSE = {np.sqrt((e**2).mean()):.1f} psu\nbias = {e.mean():+.1f} psu',va='top',fontsize=7)
cb=fig.colorbar(h,ax=a,fraction=0.045,pad=0.02); cb.set_label('station-days',fontsize=7)
a=ax[0,1]; REG=['Shark River Slough','Taylor Slough','C-111','Florida Bay','Big Cypress']
for i,rg in enumerate(REG):
    s=ST[ST.region==rg]; x=np.full(len(s),i)+np.random.default_rng(i).uniform(-0.12,0.12,len(s))
    a.scatter(x-0.17,s.r,s=9,color=BLUE,lw=0,alpha=0.8,label='r' if i==0 else None); a.scatter(x+0.17,s.KGE,s=9,color=ORANGE,lw=0,alpha=0.8,label='KGE' if i==0 else None)
    a.plot([i-0.32,i-0.02],[s.r.median()]*2,color=INK,lw=1.2); a.plot([i+0.02,i+0.32],[s.KGE.median()]*2,color=INK,lw=1.2)
SH={'Shark River Slough':'Shark R.\nSlough','Taylor Slough':'Taylor\nSlough','C-111':'C-111\n','Florida Bay':'Florida\nBay','Big Cypress':'Big\nCypress'}
a.set_xticks(range(len(REG))); a.set_xticklabels([f'{SH[r]}\n({(ST.region==r).sum()})' for r in REG],fontsize=6.5)
a.set_ylim(0,1.02); a.legend(frameon=False,fontsize=6.5,loc='lower left',ncol=2); a.set_ylabel('per-station skill')
a.set_title('b  Station correlation and KGE by region (stations)')
a=ax[1,0]; x=np.arange(len(P)); a.bar(x,P.R2,color=BLUE,width=0.6)
for i,(r2,b) in enumerate(zip(P.R2,P.bias)): a.text(i,r2+0.01,f'{r2:.2f}\nbias '+(f'{b:+.1f}' if abs(b)>=0.05 else '0.0'),ha='center',va='bottom',fontsize=6.5)
a.set_xticks(x); a.set_xticklabels(['2000–03','2004–07','2008–\nMar 2013','Apr 2013–\n2024','2025'],fontsize=6.5); a.set_ylim(0,1.12); a.set_ylabel('R² (held out in space)')
a.set_title('c  Skill by period')
a=ax[1,1]; labs=['0–5','5–15','15–30','> 30']; vals=[V[f'bin{k}'] for k in range(4)]
a.bar(range(4),[100*v for v in vals],color=BLUE,width=0.6); a.axhline(90,color=INK,lw=0.8,ls='--')
for i,v in enumerate(vals): a.text(i,100*v-6,f'{100*v:.1f} %',ha='center',fontsize=6.5,color='white')
a.set_xticks(range(4)); a.set_xticklabels(labs); a.set_xlabel('predicted salinity (psu)'); a.set_ylabel('observations inside 90 % interval (%)'); a.set_ylim(0,105)
a.set_title('d  90 % interval coverage')
for a in ax.ravel(): [a.spines[s].set_visible(False) for s in ('top','right')]
fig.savefig('Fig03_validation.png',dpi=300,bbox_inches='tight',facecolor='white'); print('ok')
