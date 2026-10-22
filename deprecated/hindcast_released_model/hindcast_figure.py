import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
INK='#3d3d3a'; MUTED='#8a897f'; GRIDC='#e6e5df'; OBS='#e0592a'; MOD='#1f5fae'
plt.rcParams.update({'font.size':7.5,'axes.edgecolor':GRIDC,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'text.color':INK,'axes.titlecolor':INK,'axes.titlesize':8,'axes.titlelocation':'left','legend.frameon':False})
T=pd.read_parquet('hindcast_station_predictions.parquet'); o=T[T.salinity.notna()].copy(); o['y']=o.date.dt.year
A=o.groupby(['station','y']).agg(n=('salinity','size'),obs=('salinity','mean'),cv=('p_cv','mean')).reset_index(); A=A[A.n>=300]
yr=A.groupby('y').agg(obs=('obs','mean'),cv=('cv','mean'),k=('station','size'))
R=pd.read_csv('hindcast_skill_by_era.csv'); R=R[R.model.str.startswith('B')]
fig=plt.figure(figsize=(7.09,6.6)); gs=fig.add_gridspec(3,2,height_ratios=[1,1,1],hspace=0.55,wspace=0.28)
a=fig.add_subplot(gs[0,0]); a.axvspan(2000.5,2012.5,color='#f1f0ec',lw=0)
a.plot(yr.index,yr.obs,color=OBS,lw=1.6,marker='o',ms=3,label='Gauges'); a.plot(yr.index,yr.cv,color=MOD,lw=1.6,marker='s',ms=3,label='Model, station held out')
a.text(2006.5,a.get_ylim()[1]*0.98,'hindcast',ha='center',va='top',color=MUTED,fontsize=7); a.set_ylabel('Mean salinity (psu)'); a.grid(axis='y',color=GRIDC,lw=0.5)
a.legend(loc='upper left',bbox_to_anchor=(0,-0.13),ncol=2,fontsize=6.5); a.set_title('a  Annual mean, sonde stations (≥ 300 days/yr)')
for s in ('top','right'): a.spines[s].set_visible(False)
a=fig.add_subplot(gs[0,1]); x=np.arange(len(R))
a.bar(x,R.R2,color=[MOD if e not in ('2013-24','2025') else MUTED for e in R.era],width=0.6)
for i,(r2,b) in enumerate(zip(R.R2,R.bias)): a.text(i,r2+0.01,f'{r2:.2f}\nbias {b:+.1f}',ha='center',va='bottom',fontsize=6.5)
a.set_xticks(x); a.set_xticklabels(['2000–03','2004–07','2008–13','2013–24','2025']); a.set_ylim(0,1.05); a.set_ylabel('R² (held out in space)')
a.set_title('b  Skill by period: hindcast (blue) vs training years'); a.grid(axis='y',color=GRIDC,lw=0.5)
for s in ('top','right'): a.spines[s].set_visible(False)
for i,stn in enumerate(['ENPNR','ENPLB']):
    a=fig.add_subplot(gs[1+i,:]); g=o[o.station==stn].sort_values('date'); full=T[T.station==stn].set_index('date').sort_index()
    a.axvspan(pd.Timestamp('2000-04-01'),pd.Timestamp('2013-03-31'),color='#f1f0ec',lw=0)
    a.plot(full.index,full.p_cv,color=MOD,lw=0.6,label='Model, station held out')
    a.scatter(g.date,g.salinity,s=0.3,lw=0,color=OBS,label='Gauge (daily mean)',rasterized=True)
    pre=g[g.date<'2013-04-01']; y,p=pre.salinity.values,pre.p_cv.values
    a.set_title(f'{"cd"[i]}  {stn}   2000–13: r {np.corrcoef(y,p)[0,1]:.2f} · RMSE {np.sqrt(((y-p)**2).mean()):.1f} · bias {(p-y).mean():+.1f} psu',fontsize=7.5)
    a.set_ylabel('psu'); a.grid(axis='y',color=GRIDC,lw=0.5); a.set_xlim(pd.Timestamp('2000-01-01'),pd.Timestamp('2026-01-01'))
    for s in ('top','right'): a.spines[s].set_visible(False)
    if i==0: a.legend(loc='upper right',fontsize=6.5,markerscale=8,ncol=2)
fig.savefig('Fig_hindcast_2000_2012.png',dpi=300,bbox_inches='tight',facecolor='white'); print('ok')
