"""Daily out-of-block estimates (final blend = 0.7*two-step K1 + 0.3*v3) vs gauge salinity at every station."""
import pickle, numpy as np, pandas as pd, matplotlib, os; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.dates as mdates
d=pickle.load(open('best_two_step_oof.pkl','rb')); P=d['P']
W={'K1 two-step':0.7,'v3 one-step':0.3}
A=d['keyA'].copy(); B=d['keyB'].copy()
A['model']=sum(w*P[k][0] for k,w in W.items()); B['model']=sum(w*P[k][1] for k,w in W.items())
A['set']='2013–2024 (held-out block)'; B['set']='2025 (held-out block + year)'
D=pd.concat([A,B]).dropna(subset=['model']).sort_values(['station','date'])
def stats(g):
    y=g.salinity.values; p=g.model.values
    if len(y)<10: return pd.Series(dict(n=len(y),obs_mean=np.nan,model_mean=np.nan,bias=np.nan,rmse=np.nan,r=np.nan,nse=np.nan,kge=np.nan))
    r=np.corrcoef(y,p)[0,1] if y.std()>0 else np.nan
    kge=1-np.sqrt((r-1)**2+(p.std()/y.std()-1)**2+(p.mean()/y.mean()-1)**2) if y.std()>0 and y.mean()>0 else np.nan
    return pd.Series(dict(n=len(y),obs_mean=y.mean(),model_mean=p.mean(),bias=p.mean()-y.mean(),rmse=np.sqrt(((y-p)**2).mean()),r=r,
                          nse=1-((y-p)**2).sum()/((y-y.mean())**2).sum() if y.std()>0 else np.nan,kge=kge))
ST=D.groupby('station').apply(stats); ST['basin']=D.groupby('station').basin.first()
ST25=D[D.date.dt.year==2025].groupby('station').apply(stats)
ST=ST.join(ST25[['n','r','rmse','bias','kge']].add_suffix('_2025'))
ST.round(3).to_csv('station_timeseries_metrics.csv')
D[['station','date','basin','salinity','model']].to_csv('station_daily_model_vs_gauge.csv.gz',index=False)
ok=ST.n>=50
summ=dict(stations=int(ok.sum()),median_r=ST.r[ok].median(),median_kge=ST.kge[ok].median(),median_rmse=ST.rmse[ok].median(),median_abs_bias=ST.bias[ok].abs().median(),
          frac_r_gt_0_7=(ST.r[ok]>0.7).mean(),frac_kge_gt_0_5=(ST.kge[ok]>0.5).mean(),frac_kge_gt_0=(ST.kge[ok]>0).mean(),
          median_r_2025=ST.r_2025[ST.n_2025>=30].median(),median_kge_2025=ST.kge_2025[ST.n_2025>=30].median())
print({k:round(v,3) for k,v in summ.items()}); pd.Series(summ).to_csv('station_timeseries_summary.csv')
print(ST[ok].sort_values('kge').head(10)[['basin','n','obs_mean','bias','rmse','r','kge']].round(2))
# ---- figures
OBS='#eb6834'; MOD='#2a78d6'; INK='#3d3d3a'; GRID='#e6e5df'
os.makedirs('ts_figs',exist_ok=True)
plt.rcParams.update({'font.size':8,'axes.edgecolor':GRID,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'axes.titlecolor':INK})
for bas,g in ST[ok].sort_values('obs_mean').groupby('basin'):
    sts=list(g.index); nc=4; nr=int(np.ceil(len(sts)/nc))
    fig,axs=plt.subplots(nr,nc,figsize=(16,2.1*nr),sharex=True,squeeze=False)
    for ax,s in zip(axs.ravel(),sts):
        x=D[D.station==s]
        ax.scatter(x.date,x.salinity,s=2,color=OBS,lw=0,label='Gauge (daily mean)')
        # break the model line across gaps > 10 days
        xx=x.set_index('date').model; dd=xx.index.to_series().diff().dt.days; sparse=dd.median()>10; gap=(dd>10)&(~sparse); seg=gap.cumsum()
        if sparse: ax.scatter(xx.index,xx.values,s=5,color=MOD,lw=0,marker='D')
        for _,q in xx.groupby(seg.values): ax.plot(q.index,q.values,color=MOD,lw=0.5 if sparse else 0.9,label='Model (held out)' if _==0 else None)
        ax.axvspan(pd.Timestamp('2025-01-01'),pd.Timestamp('2025-12-31'),color=GRID,alpha=0.6,lw=0)
        m=ST.loc[s]; ax.set_title(f"{s}   r {m.r:.2f} · KGE {m.kge:.2f} · RMSE {m.rmse:.1f} · bias {m.bias:+.1f}",fontsize=7.5,loc='left')
        ax.grid(axis='y',color=GRID,lw=0.5); ax.tick_params(length=2)
        for sp in ('top','right'): ax.spines[sp].set_visible(False)
    for ax in axs.ravel()[len(sts):]: ax.axis('off')
    for ax in axs[:,0]: ax.set_ylabel('Salinity (psu)')
    axs[0,0].legend(loc='upper left',fontsize=7,frameon=False,markerscale=3)
    axs[-1,0].xaxis.set_major_locator(mdates.YearLocator(2)); axs[-1,0].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    fig.suptitle(f'{bas} basin — daily salinity, gauge vs model (stations held out in spatial-block CV; 2025 shaded = held-out year)',x=0.01,ha='left',fontsize=10,color=INK)
    fig.tight_layout(rect=(0,0,1,0.97)); fig.savefig(f'ts_figs/timeseries_{bas}.png',dpi=90); plt.close(fig)
# summary figure: pooled obs vs model + per-station r and KGE
fig,ax=plt.subplots(1,3,figsize=(15,4.3))
hb=ax[0].hexbin(D.salinity,D.model,gridsize=70,bins='log',cmap='Blues',mincnt=1,extent=(0,70,0,70)); ax[0].plot([0,70],[0,70],color=INK,lw=0.8)
ax[0].set_xlabel('Gauge salinity (psu)'); ax[0].set_ylabel('Model salinity (psu)'); ax[0].set_title('All held-out station-days (log count)',loc='left')
ax[1].hist(ST.r[ok],bins=np.arange(0,1.01,0.05),color=MOD,edgecolor='white'); ax[1].axvline(ST.r[ok].median(),color=INK,lw=1)
ax[1].set_xlabel('Within-station correlation r'); ax[1].set_ylabel('Stations'); ax[1].set_title(f'Timing: median r = {ST.r[ok].median():.2f}',loc='left')
ax[2].hist(ST.kge[ok].clip(-1,1),bins=np.arange(-1,1.01,0.1),color=MOD,edgecolor='white'); ax[2].axvline(ST.kge[ok].median(),color=INK,lw=1)
ax[2].set_xlabel('KGE (clipped at −1)'); ax[2].set_title(f'Overall skill: median KGE = {ST.kge[ok].median():.2f}',loc='left')
for a in ax:
    for sp in ('top','right'): a.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig('ts_figs/summary_model_vs_gauge.png',dpi=90)
print('figs done')
