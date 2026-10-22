"""Make every figure in the methods/results document (Figures 1-7) from saved outputs.

    python make_figures.py            # all figures
    python make_figures.py 4 5        # only Figures 4 and 5

Colour convention: wherever colour encodes salinity, red = more saline and blue = fresher
(maps, tables, anomalies, trends, and the ecosystem/region line colours, which are ordered fresh -> saline).

Inputs (paths relative to the working directory; see README):
  best_two_step_oof.pkl                       held-out CV predictions (Figs 1-3)
  daily_maps/median_int16.npy                 250 m daily maps (Fig 7)
  map30/monthly30_int16.npy, map30/C_int16.npy  30 m monthly maps and 30 m detail layers (Figs 5, 7)
  px30_index.npy                              30 m in-domain cell index
  change30/*                                  30 m analysis outputs (Figs 4-6)
  site_variables/ENP_WaterFrequency_L30.tif   30 m grid definition
Outputs: figures/fig1_... to figures/fig7_...png
"""
import sys, os, pickle, numpy as np, pandas as pd, rasterio, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.dates as mdates
from matplotlib.colors import TwoSlopeNorm, to_rgb

GRID_TIF=os.environ.get('GRID_TIF','/mnt/user-data/uploads/ENP_Landsat_Salinity_Model/site_variables/ENP_WaterFrequency_L30.tif')
D='change30'; OUT='figures'; os.makedirs(OUT,exist_ok=True)
H,W=3839,4171; D30=0.00027092271166751034
INK='#3d3d3a'; MUTED='#8a897f'; GRIDC='#e6e5df'
SAL='RdYlBu_r'      # salinity: blue fresh -> red saline
DIV='RdBu_r'        # anomalies/trends: red = saltier
# categorical colours ordered fresh -> saline (validated for CVD separation; direct labels/legends added)
RAMP=['#1f5fae','#56b1e3','#2f9e6b','#e7a33a','#e0592a','#a3172c']
ECO=['Freshwater marsh & prairie','Swamp forest & scrub','Salt marsh','Mangrove','Inland & estuarine water','Open bay']
REG=['Big Cypress','Shark River Slough','Taylor Slough','C-111','Gulf coastal waters','Florida Bay']
REG_BY_SAL=['Taylor Slough','Shark River Slough','Big Cypress','C-111','Gulf coastal waters','Florida Bay']   # mean 7.3 ... 34.5 psu
COL_E=dict(zip(ECO,RAMP)); COL_R=dict(zip(REG_BY_SAL,RAMP))
plt.rcParams.update({'font.size':9,'axes.edgecolor':GRIDC,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'axes.titlecolor':INK,'text.color':INK})

def clean(a):
    for s in ('top','right'): a.spines[s].set_visible(False)
def ink_on(cmap,norm,v):
    r,g,b=to_rgb(plt.get_cmap(cmap)(norm(v))); return 'white' if 0.299*r+0.587*g+0.114*b<0.5 else INK

_cache={}
def base():
    if 'ext' not in _cache:
        wfr=rasterio.open(GRID_TIF); b=wfr.bounds
        _cache['ext']=[b.left,b.right,b.bottom,b.top]; _cache['idx']=np.load('px30_index.npy')
        reg=grid(np.load(f'{D}/region30.npy').astype('float32')); reg[reg<0]=np.nan; _cache['reg']=reg
    return _cache['ext']
def grid(v):
    idx=_cache.get('idx');
    if idx is None: idx=np.load('px30_index.npy'); _cache['idx']=idx
    a=np.full(H*W,np.nan,'float32'); a[idx]=v; return a.reshape(H,W)
def mapax(a):
    ext=base(); a.set_xticks([]); a.set_yticks([]); [a.spines[s].set_visible(False) for s in a.spines]
    for k in range(6):
        m=(_cache['reg']==k).astype('float32')[::4,::4]; a.contour(m,levels=[0.5],colors=INK,linewidths=0.5,extent=tuple(ext),origin='upper')

# ---------------------------------------------------------------- Figures 1-3: validation
def validation_data():
    d=pickle.load(open('best_two_step_oof.pkl','rb')); P=d['P']; Wt={'K1 two-step':0.7,'v3 one-step':0.3}
    A=d['keyA'].copy(); B=d['keyB'].copy()
    A['model']=sum(w*P[k][0] for k,w in Wt.items()); B['model']=sum(w*P[k][1] for k,w in Wt.items())
    Dv=pd.concat([A,B]).dropna(subset=['model']).sort_values(['station','date'])
    def stats(g):
        y=g.salinity.values; p=g.model.values
        if len(y)<10 or y.std()==0: return pd.Series(dict(n=len(y),obs_mean=y.mean(),bias=np.nan,rmse=np.nan,r=np.nan,kge=np.nan))
        r=np.corrcoef(y,p)[0,1]; kge=1-np.sqrt((r-1)**2+(p.std()/y.std()-1)**2+(p.mean()/y.mean()-1)**2) if y.mean()>0 else np.nan
        return pd.Series(dict(n=len(y),obs_mean=y.mean(),bias=p.mean()-y.mean(),rmse=np.sqrt(((y-p)**2).mean()),r=r,kge=kge))
    ST=Dv.groupby('station').apply(stats); ST['basin']=Dv.groupby('station').basin.first()
    return Dv,ST

OBS='#e0592a'; MOD='#1f5fae'
def fig1(Dv,ST):
    ok=ST.n>=50
    fig,ax=plt.subplots(1,3,figsize=(15,4.3))
    ax[0].hexbin(Dv.salinity,Dv.model,gridsize=70,bins='log',cmap='Greys',mincnt=1,extent=(0,70,0,70)); ax[0].plot([0,70],[0,70],color=OBS,lw=0.9)
    ax[0].set_xlabel('Gauge salinity (psu)'); ax[0].set_ylabel('Model salinity (psu)'); ax[0].set_title('a  All held-out station-days (log count)',loc='left')
    ax[1].hist(ST.r[ok],bins=np.arange(0,1.01,0.05),color=MOD,edgecolor='white'); ax[1].axvline(ST.r[ok].median(),color=INK,lw=1)
    ax[1].set_xlabel('Within-station correlation r'); ax[1].set_ylabel('Stations'); ax[1].set_title(f'b  Timing: median r = {ST.r[ok].median():.2f}',loc='left')
    ax[2].hist(ST.kge[ok].clip(-1,1),bins=np.arange(-1,1.01,0.1),color=MOD,edgecolor='white'); ax[2].axvline(ST.kge[ok].median(),color=INK,lw=1)
    ax[2].set_xlabel('KGE (clipped at −1)'); ax[2].set_title(f'c  Overall skill: median KGE = {ST.kge[ok].median():.2f}',loc='left')
    for a in ax: clean(a)
    fig.tight_layout(); fig.savefig(f'{OUT}/fig1_model_vs_gauge_summary.png',dpi=110); plt.close(fig)

def fig_ts(Dv,ST,basin,num,label):
    g=ST[(ST.n>=50)&(ST.basin==basin)].sort_values('obs_mean'); sts=list(g.index); nc=4; nr=int(np.ceil(len(sts)/nc))
    fig,axs=plt.subplots(nr,nc,figsize=(16,2.1*nr),sharex=True,squeeze=False)
    for ax,s in zip(axs.ravel(),sts):
        x=Dv[Dv.station==s]
        ax.scatter(x.date,x.salinity,s=2,color=OBS,lw=0,label='Gauge (daily mean)')
        xx=x.set_index('date').model; dd=xx.index.to_series().diff().dt.days; sparse=dd.median()>10; seg=((dd>10)&(~sparse)).cumsum()
        if sparse: ax.scatter(xx.index,xx.values,s=5,color=MOD,lw=0,marker='D')
        for i,(_,q) in enumerate(xx.groupby(seg.values)): ax.plot(q.index,q.values,color=MOD,lw=0.5 if sparse else 0.9,label='Model (held out)' if i==0 else None)
        ax.axvspan(pd.Timestamp('2025-01-01'),pd.Timestamp('2025-12-31'),color=GRIDC,alpha=0.6,lw=0)
        m=ST.loc[s]; ax.set_title(f"{s}   r {m.r:.2f} · KGE {m.kge:.2f} · RMSE {m.rmse:.1f} · bias {m.bias:+.1f}",fontsize=7.5,loc='left')
        ax.grid(axis='y',color=GRIDC,lw=0.5); ax.tick_params(length=2); clean(ax)
    for ax in axs.ravel()[len(sts):]: ax.axis('off')
    for ax in axs[:,0]: ax.set_ylabel('Salinity (psu)')
    axs[0,0].legend(loc='upper left',fontsize=7,frameon=False,markerscale=3)
    axs[-1,0].xaxis.set_major_locator(mdates.YearLocator(2)); axs[-1,0].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    fig.suptitle(f'{label} basin: daily salinity, gauge vs model (stations held out in spatial-block CV; 2025 shaded = held-out year)',x=0.01,ha='left',fontsize=10)
    fig.tight_layout(rect=(0,0,1,0.97)); fig.savefig(f'{OUT}/fig{num}_timeseries_{basin}.png',dpi=90); plt.close(fig)

# ---------------------------------------------------------------- Figure 4: regions and ecosystems
def fig4():
    ext=base(); R=np.load(f'{D}/pixel_results30.npz')
    T1=pd.read_csv(f'{D}/mean_salinity_region_x_ecosystem.csv',index_col=0).reindex(index=ECO,columns=REG); A1=pd.read_csv(f'{D}/area_km2_region_x_ecosystem.csv',index_col=0).reindex(index=ECO,columns=REG)
    CL=pd.read_csv(f'{D}/monthly_climatology_by_ecosystem.csv'); MS=pd.read_csv(f'{D}/monthly_series_by_ecosystem.csv',index_col=0,parse_dates=True)
    fig=plt.figure(figsize=(15,10.5)); gs=fig.add_gridspec(2,3,width_ratios=[1.05,1,1],hspace=0.32,wspace=0.28)
    a=fig.add_subplot(gs[0,0]); im=a.imshow(grid(R['mean']),extent=ext,cmap=SAL,vmin=0,vmax=40,interpolation='antialiased'); mapax(a)
    for nm,(x,y) in {'Big Cypress':(-81.33,25.72),'Shark River\nSlough':(-80.95,25.62),'Taylor\nSlough':(-80.66,25.40),'C-111':(-80.47,25.33),'Gulf\nwaters':(-81.40,25.35),'Florida Bay':(-80.80,24.98)}.items():
        a.text(x,y,nm,fontsize=8,ha='center',va='center',bbox=dict(boxstyle='round,pad=0.15',fc='white',ec='none',alpha=0.8))
    fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal').set_label('Mean salinity 2013–2025 (psu), 30 m'); a.set_title('a  Mean salinity and regions',loc='left',fontsize=10)
    a=fig.add_subplot(gs[0,1:]); v=np.where(A1.values>=5,T1.values,np.nan); norm=plt.Normalize(0,40); a.imshow(v,cmap=SAL,norm=norm,aspect='auto')
    a.set_xticks(range(len(REG))); a.set_xticklabels([r.replace(' ','\n',1) for r in REG]); a.set_yticks(range(len(ECO))); a.set_yticklabels(ECO)
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            if np.isfinite(v[i,j]): a.text(j,i,f'{v[i,j]:.1f}\n{A1.values[i,j]:,.0f} km²',ha='center',va='center',fontsize=8,color=ink_on(SAL,norm,v[i,j]))
            else: a.text(j,i,'–',ha='center',va='center',fontsize=8,color=MUTED)
    a.set_title('b  Mean salinity (psu) and area by region × ecosystem type (30 m cells)',loc='left',fontsize=10); [a.spines[s].set_visible(False) for s in a.spines]; a.tick_params(length=0)
    a=fig.add_subplot(gs[1,0])
    for e in ECO:
        c=CL[CL.ecosystem==e]; a.plot(c.month,c['mean'],color=COL_E[e],lw=2,marker='o',ms=4); a.text(12.3,c['mean'].iloc[-1],e,fontsize=7.5,va='center')
    a.set_xticks(range(1,13)); a.set_xticklabels(list('JFMAMJJASOND')); a.set_xlim(0.6,16.5); a.set_ylabel('Salinity (psu)'); a.grid(axis='y',color=GRIDC,lw=0.6); clean(a)
    a.set_title('c  Seasonal cycle by ecosystem (2013–2025)',loc='left',fontsize=10)
    a=fig.add_subplot(gs[1,1:])
    for e in ECO: a.plot(MS.index,MS[e],color=COL_E[e],lw=1.4,label=e)
    a.set_ylabel('Salinity (psu)'); a.grid(axis='y',color=GRIDC,lw=0.6); clean(a); a.legend(frameon=False,fontsize=8,ncol=6,loc='upper center',bbox_to_anchor=(0.5,-0.08))
    for y0,y1,lab in [('2015-01-01','2015-08-31','2015 drought'),('2025-01-01','2025-06-30','2025 dry')]:
        a.axvspan(pd.Timestamp(y0),pd.Timestamp(y1),color=GRIDC,alpha=0.7,lw=0); a.text(pd.Timestamp(y0),43,lab,fontsize=7.5,color=MUTED,va='bottom')
    a.set_ylim(-1,45); a.set_title('d  Monthly mean salinity by ecosystem, 2013–2025',loc='left',fontsize=10)
    fig.savefig(f'{OUT}/fig4_salinity_regions_ecosystems_30m.png',dpi=130,bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- Figure 5: change
def fig5():
    ext=base(); R=np.load(f'{D}/pixel_results30.npz'); TM=np.load(f'{D}/pixel_trends_monthly_anomaly30.npz')
    MO=np.load('map30/monthly30_int16.npy',mmap_mode='r'); yr=np.repeat(np.arange(2013,2026),12)
    def annual(y):
        x=np.asarray(MO[yr==y],dtype='float32')/10; x[x<0]=np.nan; return np.nanmean(x,0)
    mean=grid(R['mean']); an15=annual(2015)-mean; an24=annual(2024)-mean
    TS=pd.read_csv(f'{D}/annual_seasonal_means_by_group.csv')
    fig=plt.figure(figsize=(16,10.5)); gs=fig.add_gridspec(2,4,height_ratios=[1,0.85],hspace=0.28,wspace=0.1)
    a=fig.add_subplot(gs[0,0]); im=a.imshow(grid(R['sd_annual']),extent=ext,cmap='Purples',vmin=0,vmax=6); mapax(a)
    fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal').set_label('SD of annual means (psu)'); a.set_title('a  Year-to-year variability',loc='left',fontsize=10)
    for k,(arr,t,y) in enumerate([(an15,'b  2015 drought',2015),(an24,'c  2024 wet year',2024)]):
        a=fig.add_subplot(gs[0,k+1]); im=a.imshow(arr,extent=ext,cmap=DIV,norm=TwoSlopeNorm(0,-10,10)); mapax(a)
        fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal').set_label(f'{y} anomaly vs 2013–2025 mean (psu; red = saltier)'); a.set_title(t,loc='left',fontsize=10)
    a=fig.add_subplot(gs[0,3]); im=a.imshow(grid(TM['slope_monthly']),extent=ext,cmap=DIV,norm=TwoSlopeNorm(0,-0.5,0.5)); mapax(a)
    fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal').set_label('Trend in monthly anomalies (psu yr⁻¹)')
    nsig=int(np.nansum(TM['q_monthly']<0.05)); a.set_title(f'd  Monthly-anomaly trend ({"none" if nsig==0 else f"{nsig:,} cells"} significant)',loc='left',fontsize=10)
    for k,(sea,t) in enumerate([('dry','e  Dry-season (Dec–May) anomaly'),('wet','f  Wet-season (Jun–Nov) anomaly')]):
        a=fig.add_subplot(gs[1,2*k:2*k+2]); an=TS[(TS.group=='ecosystem')&(TS.season==sea)].pivot(index='name',columns='year',values='mean').reindex(ECO); an=an.sub(an.mean(axis=1),axis=0)
        norm=TwoSlopeNorm(0,-6,6); a.imshow(an.values,cmap=DIV,norm=norm,aspect='auto')
        for i in range(an.shape[0]):
            for j in range(an.shape[1]): v=an.values[i,j]; a.text(j,i,f'{v:+.1f}',ha='center',va='center',fontsize=7,color=ink_on(DIV,norm,v))
        a.set_yticks(range(len(ECO))); a.set_yticklabels(ECO if k==0 else ['']*len(ECO),fontsize=8.5); a.set_xticks(range(an.shape[1])); a.set_xticklabels(an.columns,fontsize=8)
        [a.spines[s_].set_visible(False) for s_ in a.spines]; a.tick_params(length=0); a.set_title(t+' by ecosystem (psu vs 2013–2025 mean)',loc='left',fontsize=10)
    fig.savefig(f'{OUT}/fig5_salinity_change_30m.png',dpi=130,bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- Figure 6: high-salinity events
def fig6():
    ext=base(); EV=np.load(f'{D}/pixel_events30.npz')
    RE=pd.read_csv(f'{D}/high_salinity_events_by_region.csv').set_index('group').reindex(REG_BY_SAL); Y=pd.read_csv(f'{D}/high_salinity_events_by_region_year.csv'); yrs=np.arange(2013,2026)
    ED=Y.pivot(index='group',columns='year',values='event_days_per_pixel').reindex(index=REG_BY_SAL,columns=yrs); HY=Y.pivot(index='group',columns='year',values='hyper_days_mean').reindex(index=REG_BY_SAL,columns=yrs)
    GV=pd.read_csv(f'{D}/event_days_gauge_vs_model_30m.csv'); V=GV.groupby('year')[['obs','mod']].mean()
    fig=plt.figure(figsize=(16,11.5)); gs=fig.add_gridspec(2,3,hspace=0.32,wspace=0.25)
    a=fig.add_subplot(gs[0,0]); im=a.imshow(grid(EV['mean_excess']),extent=ext,cmap='Reds',vmin=0,vmax=20); mapax(a)
    fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal').set_label('Mean event peak above local median (psu)'); a.set_title('a  Event magnitude',loc='left',fontsize=10)
    a=fig.add_subplot(gs[0,1]); im=a.imshow(grid(EV['hyper_days_per_yr']),extent=ext,cmap='Reds',vmin=0,vmax=20); mapax(a)
    fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal').set_label('Days per year above 40 psu (2013–2025 mean)'); a.set_title('b  Hypersaline days',loc='left',fontsize=10)
    a=fig.add_subplot(gs[0,2]); y=np.arange(len(REG_BY_SAL))[::-1]
    for yi,r in zip(y,REG_BY_SAL):
        a.plot([RE.loc[r,'median_excess_psu'],RE.loc[r,'p90_excess_psu']],[yi+0.15]*2,color=COL_R[r],lw=2); a.plot(RE.loc[r,'median_excess_psu'],yi+0.15,'o',color=COL_R[r],ms=7)
        a.plot([RE.loc[r,'median_peak_psu'],RE.loc[r,'p90_peak_psu']],[yi-0.15]*2,color=COL_R[r],lw=2,alpha=0.5); a.plot(RE.loc[r,'median_peak_psu'],yi-0.15,'s',color=COL_R[r],ms=6,alpha=0.6)
    a.set_yticks(y); a.set_yticklabels(REG_BY_SAL); a.set_xlabel('psu\n● peak above local median   ■ peak salinity   (line: median to 90th pct)',fontsize=7.5); a.grid(axis='x',color=GRIDC,lw=0.6); clean(a)
    a.set_title('c  Event magnitude by region',loc='left',fontsize=10)
    a=fig.add_subplot(gs[1,0:2])
    for r in REG_BY_SAL: a.plot(yrs,ED.loc[r].values,color=COL_R[r],lw=2,marker='o',ms=4,label=r)
    a.set_xlim(2012.6,2025.4); a.set_xticks(yrs); a.set_ylabel('Event days per pixel per year'); a.grid(axis='y',color=GRIDC,lw=0.6); clean(a)
    a.legend(frameon=False,fontsize=8,ncol=6,loc='upper center',bbox_to_anchor=(0.5,-0.07))
    ax2=a.inset_axes([0.47,0.52,0.3,0.44]); ax2.plot(V.index,V.obs,color=INK,lw=1.5,marker='o',ms=3,label='gauges (observed)'); ax2.plot(V.index,V['mod'],color=OBS,lw=1.5,marker='s',ms=3,label='30 m model at gauges')
    ax2.set_title(f'Check: event days per station-year, {GV.station.nunique()} gauges',fontsize=7.5,loc='left'); ax2.tick_params(labelsize=7); ax2.legend(fontsize=7,frameon=False); clean(ax2)
    a.set_title('d  Frequency: days in high-salinity events (≥ 7 d above the local 95th percentile)',loc='left',fontsize=10)
    a=fig.add_subplot(gs[1,2]); norm=plt.Normalize(0,60); a.imshow(HY.values,cmap='Reds',norm=norm,aspect='auto')
    for i in range(HY.shape[0]):
        for j in range(HY.shape[1]): v=HY.values[i,j]; a.text(j,i,f'{v:.0f}' if v>=0.5 else '·',ha='center',va='center',fontsize=7,color=ink_on('Reds',norm,v))
    a.set_yticks(range(len(REG_BY_SAL))); a.set_yticklabels(REG_BY_SAL,fontsize=8); a.set_xticks(range(len(yrs))); a.set_xticklabels([str(y)[2:] for y in yrs],fontsize=8)
    [a.spines[s].set_visible(False) for s in a.spines]; a.tick_params(length=0); a.set_title("e  Days above 40 psu per pixel ('13–'25)",loc='left',fontsize=10)
    fig.savefig(f'{OUT}/fig6_high_salinity_events_30m.png',dpi=130,bbox_inches='tight'); plt.close(fig)

# ---------------------------------------------------------------- Figure 7: 250 m vs 30 m
def fig7():
    ext=base(); b=ext; m=12*7+4
    MED=np.load('daily_maps/median_int16.npy',mmap_mode='r'); ed=pd.date_range('2013-01-01','2025-12-31')
    dd=np.flatnonzero((ed>='2020-05-01')&(ed<='2020-05-31')); g=np.asarray(MED[dd[0]:dd[-1]+1],dtype='float32')/10; g[g<0]=np.nan
    with np.errstate(all='ignore'): g250=np.nanmean(g,0)
    MO=np.load('map30/monthly30_int16.npy',mmap_mode='r'); m30=np.asarray(MO[m],dtype='float32')/10; m30[m30<0]=np.nan
    C=np.load('map30/C_int16.npy',mmap_mode='r'); c=np.asarray(C[m],dtype='float32'); c[c==-32768]=np.nan; c/=10
    ext250=[b[0]-3.5*D30,b[0]+(521*8-3.5)*D30,b[3]-(480*8-3.5)*D30,b[3]+3.5*D30]; box=(-81.20,-80.90,25.28,25.47)
    fig,ax=plt.subplots(1,3,figsize=(15,4.6),constrained_layout=True)
    im=ax[0].imshow(g250,extent=ext250,cmap=SAL,vmin=0,vmax=40,interpolation='nearest'); ax[0].set_title('a  250 m map, May 2020 mean',loc='left',fontsize=10)
    ax[1].imshow(m30,extent=ext,cmap=SAL,vmin=0,vmax=40,interpolation='nearest'); ax[1].set_title('b  30 m map, May 2020 mean',loc='left',fontsize=10)
    im2=ax[2].imshow(c,extent=ext,cmap=DIV,norm=TwoSlopeNorm(0,-3,3),interpolation='nearest'); ax[2].set_title('c  30 m detail: 30 m model − 250 m map, 15 May 2020',loc='left',fontsize=10)
    for a in ax: a.set_xlim(box[0],box[1]); a.set_ylim(box[2],box[3]); a.set_xticks([]); a.set_yticks([]); a.set_aspect(1/np.cos(np.deg2rad(25.37)))
    fig.colorbar(im,ax=ax[:2],shrink=0.85,pad=0.01).set_label('Salinity (psu)'); fig.colorbar(im2,ax=ax[2],shrink=0.85,pad=0.01).set_label('psu (red = saltier at 30 m)')
    fig.savefig(f'{OUT}/fig7_250m_vs_30m.png',dpi=140,bbox_inches='tight'); plt.close(fig)

if __name__=='__main__':
    want={int(a) for a in sys.argv[1:]} or set(range(1,8))
    if want&{1,2,3}:
        Dv,ST=validation_data()
        if 1 in want: fig1(Dv,ST)
        if 2 in want: fig_ts(Dv,ST,'TaylorSlough',2,'Taylor Slough')
        if 3 in want: fig_ts(Dv,ST,'SRS',3,'Shark River Slough')
    for k,f in {4:fig4,5:fig5,6:fig6,7:fig7}.items():
        if k in want: f()
    print('figures written to',OUT,sorted(want))
