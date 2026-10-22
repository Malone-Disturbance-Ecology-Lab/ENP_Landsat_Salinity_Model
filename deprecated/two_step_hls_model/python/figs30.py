"""Figures 4-7 from the 30 m products."""
import numpy as np, pandas as pd, rasterio, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
D='change30'; H,W=3839,4171
wfr=rasterio.open('/mnt/user-data/uploads/ENP_Landsat_Salinity_Model/site_variables/ENP_WaterFrequency_L30.tif'); b=wfr.bounds; ext=[b.left,b.right,b.bottom,b.top]
idx=np.load('px30_index.npy'); R=np.load(f'{D}/pixel_results30.npz'); EV=np.load(f'{D}/pixel_events30.npz')
def grid(v):
    a=np.full(H*W,np.nan,'float32'); a[idx]=v; return a.reshape(H,W)
reg=grid(np.load(f'{D}/region30.npy').astype('float32')); reg[reg<0]=np.nan
ECO=['Freshwater marsh & prairie','Swamp forest & scrub','Salt marsh','Mangrove','Inland & estuarine water','Open bay']
REG=['Big Cypress','Shark River Slough','Taylor Slough','C-111','Gulf coastal waters','Florida Bay']
COL={'Open bay':'#2a78d6','Mangrove':'#eb6834','Inland & estuarine water':'#1baf7a','Salt marsh':'#eda100','Swamp forest & scrub':'#e87ba4','Freshwater marsh & prairie':'#008300'}
CR=dict(zip(REG,['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300']))
INK='#3d3d3a'; MUTED='#8a897f'; GRID='#e6e5df'
plt.rcParams.update({'font.size':9,'axes.edgecolor':GRID,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'axes.titlecolor':INK})
def mapax(a,outl=True):
    a.set_xticks([]); a.set_yticks([]); [a.spines[s].set_visible(False) for s in a.spines]
    if outl:
        for k in range(6):
            m=(reg==k).astype('float32')[::4,::4]; a.contour(m,levels=[0.5],colors=INK,linewidths=0.5,extent=(ext[0],ext[1],ext[2],ext[3]),origin='upper')
def clean(a):
    for s in ('top','right'): a.spines[s].set_visible(False)
AREA30=(0.00027092271*111.0)**2*np.cos(np.deg2rad(25.3))
# ================= Figure 4
T1=pd.read_csv(f'{D}/mean_salinity_region_x_ecosystem.csv',index_col=0).reindex(index=ECO,columns=REG); A1=pd.read_csv(f'{D}/area_km2_region_x_ecosystem.csv',index_col=0).reindex(index=ECO,columns=REG)
CL=pd.read_csv(f'{D}/monthly_climatology_by_ecosystem.csv'); MS=pd.read_csv(f'{D}/monthly_series_by_ecosystem.csv',index_col=0,parse_dates=True)
fig=plt.figure(figsize=(15,10.5)); gs=fig.add_gridspec(2,3,width_ratios=[1.05,1,1],hspace=0.32,wspace=0.28)
a=fig.add_subplot(gs[0,0]); im=a.imshow(grid(R['mean']),extent=ext,cmap='Blues',vmin=0,vmax=40,interpolation='antialiased'); mapax(a)
for nm,(x,y) in {'Big Cypress':(-81.33,25.72),'Shark River\nSlough':(-80.95,25.62),'Taylor\nSlough':(-80.66,25.40),'C-111':(-80.47,25.33),'Gulf\nwaters':(-81.40,25.35),'Florida Bay':(-80.80,24.98)}.items():
    a.text(x,y,nm,fontsize=8,color=INK,ha='center',va='center',bbox=dict(boxstyle='round,pad=0.15',fc='white',ec='none',alpha=0.75))
cb=fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal'); cb.set_label('Mean salinity 2013–2025 (psu), 30 m'); a.set_title('a  Mean salinity and regions',loc='left',fontsize=10)
a=fig.add_subplot(gs[0,1:]); v=np.where(A1.values>=5,T1.values,np.nan); a.imshow(v,cmap='Blues',vmin=0,vmax=40,aspect='auto')
a.set_xticks(range(len(REG))); a.set_xticklabels([r.replace(' ','\n',1) for r in REG]); a.set_yticks(range(len(ECO))); a.set_yticklabels(ECO)
for i in range(v.shape[0]):
    for j in range(v.shape[1]):
        a.text(j,i,f'{v[i,j]:.1f}\n{A1.values[i,j]:,.0f} km²' if np.isfinite(v[i,j]) else '–',ha='center',va='center',fontsize=8,color=('white' if np.isfinite(v[i,j]) and v[i,j]>22 else (INK if np.isfinite(v[i,j]) else MUTED)))
a.set_title('b  Mean salinity (psu) and area by region × ecosystem type (30 m cells)',loc='left',fontsize=10); [a.spines[s].set_visible(False) for s in a.spines]; a.tick_params(length=0)
a=fig.add_subplot(gs[1,0])
for e in ECO:
    c=CL[CL.ecosystem==e]; a.plot(c.month,c['mean'],color=COL[e],lw=2,marker='o',ms=4); a.text(12.3,c['mean'].iloc[-1],e,color=INK,fontsize=7.5,va='center')
a.set_xticks(range(1,13)); a.set_xticklabels(list('JFMAMJJASOND')); a.set_xlim(0.6,16.5); a.set_ylabel('Salinity (psu)'); a.grid(axis='y',color=GRID,lw=0.6); clean(a); a.set_title('c  Seasonal cycle by ecosystem (2013–2025)',loc='left',fontsize=10)
a=fig.add_subplot(gs[1,1:])
for e in ECO: a.plot(MS.index,MS[e],color=COL[e],lw=1.4,label=e)
a.set_ylabel('Salinity (psu)'); a.grid(axis='y',color=GRID,lw=0.6); clean(a); a.legend(frameon=False,fontsize=8,ncol=6,loc='upper center',bbox_to_anchor=(0.5,-0.08))
for y0,y1,lab in [('2015-01-01','2015-08-31','2015 drought'),('2025-01-01','2025-06-30','2025 dry')]:
    a.axvspan(pd.Timestamp(y0),pd.Timestamp(y1),color=GRID,alpha=0.7,lw=0); a.text(pd.Timestamp(y0),-1.5,lab,fontsize=7.5,color=MUTED,va='top')
a.set_title('d  Monthly mean salinity by ecosystem, 2013–2025',loc='left',fontsize=10)
fig.savefig(f'{D}/fig4_salinity_regions_ecosystems_30m.png',dpi=130,bbox_inches='tight'); plt.close(fig)
# ================= Figure 5
MO=np.load('map30/monthly30_int16.npy',mmap_mode='r'); yr=np.repeat(np.arange(2013,2026),12)
def annual(y):
    x=np.asarray(MO[yr==y],dtype='float32')/10; x[x<0]=np.nan; return np.nanmean(x,0)
mean=grid(R['mean']); an15=annual(2015)-mean; an24=annual(2024)-mean
TS=pd.read_csv(f'{D}/annual_seasonal_means_by_group.csv'); div=plt.get_cmap('RdBu_r')
fig=plt.figure(figsize=(16,10.5)); gs=fig.add_gridspec(2,4,height_ratios=[1,0.85],hspace=0.28,wspace=0.1)
a=fig.add_subplot(gs[0,0]); im=a.imshow(grid(R['sd_annual']),extent=ext,cmap='Oranges',vmin=0,vmax=6); mapax(a); cb=fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal'); cb.set_label('SD of annual means (psu)'); a.set_title('a  Year-to-year variability',loc='left',fontsize=10)
for k,(arr,t,y) in enumerate([(an15,'b  2015 drought',2015),(an24,'c  2024 wet year',2024)]):
    a=fig.add_subplot(gs[0,k+1]); im=a.imshow(arr,extent=ext,cmap=div,norm=TwoSlopeNorm(0,-10,10)); mapax(a); cb=fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal'); cb.set_label(f'{y} anomaly vs 2013–2025 mean (psu)'); a.set_title(t,loc='left',fontsize=10)
a=fig.add_subplot(gs[0,3]); im=a.imshow(grid(R['slope_annual']),extent=ext,cmap=div,norm=TwoSlopeNorm(0,-1,1)); mapax(a); cb=fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal'); cb.set_label('Theil–Sen trend (psu yr⁻¹)')
nsig=int(np.nansum(R['q_annual']<0.05)); a.set_title(f'd  2013–2025 trend ({"none" if nsig==0 else f"{nsig:,} cells"} significant)',loc='left',fontsize=10)
for k,(sea,t) in enumerate([('dry','e  Dry-season (Dec–May) anomaly'),('wet','f  Wet-season (Jun–Nov) anomaly')]):
    a=fig.add_subplot(gs[1,2*k:2*k+2]); an=TS[(TS.group=='ecosystem')&(TS.season==sea)].pivot(index='name',columns='year',values='mean').reindex(ECO); an=an.sub(an.mean(axis=1),axis=0)
    a.imshow(an.values,cmap=div,norm=TwoSlopeNorm(0,-6,6),aspect='auto')
    for i in range(an.shape[0]):
        for j in range(an.shape[1]): v=an.values[i,j]; a.text(j,i,f'{v:+.1f}',ha='center',va='center',fontsize=7,color='white' if abs(v)>3.5 else INK)
    a.set_yticks(range(len(ECO))); a.set_yticklabels(ECO if k==0 else ['']*len(ECO),fontsize=8.5); a.set_xticks(range(an.shape[1])); a.set_xticklabels(an.columns,fontsize=8)
    [a.spines[s_].set_visible(False) for s_ in a.spines]; a.tick_params(length=0); a.set_title(t+' by ecosystem (psu vs 2013–2025 mean)',loc='left',fontsize=10)
fig.savefig(f'{D}/fig5_salinity_change_30m.png',dpi=130,bbox_inches='tight'); plt.close(fig)
# ================= Figure 6
RE=pd.read_csv(f'{D}/high_salinity_events_by_region.csv').set_index('group').reindex(REG); Y=pd.read_csv(f'{D}/high_salinity_events_by_region_year.csv'); yrs=np.arange(2013,2026)
ED=Y.pivot(index='group',columns='year',values='event_days_per_pixel').reindex(index=REG,columns=yrs); HY=Y.pivot(index='group',columns='year',values='hyper_days_mean').reindex(index=REG,columns=yrs)
V=pd.read_csv(f'{D}/event_days_gauge_vs_model_30m.csv').groupby('year')[['obs','mod']].mean()
fig=plt.figure(figsize=(16,11.5)); gs=fig.add_gridspec(2,3,hspace=0.32,wspace=0.25)
a=fig.add_subplot(gs[0,0]); im=a.imshow(grid(EV['mean_excess']),extent=ext,cmap='Reds',vmin=0,vmax=20); mapax(a); cb=fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal'); cb.set_label('Mean event peak above local median (psu)'); a.set_title('a  Event magnitude',loc='left',fontsize=10)
a=fig.add_subplot(gs[0,1]); im=a.imshow(grid(EV['hyper_days_per_yr']),extent=ext,cmap='Reds',vmin=0,vmax=20); mapax(a); cb=fig.colorbar(im,ax=a,shrink=0.8,pad=0.03,orientation='horizontal'); cb.set_label('Days per year above 40 psu (2013–2025 mean)'); a.set_title('b  Hypersaline days',loc='left',fontsize=10)
a=fig.add_subplot(gs[0,2]); y=np.arange(len(REG))[::-1]
for yi,r in zip(y,REG):
    a.plot([RE.loc[r,'median_excess_psu'],RE.loc[r,'p90_excess_psu']],[yi+0.15]*2,color=CR[r],lw=2); a.plot(RE.loc[r,'median_excess_psu'],yi+0.15,'o',color=CR[r],ms=7)
    a.plot([RE.loc[r,'median_peak_psu'],RE.loc[r,'p90_peak_psu']],[yi-0.15]*2,color=CR[r],lw=2,alpha=0.5); a.plot(RE.loc[r,'median_peak_psu'],yi-0.15,'s',color=CR[r],ms=6,alpha=0.6)
a.set_yticks(y); a.set_yticklabels(REG); a.set_xlabel('psu\n● peak above local median   ■ peak salinity   (line: median to 90th pct)',fontsize=7.5); a.grid(axis='x',color=GRID,lw=0.6); clean(a); a.set_title('c  Event magnitude by region',loc='left',fontsize=10)
a=fig.add_subplot(gs[1,0:2])
for r in REG: a.plot(yrs,ED.loc[r].values,color=CR[r],lw=2,marker='o',ms=4,label=r)
a.set_xlim(2012.6,2025.4); a.set_xticks(yrs); a.set_ylabel('Event days per pixel per year'); a.grid(axis='y',color=GRID,lw=0.6); clean(a); a.legend(frameon=False,fontsize=8,ncol=6,loc='upper center',bbox_to_anchor=(0.5,-0.07))
ax2=a.inset_axes([0.47,0.52,0.3,0.44]); ax2.plot(V.index,V.obs,color=INK,lw=1.5,marker='o',ms=3,label='gauges (observed)'); ax2.plot(V.index,V['mod'],color='#2a78d6',lw=1.5,marker='s',ms=3,label='30 m model at gauges')
ax2.set_title(f'Check: event days per station-year, {pd.read_csv(f"{D}/event_days_gauge_vs_model_30m.csv").station.nunique()} gauges',fontsize=7.5,loc='left'); ax2.tick_params(labelsize=7); ax2.legend(fontsize=7,frameon=False); clean(ax2)
a.set_title('d  Frequency: days in high-salinity events (≥ 7 d above the local 95th percentile)',loc='left',fontsize=10)
a=fig.add_subplot(gs[1,2]); a.imshow(HY.values,cmap='Reds',vmin=0,vmax=60,aspect='auto')
for i in range(HY.shape[0]):
    for j in range(HY.shape[1]): v=HY.values[i,j]; a.text(j,i,f'{v:.0f}' if v>=0.5 else '·',ha='center',va='center',fontsize=7,color='white' if v>35 else INK)
a.set_yticks(range(len(REG))); a.set_yticklabels(REG,fontsize=8); a.set_xticks(range(len(yrs))); a.set_xticklabels([str(y)[2:] for y in yrs],fontsize=8); [a.spines[s].set_visible(False) for s in a.spines]; a.tick_params(length=0)
a.set_title("e  Days above 40 psu per pixel ('13–'25)",loc='left',fontsize=10)
fig.savefig(f'{D}/fig6_high_salinity_events_30m.png',dpi=130,bbox_inches='tight'); plt.close(fig)
# ================= Figure 7: 250 m vs 30 m (Shark River estuary / Whitewater Bay), May 2020 monthly mean
m=12*7+4; MED=np.load('daily_maps/median_int16.npy',mmap_mode='r'); ed=pd.date_range('2013-01-01','2025-12-31')
dd=np.flatnonzero((ed>='2020-05-01')&(ed<='2020-05-31')); g=np.asarray(MED[dd[0]:dd[-1]+1],dtype='float32')/10; g[g<0]=np.nan; g250=np.nanmean(g,0)
m30=np.asarray(MO[m],dtype='float32')/10; m30[m30<0]=np.nan
ext250=[b.left,b.left+521*8*0.00027092271166751034,b.top-480*8*0.0002709227116675095,b.top]
box=(-81.25,-80.85,25.25,25.50)
fig,ax=plt.subplots(1,3,figsize=(16,5.2))
for a,(arr,e_,t) in zip(ax,[(g250,ext250,'a  250 m map (May 2020 mean)'),(m30,ext,'b  30 m map (May 2020 mean)'),(m30-np.nan_to_num(0),ext,'')]): pass
im=ax[0].imshow(g250,extent=ext250,cmap='Blues',vmin=0,vmax=40,interpolation='nearest'); ax[0].set_title('a  250 m (May 2020 monthly mean)',loc='left',fontsize=10)
ax[1].imshow(m30,extent=ext,cmap='Blues',vmin=0,vmax=40,interpolation='nearest'); ax[1].set_title('b  30 m (May 2020 monthly mean)',loc='left',fontsize=10)
C=np.load('map30/C_int16.npy',mmap_mode='r'); cdet=np.asarray(C[m],dtype='float32'); cdet[cdet==-32768]=np.nan; cdet/=10
im2=ax[2].imshow(cdet,extent=ext,cmap=div,norm=TwoSlopeNorm(0,-5,5),interpolation='nearest'); ax[2].set_title('c  30 m detail (30 m model − 250 m map), 15 May 2020',loc='left',fontsize=10)
for a in ax: a.set_xlim(box[0],box[1]); a.set_ylim(box[2],box[3]); a.set_xticks([]); a.set_yticks([])
fig.colorbar(im,ax=ax[:2],shrink=0.8,pad=0.02).set_label('Salinity (psu)'); fig.colorbar(im2,ax=ax[2],shrink=0.8,pad=0.02).set_label('psu')
fig.savefig(f'{D}/fig7_250m_vs_30m.png',dpi=130,bbox_inches='tight'); print('figs done')
