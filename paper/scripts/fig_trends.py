"""Trends 2000-2025 (final model, surface-water mask): (a) seasonal Theil-Sen slope of monthly anomalies, (b) cells significant after FDR,
(c) slope distribution by ecosystem, (d) gauge vs model trends at long-record gauges (held out in space; maps at gauge cells)."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
INK='#3d3d3a'; MUTED='#77756e'; GRID='#e6e5df'; BLUE='#2a78d6'; ORANGE='#eb6834'; AQUA='#1baf7a'
plt.rcParams.update({'font.size':7.5,'axes.edgecolor':GRID,'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.titlelocation':'left','axes.titlesize':8})
U='/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model/best_two_step'
z=np.load(f'{U}/long_masked_final/results/pixel_long_masked_2000_2025.npz'); IDX=np.load('work/map/in/px30_index.npy'); ECO=np.load(f'{U}/events_masked/in/eco30.npy')
H,W=3839,4171; st=4
def grid(v):
    g=np.full(H*W,np.nan,'float32'); g[IDX]=v; return g.reshape(H,W)[::st,::st]
sl=z['slope_monthly']; q=z['q_monthly']
fig=plt.figure(figsize=(7.09,6.6)); gs=fig.add_gridspec(2,2,hspace=0.32,wspace=0.25)
a=fig.add_subplot(gs[0,0]); im=a.imshow(grid(sl),cmap='RdBu_r',norm=TwoSlopeNorm(0,-0.15,0.15),interpolation='nearest'); a.set_xticks([]); a.set_yticks([])
cb=fig.colorbar(im,ax=a,fraction=0.04,pad=0.01,extend='both'); cb.set_label('psu yr⁻¹',fontsize=7); a.set_title('a  Trend in monthly salinity anomalies, 2000–2025')
sig=np.where(np.isfinite(q),np.where((q<0.05)&(sl>0),2,np.where((q<0.05)&(sl<0),0,1)),np.nan)
a=fig.add_subplot(gs[0,1]); a.imshow(grid(sig),cmap=ListedColormap([BLUE,'#e6e5df',ORANGE]),vmin=0,vmax=2,interpolation='nearest'); a.set_xticks([]); a.set_yticks([])
from matplotlib.patches import Patch
a.legend(handles=[Patch(color=ORANGE,label=f'significant increase ({100*np.nanmean(sig==2):.1f} % of cells)'),Patch(color='#e6e5df',label='not significant')],loc='lower left',fontsize=6.3,frameon=False)
a.set_title('b  Significant after false-discovery correction (q < 0.05)')
for ax_ in fig.axes: [ax_.spines[s].set_visible(False) for s in ax_.spines]
a=fig.add_subplot(gs[1,0]); EN=['Freshwater marsh & prairie','Swamp forest & scrub','Salt marsh','Mangrove','Inland & estuarine water','Open bay']
data=[sl[(ECO==i)&np.isfinite(sl)] for i in range(6)]
bp=a.boxplot(data,vert=False,whis=(5,95),showfliers=False,widths=0.55,patch_artist=True,medianprops=dict(color=INK,lw=1.2))
for b_ in bp['boxes']: b_.set(facecolor='#cfe0f5',edgecolor=BLUE)
a.axvline(0,color=MUTED,lw=0.8); a.set_yticks(range(1,7)); a.set_yticklabels([e.replace(' & prairie','').replace(' & scrub','').replace('Inland & estuarine water','Inland/estuarine') for e in EN],fontsize=6.8)
a.set_xlabel('trend (psu yr⁻¹); box 25–75 %, whiskers 5–95 %'); a.grid(axis='x',color=GRID,lw=0.5); a.invert_yaxis(); [a.spines[s].set_visible(False) for s in ('top','right')]
a.set_title('c  Trends by ecosystem')
a=fig.add_subplot(gs[1,1]); T=pd.read_csv('work/final/gauge/gauge_trends_2000_2025.csv'); x=T.obs_slope
lim=0.36; a.axhline(0,color=GRID,lw=0.8); a.axvline(0,color=GRID,lw=0.8); a.plot([-lim,lim],[-lim,lim],color=MUTED,lw=0.8,ls='--')
t=T.dropna(subset=['model_salinity_30m_onestep_masked_slope'])
a.scatter(t.obs_slope,t.model_salinity_30m_onestep_masked_slope,s=12,color=BLUE,alpha=0.75,lw=0,label=f'30 m maps at gauge cell (ρ = {t.obs_slope.corr(t.model_salinity_30m_onestep_masked_slope,method="spearman"):.2f})')
a.scatter(x,T.pred_cv_slope,s=14,facecolor='none',edgecolor=AQUA,lw=0.9,label=f'model, station held out (ρ = {x.corr(T.pred_cv_slope,method="spearman"):.2f})')
a.set_xlim(-0.15,lim); a.set_ylim(-0.15,lim); a.set_aspect('equal'); a.legend(frameon=False,fontsize=6,loc='upper left',handletextpad=0.2)
a.set_xlabel('gauge trend (psu yr⁻¹)'); a.set_ylabel('model trend (psu yr⁻¹)'); [a.spines[s].set_visible(False) for s in ('top','right')]
a.set_title(f'd  Trends at {len(T)} long-record gauges')
fig.savefig('Fig08_trends.png',dpi=300,bbox_inches='tight',facecolor='white'); print('ok', [round(float(np.median(d)),4) for d in data])
