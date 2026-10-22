"""Drivers of modelled salinity (final model): (a) held-out RMSE rise when a predictor group is shuffled, by salinity zone;
(b) group response curves (each station's own percentiles) in the 5-25 psu zone; (c) dominant driver group of day-to-day
variation at 500 m nodes; (d) share of day-to-day variance carried by the anchor stations vs hydroclimate."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
INK='#3d3d3a'; MUTED='#77756e'; GRID='#e6e5df'; CAT=['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300']
plt.rcParams.update({'font.size':7.5,'axes.edgecolor':GRID,'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.titlelocation':'left','axes.titlesize':8})
COL={'anchor stations':CAT[1],'water level':CAT[0],'climate':CAT[2],'sea level':CAT[4],'flows':CAT[3],'season':'#9a9890','site variables':'#5d5b55'}
fig=plt.figure(figsize=(7.09,6.8)); gs=fig.add_gridspec(2,2,hspace=0.38,wspace=0.3)
# (a)
P=pd.concat([pd.read_csv('permutation_importance_groups.csv'),pd.read_csv('permutation_importance_combined.csv')])
P['group']=P.group.replace({'hydroclimate (water level, flows, climate, sea level, season)':'hydroclimate combined','all daily drivers':'all daily drivers combined'})
order=['site variables','anchor stations','climate','water level','sea level','flows','season','hydroclimate combined','all daily drivers combined']
Z=['< 5 psu','5–25 psu','> 25 psu']; a=fig.add_subplot(gs[0,0]); y=np.arange(len(order))
for k,z in enumerate(Z):
    v=P[P.zone==z].set_index('group').loc[order,'rmse_rise']; a.barh(y+(k-1)*0.26,v,height=0.24,color=['#b9d3f2','#2a78d6','#123f73'][k],label=f'stations {z}')
a.set_yticks(y); a.set_yticklabels(order,fontsize=6.6); a.invert_yaxis(); a.axvline(0,color=MUTED,lw=0.6); a.axhline(6.5,color=GRID,lw=0.8)
a.set_xlabel('rise in held-out RMSE when shuffled (psu)'); a.legend(frameon=False,fontsize=6.2,loc='center right',bbox_to_anchor=(1.0,0.42)); a.grid(axis='x',color=GRID,lw=0.5)
a.set_title('a  Importance of predictor groups'); [a.spines[s].set_visible(False) for s in ('top','right')]
# (b)
R=pd.read_csv('group_response.csv'); a=fig.add_subplot(gs[0,1]); r=R[R.zone=='5–25 psu']
for g in ['anchor stations','climate','water level','sea level','flows']:
    s=r[r.group==g]; a.plot(s.q,s.salinity,color=COL[g],lw=2); a.text(96,s.salinity.iloc[-1]+{'water level':0.5,'climate':-0.6}.get(g,0),g,fontsize=6.4,va='center',color=INK)
a.set_xlim(3,124); a.set_xticks([5,25,50,75,95]); a.set_xlabel("driver group at each station's own percentile"); a.set_ylabel('mean predicted salinity (psu)')
a.grid(axis='y',color=GRID,lw=0.5); a.set_title('b  Response to each driver group, stations 5–25 psu'); [a.spines[s].set_visible(False) for s in ('top','right')]
# (c, d)
z=np.load('node_attribution.npz'); rc=np.load('node_rc.npy'); H2,W2=480,521
grp=['water level','flows','climate','anchor stations','sea level','season']; D=np.stack([z[f'var_drop_{g}'] for g in grp]); V=z['var_full']
dom=np.argmax(D,0).astype(float); dom[V<0.25]=np.nan
def grid(v):
    g=np.full((H2,W2),np.nan,'float32'); g[rc[0],rc[1]]=v; return g[::2,::2]
a=fig.add_subplot(gs[1,0]); cm=ListedColormap([COL[g] for g in grp])
a.imshow(grid(dom),cmap=cm,vmin=-0.5,vmax=5.5,interpolation='nearest'); a.set_xticks([]); a.set_yticks([]); [a.spines[s].set_visible(False) for s in a.spines]
from matplotlib.patches import Patch
okn=np.isfinite(dom); sh={g:100*np.mean(dom[okn]==i) for i,g in enumerate(grp)}
a.legend(handles=[Patch(color=COL[g],label=f'{g} ({sh[g]:.0f} %)') for g in grp if sh[g]>=0.5],loc='lower left',fontsize=6.2,frameon=False)
a.set_title('c  Driver group with the largest share of\n    day-to-day variation (SD ≥ 0.5 psu)')
a=fig.add_subplot(gs[1,1]); hyd=D[[0,1,2,4]].sum(0); anc=D[3]; frac=np.where(V>=0.25,anc/np.maximum(anc+hyd,1e-6),np.nan)
im=a.imshow(grid(np.clip(frac,0,1)),cmap='RdBu_r',vmin=0,vmax=1,interpolation='nearest'); a.set_xticks([]); a.set_yticks([]); [a.spines[s].set_visible(False) for s in a.spines]
cb=fig.colorbar(im,ax=a,fraction=0.045,pad=0.01); cb.set_label('anchor-station share',fontsize=6.5); cb.ax.tick_params(labelsize=6)
a.set_title('d  Anchor stations vs hydroclimate\n    (share of attributed variance)')
fig.savefig('Fig10_drivers.png',dpi=300,bbox_inches='tight',facecolor='white'); print('ok')
