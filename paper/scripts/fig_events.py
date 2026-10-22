"""High-salinity events 2000-2025 (final model, surface-water mask): (a) events per year per cell, (b) event days per pixel by region and year,
(c) annual event days at gauges vs the maps at the gauge cells."""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
INK='#3d3d3a'; MUTED='#77756e'; GRID='#e6e5df'; BLUE='#2a78d6'; ORANGE='#eb6834'
plt.rcParams.update({'font.size':7.5,'axes.edgecolor':GRID,'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.titlelocation':'left','axes.titlesize':8})
U=f'{ROOT}/best_two_step/long_masked_final/results'
z=np.load(f'{U}/pixel_long_masked_2000_2025.npz'); IDX=np.load(f'{WORK}/map/in/px30_index.npy'); H,W=3839,4171
fig=plt.figure(figsize=(7.09,5.6)); gs=fig.add_gridspec(2,1,height_ratios=[1.15,1],hspace=0.35)
a=fig.add_subplot(gs[1,0]); R=pd.read_csv(f'{WORK}/final/gauge/gauge_event_days_2000_2025.csv'); m=R.groupby('year')[['obs','model']].mean()
a.plot(m.index,m.obs,color=ORANGE,lw=1.8,marker='o',ms=2.5,label='gauges'); a.plot(m.index,m.model,color=BLUE,lw=1.8,marker='o',ms=2.5,label='maps at gauge cells')
a.set_ylabel('event days per station-year'); a.grid(axis='y',color=GRID,lw=0.6); [a.spines[s].set_visible(False) for s in ('top','right')]
a.legend(frameon=False,fontsize=6.5,loc='upper right'); a.set_xlim(2000.5,2025.5)
a.set_title(f'b  Event days at {R.station.nunique()} gauges (station-year r = {np.corrcoef(R.obs,R.model)[0,1]:.2f})')
a=fig.add_subplot(gs[0,0]); y=pd.read_csv(f'{U}/high_salinity_events_by_group_year_2000_2025.csv'); y=y[y.group=='region']
REG=['Big Cypress','Shark River Slough','Taylor Slough','C-111','Gulf coastal waters','Florida Bay']
P=y.pivot(index='name',columns='year',values='event_days_per_pixel').loc[REG]; im=a.imshow(P.values,cmap='Reds',vmin=0,vmax=110,aspect='auto')
a.set_yticks(range(len(REG))); a.set_yticklabels(REG); a.set_xticks(range(P.shape[1])); a.set_xticklabels([f"'{str(c)[2:]}" for c in P.columns])
for i in range(P.shape[0]):
    for j in range(P.shape[1]):
        v=P.values[i,j]; a.text(j,i,f'{v:.0f}',ha='center',va='center',fontsize=6,color='white' if v>60 else INK)
a.tick_params(length=0); [a.spines[s].set_visible(False) for s in a.spines]
cb=fig.colorbar(im,ax=a,fraction=0.02,pad=0.01); cb.set_label('event days per pixel',fontsize=7); a.set_title('a  Days in high-salinity events by region and year (2000 = April–December)')
fig.savefig('Fig09_events.png',dpi=300,bbox_inches='tight',facecolor='white')
