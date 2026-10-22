"""Daily held-out model (station's spatial block left out) and gauge salinity at six stations, 2000-2025, with the 90 % interval."""
import json, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
INK='#3d3d3a'; MUTED='#77756e'; GRID='#e6e5df'; BLUE='#2a78d6'; ORANGE='#eb6834'
plt.rcParams.update({'font.size':7.5,'axes.edgecolor':GRID,'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.titlelocation':'left','axes.titlesize':8})
C=pd.read_parquet('work/final/cv_predictions_final.parquet'); Q=json.load(open('work/final/q_levels.json')); bins=np.array(Q['bins']); q=np.array(Q['q'])
ST=pd.read_csv('station_skill_final.csv').set_index('station')
PICK=[('22908295','freshwater marsh, Shark River Slough'),('ENPNR','Taylor Slough'),('ENPGI','Shark River Slough estuary'),('ENPBS','C-111 basin'),('ENPBA','Florida Bay'),('ENPCA','Big Cypress coast')]
fig,ax=plt.subplots(len(PICK),1,figsize=(7.09,9.2),sharex=True); plt.subplots_adjust(hspace=0.55)
for i,(s,lab) in enumerate(PICK):
    a=ax[i]; g=C[C.station==s].set_index('date').sort_index(); p=g.pred_cv; hw=q[np.clip(np.digitize(p,bins[1:-1]),0,3)]
    a.fill_between(p.index,np.clip(p-hw,0,None),p+hw,color=BLUE,alpha=0.12,lw=0,label='90 % interval')
    a.plot(p.index,p,color=BLUE,lw=0.6,label='model, station held out'); o=g.salinity.dropna()
    a.scatter(o.index,o.values,s=0.6,color=ORANGE,lw=0,label='gauge (daily mean)',rasterized=True)
    r=ST.loc[s]; a.set_title(f'{"abcdef"[i]}  {s} · {lab}    r {r.r:.2f} · KGE {r.KGE:.2f} · RMSE {r.RMSE:.1f} · bias {r.bias:+.1f} psu',fontsize=7.5)
    a.set_ylabel('psu'); a.grid(axis='y',color=GRID,lw=0.5); a.set_xlim(pd.Timestamp('2000-03-01'),pd.Timestamp('2026-01-31')); a.set_ylim(bottom=0)
    [a.spines[k].set_visible(False) for k in ('top','right')]
    if i==0: a.legend(loc='upper right',fontsize=6.5,markerscale=8,ncol=3,frameon=False)
fig.savefig('Fig04_station_series.png',dpi=300,bbox_inches='tight',facecolor='white'); print('ok')
