"""ESSD Fig. 8: the daily surface-water mask. (a,b) codes on 15 May and 15 Oct 2015; (c) share of April days with surface water;
(d) mask vs EDEN gauge dry flags."""
import numpy as np, pandas as pd, rasterio, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
INK='#3d3d3a'; MUTED='#77756e'; GRID='#e6e5df'
plt.rcParams.update({'font.size':7.5,'axes.edgecolor':GRID,'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.titlelocation':'left','axes.titlesize':8,'font.family':'DejaVu Sans'})
cols={0:'#c8553d',1:'#2f6fb5',2:'#f0c75e',3:'#b9d7ea',4:'#5aa469',5:'#d9d9d9',9:'#6b4f3a',255:'#ffffff'}
keys=[0,1,2,3,4,5,9,255]; cmap=ListedColormap([cols[k] for k in keys]); norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5,4.5,6,10,256],cmap.N)
lab={0:'dry',1:'wet',2:'uncertain',3:'permanent water',4:'tidal coastal wetland',5:'not assessed',9:'never wet'}
W='work/mask/troot/wetdry30m/2015'
fig=plt.figure(figsize=(7.09,6.3)); gs=fig.add_gridspec(2,2,height_ratios=[1.25,1],hspace=0.42,wspace=0.42)
for j,(f,t) in enumerate([(f'{W}/ENP_wetdry30m_daily_2015_05.tif','a  15 May 2015 (drought, end of dry season)'),(f'{W}/ENP_wetdry30m_daily_2015_10.tif','b  15 October 2015 (wet season)')]):
    ax=fig.add_subplot(gs[0,j]); A=rasterio.open(f).read(15)[::3,::3]; ax.imshow(A,cmap=cmap,norm=norm,interpolation='nearest'); ax.set_title(t); ax.axis('off')
fig.legend(handles=[Patch(color=cols[k],label=lab[k]) for k in [0,2,1,3,4,5,9]],loc='upper center',bbox_to_anchor=(0.5,0.50),ncol=7,frameon=False,fontsize=7,handlelength=1.2,columnspacing=1.0)
D=pd.read_csv('work/trend/change30_masked/monthly_pct_days_with_surface_water.csv'); D['yr']=D.month.str[:4].astype(int); D['mon']=D.month.str[-2:].astype(int)
ax=fig.add_subplot(gs[1,0]); A4=D[(D.mon==4)]
for name,c in [('Freshwater marsh & prairie','#1f5fae'),('Swamp forest & scrub','#56b1e3')]:
    s=A4[(A4.group=='ecosystem')&(A4.name==name)].sort_values('yr'); ax.plot(s.yr,s.pct_days_wet,color=c,lw=2,marker='o',ms=4,mec='white',mew=0.8)
    ax.lines[-1].set_label(name)
ax.set_ylim(0,100); ax.set_xlim(2012.6,2025.4); ax.legend(loc='lower right',frameon=False,fontsize=6.8); ax.set_xticks(range(2013,2026,2)); ax.set_ylabel('April days with surface water (%)'); ax.grid(axis='y',color=GRID,lw=0.6)
for s in ('top','right'): ax.spines[s].set_visible(False)
ax.set_title('c  End-of-dry-season wetness (April), 2013–2025')
V=pd.read_csv('work/mask/out/wetdry_gauge_validation_summary.csv',index_col=0)
ax=fig.add_subplot(gs[1,1]); rows=[('gauge dry','gauge dry'),('gauge wet','gauge wet')]
for i,(lb,colname) in enumerate(rows):
    left=0
    for k,code in [('dry',0),('uncertain',2),('wet',1)]:
        v=100*V.loc[k,colname]; ax.barh(i,v,left=left,color=cols[code],height=0.55,edgecolor='white',linewidth=2)
        if v>6: ax.text(left+v/2,i,f'{v:.0f}%',ha='center',va='center',color='white' if code!=2 else INK,fontsize=7)
        left+=v
ax.set_yticks([0,1]); ax.set_yticklabels([r[0] for r in rows]); ax.set_xlim(0,100); ax.set_xlabel('share of gauge-days called dry / uncertain / wet (%)')
for s in ('top','right','left'): ax.spines[s].set_visible(False)
ax.tick_params(axis='y',length=0); ax.set_title('d  Mask vs EDEN gauge dry flags (42 gauges, 2013–2025)')
fig.savefig('Fig8_surface_water_mask.png',dpi=300,bbox_inches='tight',facecolor='white'); fig.savefig('Fig8_surface_water_mask.pdf',bbox_inches='tight',facecolor='white'); print('ok')
