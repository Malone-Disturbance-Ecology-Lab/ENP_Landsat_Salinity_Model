"""Temporal patterns, final model with surface-water mask, 2000-2025: (a) seasonal cycle by ecosystem, (b) annual mean anomaly by region and year,
(c) May salinity vs December-May rainfall by region, (d) share of freshwater-marsh cells with surface water in April."""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
INK='#3d3d3a'; MUTED='#77756e'; GRID='#e6e5df'; CAT=['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4','#008300']
plt.rcParams.update({'font.size':7.5,'axes.edgecolor':GRID,'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.titlelocation':'left','axes.titlesize':8})
U=f'{ROOT}/best_two_step/long_masked_final/results'
M=pd.read_csv(f'{U}/monthly_series_by_group_2000_2025.csv',parse_dates=['month']); M['yr']=M.month.dt.year; M['mo']=M.month.dt.month
DR=pd.read_csv('annual_drivers.csv',index_col=0)
ECO=['Freshwater marsh & prairie','Swamp forest & scrub','Salt marsh','Mangrove','Inland & estuarine water','Open bay']
REG=['Shark River Slough','Taylor Slough','C-111','Big Cypress','Gulf coastal waters','Florida Bay']
fig=plt.figure(figsize=(7.09,7.4)); gs=fig.add_gridspec(2,2,hspace=0.42,wspace=0.55,height_ratios=[1,1])
a=fig.add_subplot(gs[0,0]); E=M[M.group=='ecosystem']; cl=E.groupby(['name','mo']).mean_salinity_wet.mean().unstack()
for i,e in enumerate(ECO):
    a.plot(range(1,13),cl.loc[e],color=CAT[i],lw=2,marker='o',ms=3); off={'Freshwater marsh & prairie':-1.6,'Swamp forest & scrub':1.2}.get(e,0)
    a.text(12.3,cl.loc[e,12]+off,{'Freshwater marsh & prairie':'Freshwater marsh','Swamp forest & scrub':'Swamp forest','Inland & estuarine water':'Inland/estuarine'}.get(e,e),fontsize=6.3,va='center',color=INK)
a.set_xticks(range(1,13)); a.set_xticklabels(list('JFMAMJJASOND')); a.set_xlim(0.7,15.6); a.set_ylabel('mean salinity over wet days (psu)'); a.grid(axis='y',color=GRID,lw=0.5)
a.set_title('a  Seasonal cycle by ecosystem, 2000–2025')
a=fig.add_subplot(gs[0,1]); R=M[(M.group=='region')&(M.yr>=2001)]; an=R.groupby(['name','yr']).mean_salinity_wet.mean().unstack().loc[REG]; an=an.sub(an.mean(1),axis=0)
im=a.imshow(an.values,cmap='RdBu_r',norm=TwoSlopeNorm(0,-8,8),aspect='auto'); a.set_yticks(range(len(REG))); a.set_yticklabels([r.replace('Gulf coastal waters','Gulf coast') for r in REG],fontsize=6.5)
a.set_xticks(range(0,25,4)); a.set_xticklabels([str(2001+k) for k in range(0,25,4)]); a.tick_params(length=0); [a.spines[s].set_visible(False) for s in a.spines]
cb=fig.colorbar(im,ax=a,fraction=0.05,pad=0.03); cb.set_label('annual anomaly (psu)',fontsize=7)
a.set_title('b  Annual mean anomaly by region, 2001–2025')
a=fig.add_subplot(gs[1,0]); my=M[(M.group=='region')&(M.mo==5)].pivot(index='yr',columns='name',values='mean_salinity_wet')
for i,rg in enumerate(['Taylor Slough','C-111','Florida Bay']):
    y=my[rg]; x=DR.loc[y.index,'rain_DecMay_mm']; a.scatter(x,y,s=16,color=CAT[i],lw=0,alpha=0.85); b=np.polyfit(x,y,1); xx=np.array([x.min(),x.max()]); a.plot(xx,np.polyval(b,xx),color=CAT[i],lw=1.2)
    a.text(xx[1]+8,np.polyval(b,xx[1]),f'{rg}\n(r = {np.corrcoef(x,y)[0,1]:.2f})',fontsize=6.3,va='center')
a.set_xlim(190,760); a.set_xlabel('December–May rainfall, park mean (mm)'); a.set_ylabel('May mean salinity (psu)'); a.grid(axis='y',color=GRID,lw=0.5)
a.set_title('c  End-of-dry-season salinity vs dry-season rain')
a=fig.add_subplot(gs[1,1]); fw=E[(E.name=='Freshwater marsh & prairie')&(E.mo==4)].set_index('yr').pct_cells_with_surface_water
a.bar(fw.index,fw.values,color=CAT[0],width=0.7); a.set_ylim(0,105); a.set_ylabel('cells with surface water in April (%)'); a.grid(axis='y',color=GRID,lw=0.5)
for yr in fw.nsmallest(4).index: a.text(yr,fw[yr]+2,str(yr),ha='center',fontsize=6,color=INK)
a.set_title('d  Freshwater marsh with surface water, April')
for ax_ in fig.axes[:1]+fig.axes[3:]: [ax_.spines[s].set_visible(False) for s in ('top','right')]
fig.savefig('Fig07_temporal_patterns.png',dpi=300,bbox_inches='tight',facecolor='white'); print('ok')
