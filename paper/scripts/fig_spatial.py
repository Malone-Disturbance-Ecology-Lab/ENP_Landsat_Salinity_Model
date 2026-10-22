"""Spatial patterns, final model with surface-water mask, 2000-2025: climatology maps and Venice-system areas."""
import numpy as np, pandas as pd, rasterio, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
INK='#3d3d3a'; MUTED='#77756e'; GRID='#e6e5df'
plt.rcParams.update({'font.size':7.5,'axes.edgecolor':GRID,'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.titlelocation':'left','axes.titlesize':8})
U='/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model/best_two_step'
with rasterio.open(f'{U}/climatology_final/ENP_salinity30m_climatology_2000_2025.tif') as s: B=s.read().astype('float32'); B[B<0]=np.nan; B/=10
mean,dry,wet,may,octo,sd,pw,p95=B
IDX=np.load('work/map/in/px30_index.npy'); ECO=np.load(f'{U}/events_masked/in/eco30.npy'); REG=np.load(f'{U}/events_masked/in/region30.npy')
A=0.000818
# ---- statistics
cls=[0,0.5,5,18,30,40,99]; lab=['fresh (< 0.5)','oligohaline (0.5–5)','mesohaline (5–18)','polyhaline (18–30)','euhaline (30–40)','hypersaline (> 40)']
rows=[]
for nm,a in [('mean',mean),('dry season Mar–May',dry),('wet season Aug–Oct',wet),('May',may),('October',octo)]:
    v=a.ravel()[IDX]; ok=np.isfinite(v); h=np.histogram(v[ok],cls)[0]
    rows.append(dict(map=nm,cells=int(ok.sum()),**{l:round(100*c/ok.sum(),1) for l,c in zip(lab,h)}))
S=pd.DataFrame(rows); S.to_csv('area_by_salinity_class.csv',index=False); print(S.to_string())
EN=['Freshwater marsh & prairie','Swamp forest & scrub','Salt marsh','Mangrove','Inland & estuarine water','Open bay']
RN=['Big Cypress','Shark River Slough','Taylor Slough','C-111','Gulf coastal waters','Florida Bay']
G=[]
for g,cat,names in [('ecosystem',ECO,EN),('region',REG,RN)]:
    for i,n in enumerate(names):
        m=cat==i; f=lambda a: float(np.nanmedian(a.ravel()[IDX][m]))
        G.append(dict(group=g,name=n,km2=round(m.sum()*A),mean=f(mean),dry=f(dry),wet=f(wet),may=f(may),october=f(octo),interannual_sd=f(sd),pct_days_wet=f(pw),p95_monthly=f(p95)))
G=pd.DataFrame(G).round(2); G.to_csv('climatology_by_group.csv',index=False); print(G.to_string())
# ---- figure
st=3
def show(a,arr,title,cmap,vmin,vmax,label,norm=None):
    im=a.imshow(arr[::st,::st],cmap=cmap,vmin=None if norm else vmin,vmax=None if norm else vmax,norm=norm,interpolation='nearest'); a.set_xticks([]); a.set_yticks([]); [a.spines[k].set_visible(False) for k in a.spines]
    cb=fig.colorbar(im,ax=a,fraction=0.04,pad=0.01,format='%g'); cb.set_label(label,fontsize=6.5); cb.ax.tick_params(labelsize=6); a.set_title(title)
fig,ax=plt.subplots(2,3,figsize=(7.09,5.0)); plt.subplots_adjust(wspace=0.08,hspace=0.18)
vb=[0,0.5,5,10,18,25,30,35,40,45]; cm=plt.get_cmap('viridis',len(vb)-1); nrm=BoundaryNorm(vb,cm.N)
show(ax[0,0],mean,'a  Mean, 2000–2025',cm,0,45,'psu',nrm); show(ax[0,1],dry,'b  Dry season (Mar–May)',cm,0,45,'psu',nrm); show(ax[0,2],wet,'c  Wet season (Aug–Oct)',cm,0,45,'psu',nrm)
show(ax[1,0],dry-wet,'d  Dry minus wet season','Reds',0,20,'psu'); show(ax[1,1],sd,'e  Interannual SD of annual means','magma_r',0,6,'psu'); show(ax[1,2],pw,'f  Days with surface water','Blues',50,100,'%')
fig.savefig('Fig06_spatial_patterns.png',dpi=300,bbox_inches='tight',facecolor='white')
# 30 m detail inset figure data: Shark River mouth window
r0,r1,c0,c1=1500,1900,700,1200; np.save('detail_window.npy',mean[r0:r1,c0:c1])
print('domain cells with any wet day',int(np.isfinite(mean.ravel()[IDX]).sum()),'of',len(IDX))
