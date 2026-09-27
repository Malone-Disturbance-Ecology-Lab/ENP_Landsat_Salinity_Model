"""Export daily maps: monthly GeoTIFFs of daily medians, monthly mean + 90% bounds per year, station daily series, previews."""
import numpy as np, pandas as pd, rasterio, os, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT='SERVER/best_two_step/daily_maps'; os.makedirs(OUT+'/daily',exist_ok=True); os.makedirs(OUT+'/monthly',exist_ok=True)
M=np.load('daily_maps/median_int16.npy',mmap_mode='r'); qlev=np.load('daily_maps/interval_q_by_level.npy'); bins=[-1,5,15,30,1e9]
dates=pd.date_range('2013-01-01','2025-12-31'); assert int(open('daily_maps/progress.txt').read())==len(dates)
wfr=rasterio.open('SERVER/site_variables/ENP_WaterFrequency_L30.tif'); K=8
H,W=M.shape[1:]; tr=wfr.transform*rasterio.Affine.scale(K)
base=dict(driver='GTiff',width=W,height=H,crs=wfr.crs,transform=tr,compress='deflate',tiled=True,blockxsize=256,blockysize=256)
monthly={k:[] for k in ('mean','lower90','upper90')}; mkeys=[]
for (y,m),idx in pd.Series(np.arange(len(dates)),index=dates).groupby([dates.year,dates.month]):
    idx=idx.values; a=np.asarray(M[idx[0]:idx[-1]+1]); os.makedirs(f'{OUT}/daily/{y}',exist_ok=True)
    with rasterio.open(f'{OUT}/daily/{y}/ENP_salinity_daily_{y}_{m:02d}.tif','w',count=len(idx),dtype='int16',nodata=-1,predictor=2,**base) as o:
        o.write(a); [o.set_band_description(i+1,str(dates[j].date())) for i,j in enumerate(idx)]
        o.update_tags(units='salinity psu x 10 (divide by 10)',model='0.7*two-step (vegetation-spectra mean + daily anomaly) + 0.3*v3; 250 m')
    p=np.where(a>=0,a/10.0,np.nan).astype('float32'); q=qlev[np.clip(np.digitize(np.nan_to_num(p),bins[1:-1]),0,3)]
    monthly['mean'].append(np.nanmean(p,0)); monthly['lower90'].append(np.nanmean(np.clip(p-q,0,None),0)); monthly['upper90'].append(np.nanmean(p+q,0)); mkeys.append(f'{y}-{m:02d}')
    if m==12: print('exported',y,flush=True)
for k,L in monthly.items():
    A=np.stack(L)
    for y in range(2013,2026):
        ii=[i for i,s in enumerate(mkeys) if s.startswith(str(y))]
        with rasterio.open(f'{OUT}/monthly/ENP_salinity_monthly_{k}_{y}.tif','w',count=12,dtype='float32',nodata=np.nan,predictor=3,**base) as o:
            o.write(np.round(A[ii],1)); [o.set_band_description(j+1,mkeys[i]) for j,i in enumerate(ii)]
# ---- continuous daily series at every station (final model; in-sample at monitored stations)
st=pd.read_csv('sitevars/stations_coast.csv')
rows=[]
for _,s in st.iterrows():
    r,c=rasterio.transform.rowcol(tr,s.longitude,s.latitude)
    if 0<=r<H and 0<=c<W:
        v=np.asarray(M[:,r,c]).astype('float32'); v[v<0]=np.nan; rows.append(pd.DataFrame({'station':s.station,'date':dates,'model_salinity':v/10}))
S=pd.concat(rows); AF=pd.read_pickle('all_features.pkl'); AF['date']=pd.to_datetime(AF.date)
S=S.merge(AF[['station','date','salinity']].rename(columns={'salinity':'gauge_salinity'}),on=['station','date'],how='left')
S.to_csv(f'{OUT}/station_daily_model_series_2013_2025.csv.gz',index=False)
# ---- previews
ext=[wfr.bounds.left,wfr.bounds.right,wfr.bounds.bottom,wfr.bounds.top]
show=['2015-05-20','2015-10-20','2024-05-20','2024-10-20','2025-05-20','2025-10-20']
fig,ax=plt.subplots(2,3,figsize=(16,10))
for a,d in zip(ax.T.ravel(),show):
    x=np.asarray(M[dates.get_loc(pd.Timestamp(d))],dtype='float32'); x[x<0]=np.nan
    im=a.imshow(x/10,extent=ext,cmap='Blues',vmin=0,vmax=45); a.set_title(d+' (daily)',loc='left'); a.set_xticks([]); a.set_yticks([])
fig.colorbar(im,ax=ax,shrink=.6,label='Salinity (psu)'); fig.savefig(f'{OUT}/preview_daily_dry_vs_wet.png',dpi=70)
A=np.stack(monthly['mean']); fig,ax=plt.subplots(1,2,figsize=(13,5))
ax[0].imshow(np.nanmean(A,0),extent=ext,cmap='Blues',vmin=0,vmax=45); ax[0].set_title('Mean salinity 2013–2025',loc='left')
im=ax[1].imshow(np.nanstd(A,0),extent=ext,cmap='Oranges',vmin=0,vmax=12); ax[1].set_title('Month-to-month SD',loc='left'); fig.colorbar(im,ax=ax[1],shrink=.8)
for a in ax: a.set_xticks([]); a.set_yticks([])
fig.savefig(f'{OUT}/preview_mean_and_variability.png',dpi=70)
print('done')
