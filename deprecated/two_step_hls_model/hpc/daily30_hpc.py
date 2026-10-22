"""Full-model DAILY salinity at 30 m (native HLS/Landsat grid) for one month, for HPC array jobs.
Usage:  python daily30_hpc.py <month_index 0..155> <input_dir> <output_dir>
Writes <output_dir>/YYYY/ENP_salinity30m_daily_YYYY_MM.tif (one band per day; int16 psu*10; -1 = no data; EPSG:4326).
Same model and inputs as the 250 m maps and the 30 m monthly maps (see README_HPC.md). ~7.6 M pixels x ~30 days per job."""
import sys, os, json, time, pickle, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, rasterio
from rasterio.transform import Affine
m_job=int(sys.argv[1]); IN=sys.argv[2]; OUTD=sys.argv[3]; t0=time.time()
H,W=3839,4171; NM=156; TILE=int(os.environ.get('TILE_ROWS','64'))
months=pd.date_range('2013-01-01',periods=NM,freq='MS'); mid=months+pd.Timedelta(days=14); edates=pd.date_range('2013-01-01','2025-12-31')
M0=months[m_job]; days=edates[(edates>=M0)&(edates<M0+pd.offsets.MonthBegin(1))]
FM=pickle.load(open(f'{IN}/final_models.pkl','rb')); F1=list(FM['m1'].feature_names_in_); F2=list(FM['m2'].feature_names_in_); F3=list(FM['m3'].feature_names_in_)
try:        # compiled models (treelite/tl2cgen) are ~12x faster and give identical predictions
    import tl2cgen
    P={k:tl2cgen.Predictor(f'{IN}/{k}.so',nthread=int(os.environ.get('NTHREAD','4'))) for k in ('m1','m2','m3')}
    def pred(k,cols,F): return P[k].predict(tl2cgen.DMatrix(np.column_stack([cols[c] for c in F]).astype('float32'))).ravel()
except Exception as e:
    print('tl2cgen unavailable, using scikit-learn predict:',e,flush=True)
    def pred(k,cols,F): return FM[k].predict(np.column_stack([cols[c] for c in F]).astype('float32'))
idx=np.load(f'{IN}/px30_index.npy'); STK=json.load(open(f'{IN}/static30_keys.json')); STA=np.load(f'{IN}/static30.npy',mmap_mode='r')
LL=np.load(f'{IN}/lonlat30.npy',mmap_mode='r'); UT=np.load(f'{IN}/utm30.npy',mmap_mode='r')
SPK=json.load(open(f'{IN}/spec_grid_keys.json')); SPG=np.load(f'{IN}/spec_grid.npy',mmap_mode='r')
Z=np.load(f'{IN}/gridmet_filled.npz'); glat_s=Z['glat_s']; glon=Z['glon']; gd=pd.to_datetime(Z['dates']); GA={k:Z[k] for k in Z.files if k not in ('glat_s','glon','dates')}
EDM=np.load(f'{IN}/eden_daily/crop_meta.npy'); ey0,ex0,eh,ew=EDM; eh,ew=int(eh),int(ew)
EF={k:np.load(f'{IN}/eden_daily/{k}.npy',mmap_mode='r') for k in ['cell','cell30','Rcell','5km30','5km90','R5']}
DRV=pd.read_csv(f'{IN}/daily_global_drivers.csv',parse_dates=['date']).set_index('date')   # flows (+_dep), anchors
flm=[DRV.loc[months[m]:months[m]+pd.offsets.MonthEnd(0),['SRS180','TSC111180']].mean() for m in range(NM)]
gidx=pd.Index(gd)
def bil(A2,fy,fx,iy,ix):
    v=[A2[iy,ix],A2[iy,ix+1],A2[iy+1,ix],A2[iy+1,ix+1]]; w=[(1-fy)*(1-fx),(1-fy)*fx,fy*(1-fx),fy*fx]
    num=sum(np.where(np.isfinite(a),a*b,0) for a,b in zip(v,w)); den=sum(np.where(np.isfinite(a),b,0) for a,b in zip(v,w))
    return np.where(den>0,num/np.maximum(den,1e-9),np.nan).astype('float32')
OUT=np.full((len(days),H,W),-1,'int16')
msel=list(range(max(m_job-7,0),min(m_job+8,NM)))            # months needed for the 12-month centred smoothing + interpolation
RR=[int(v) for v in os.environ.get('ROWS','0,%d'%H).split(',')]   # optional row range for testing, e.g. ROWS=1500,1564
for r0 in range(RR[0],RR[1],TILE):
    r1=min(r0+TILE,RR[1]); a,b=np.searchsorted(idx,[r0*W,r1*W]); n=b-a
    if n==0: continue
    pix=idx[a:b]; rr,cc=pix//W,pix%W; st={k:np.asarray(STA[i,a:b]) for i,k in enumerate(STK)}
    lon,lat=np.asarray(LL[0,a:b]),np.asarray(LL[1,a:b]); X,Y=np.asarray(UT[0,a:b]),np.asarray(UT[1,a:b])
    giy=np.clip(np.searchsorted(glat_s,lat)-1,0,len(glat_s)-2); gix=np.clip(np.searchsorted(glon,lon)-1,0,len(glon)-2)
    gfy=np.clip((lat-glat_s[giy])/(glat_s[giy+1]-glat_s[giy]),0,1); gfx=np.clip((lon-glon[gix])/(glon[gix+1]-glon[gix]),0,1)
    ry=(Y-ey0)/400; rx=(X-ex0)/400; eiy=np.clip(np.floor(ry).astype(int),0,eh-2); eix=np.clip(np.floor(rx).astype(int),0,ew-2); efy=np.clip(ry-eiy,0,1); efx=np.clip(rx-eix,0,1)
    vi=(rr-4)/16.0; va=max(int(np.floor(vi.min()))-1,0); vb=min(int(np.floor(vi.max()))+3,SPG.shape[2])
    viy=np.clip(np.floor(vi).astype(int),va,vb-2)-va; vfy=np.clip(vi-(viy+va),0,1); vj=(cc-4)/8.0; vix=np.clip(np.floor(vj).astype(int),0,SPG.shape[3]-2); vfx=np.clip(vj-vix,0,1)
    SPt=np.asarray(SPG[:,:,va:vb,:],dtype='float32')
    Bm={}
    for m in range(max(m_job-13,0),min(m_job+14,NM)):       # enough months for centred 12-month means of the needed months
        cols=dict(st)
        for q,k in enumerate(SPK): cols[k]=bil(SPt[q,m],vfy,vfx,viy,vix)
        ti=gidx.get_loc(mid[m]); cols['P365']=bil(GA['P365'][ti],gfy,gfx,giy,gix); cols['P180']=bil(GA['P180'][ti],gfy,gfx,giy,gix)
        cols['stage_5km_cm90']=bil(np.asarray(EF['5km90'][edates.get_loc(mid[m])],'float32'),efy,efx,eiy,eix)
        cols['SRS180']=np.full(n,flm[m].SRS180,'float32'); cols['TSC111180']=np.full(n,flm[m].TSC111180,'float32'); Bm[m]=pred('m1',cols,F1)
    Bs={m:np.mean([Bm[k] for k in range(max(m-6,0),min(m+6,NM))],0) for m in msel}   # pandas rolling(12, center=True, min_periods=1)
    midn=mid.values.astype('datetime64[D]').astype(float)
    for t,Dt in enumerate(days):
        x=float(np.datetime64(Dt.date(),'D').astype(float)); j=int(np.clip(np.searchsorted(midn,x)-1,0,NM-2)); f=float(np.clip((x-midn[j])/(midn[j+1]-midn[j]),0,1))
        B=Bs[j]*(1-f)+Bs[j+1]*f
        ti=gidx.get_loc(Dt); di=edates.get_loc(Dt); cols=dict(st)
        for k in ['P30','P90','P180','P365','R30','T90','P30_dep','P90_dep','P180_dep']: cols[k]=bil(GA[k][ti],gfy,gfx,giy,gix)
        doy=Dt.dayofyear; cols['doy_cos']=np.full(n,np.cos(2*np.pi*doy/365.25),'float32'); cols['doy_sin']=np.full(n,np.sin(2*np.pi*doy/365.25),'float32')
        e={k:bil(np.asarray(EF[k][di],'float32'),efy,efx,eiy,eix) for k in EF}
        cols['stage_cell_cm']=e['cell']; cols['stage_cell_cm30']=e['cell30']; cols['stage_5km_cm30']=e['5km30']; cols['stage_5km_cm90']=e['5km90']
        cols['stage_cell_cm_dep']=e['cell']-e['Rcell']; cols['stage_cell_cm30_dep']=e['cell30']-e['Rcell']; cols['stage_5km_cm30_dep']=e['5km30']-e['R5']; cols['stage_5km_cm90_dep']=e['5km90']-e['R5']
        for c in DRV.columns: cols[c]=np.full(n,DRV.at[Dt,c],'float32')
        cols['src_flow90_w']=(st['w_SRS']+st['w_BigCypress'])*cols['SRS90']+(st['w_TaylorSlough']+st['w_C111'])*cols['TSC11190']
        p=np.clip(0.7*(B+pred('m2',cols,F2))+0.3*pred('m3',cols,F3),0,None)
        tile=np.full((r1-r0)*W,-1,'int16'); tile[pix-r0*W]=np.round(p*10).astype('int16'); OUT[t,r0:r1,:]=tile.reshape(r1-r0,W)
    print(f'rows {r0}-{r1} {round(time.time()-t0)}s',flush=True)
os.makedirs(f'{OUTD}/{M0.year}',exist_ok=True)
prof=dict(driver='GTiff',width=W,height=H,count=len(days),dtype='int16',nodata=-1,crs='EPSG:4326',
          transform=Affine(0.00027092271166751034,0,-81.51874748177192,0,-0.0002709227116675095,25.889645249658873),
          compress='deflate',predictor=2,tiled=True,blockxsize=256,blockysize=256,BIGTIFF='IF_SAFER')
with rasterio.open(f'{OUTD}/{M0.year}/ENP_salinity30m_daily_{M0.year}_{M0.month:02d}.tif','w',**prof) as o:
    o.write(OUT); [o.set_band_description(i+1,str(d.date())) for i,d in enumerate(days)]
    o.update_tags(units='daily salinity, psu x 10 (divide by 10); -1 = no data',model='0.7*(two-step: spectral 12-month mean + daily anomaly) + 0.3*one-step v3; 30 m')
print('done',M0.date(),round(time.time()-t0),'s')
