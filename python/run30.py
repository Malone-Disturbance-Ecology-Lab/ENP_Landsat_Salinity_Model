"""Full two-step model at 30 m (native HLS/Landsat grid) on the 15th of every month, 2013-2025.
All inputs resampled to 30 m: static layers native 30 m; vegetation spectra, EDEN (400 m) and gridMET (4 km) bilinear.
Output: map30/C_int16.npy  (156, 3839, 4171)  = 30 m full-model prediction minus the 250 m daily map (bilinear) on that day,
        in psu*10 (-32768 = no data). Daily 30 m salinity = bilinear(250 m daily) + C interpolated in time.
Tiles of rows are processed independently and the run resumes from map30/tiles_done.txt."""
import numpy as np, pandas as pd, json, os, time, pickle, warnings; warnings.filterwarnings('ignore')
import tl2cgen
from fillers import laplace_solver
OUT='map30'; os.makedirs(OUT,exist_ok=True); t0=time.time()
H,W=3839,4171; NM=156; TILE=96
FM=pickle.load(open('final_models.pkl','rb')); F1=list(FM['m1'].feature_names_in_); F2=list(FM['m2'].feature_names_in_); F3=list(FM['m3'].feature_names_in_)
P={k:tl2cgen.Predictor(f'./{k}.so',nthread=2) for k in ('m1','m2','m3')}
def pred(k,cols,F): return P[k].predict(tl2cgen.DMatrix(np.column_stack([cols[c] for c in F]).astype('float32'))).ravel()
idx=np.load('px30_index.npy'); STK=json.load(open('static30_keys.json')); STA=np.load('static30.npy',mmap_mode='r')
LL=np.load('lonlat30.npy',mmap_mode='r'); UT=np.load('utm30.npy',mmap_mode='r')
SPK=json.load(open('spec_grid_keys.json')); SPG=np.load('spec_grid.npy',mmap_mode='r')
MED=np.load('daily_maps/median_int16.npy',mmap_mode='r')
months=pd.date_range('2013-01-01',periods=NM,freq='MS'); mid=months+pd.Timedelta(days=14); edates=pd.date_range('2013-01-01','2025-12-31')
# ---- gridMET (identical to the 250 m run)
G=np.load('gridmet_features_2012_2025.npz'); glat=G['lat']; glon=G['lon']; gd=pd.to_datetime(G['dates'])
flip=glat[0]>glat[-1]; glat0=glat[::-1] if flip else glat; step=np.diff(glat0).mean(); npad=7; glat_s=np.r_[glat0[0]-step*np.arange(npad,0,-1),glat0]
CLIMV=['P30','P90','P180','P365','R30','T90']; GA={}; gidx=pd.Index(gd); tis=[gidx.get_loc(d) for d in mid]
for k in CLIMV:
    A_=G[k][:,::-1,:] if flip else G[k]; A_=np.concatenate([np.full((A_.shape[0],npad,A_.shape[2]),np.nan,'float32'),A_],1)
    kn=np.isfinite(A_[-1]); GA[k]=laplace_solver(kn)(np.where(kn[None],A_,0)).astype('float32'); GA[k][~np.isfinite(A_).any(axis=(1,2))]=np.nan
def cen(x,w=365,minp=60):
    v=np.isfinite(x); c=np.concatenate([np.zeros((1,)+x.shape[1:]),np.cumsum(np.where(v,x,0),0,dtype='float64')]); n=np.concatenate([np.zeros((1,)+x.shape[1:]),np.cumsum(v,0)])
    i=np.arange(len(x)); hi=np.minimum(i+w//2+1,len(x)); lo=np.maximum(i-w//2,0); return np.where(n[hi]-n[lo]>=minp,(c[hi]-c[lo])/np.maximum(n[hi]-n[lo],1),np.nan).astype('float32')
for k in ['P30','P90','P180']: GA[k+'_dep']=GA[k]-cen(GA[k])
GAm={k:v[tis] for k,v in GA.items()}; del GA,G          # mid-month fields only
# ---- EDEN daily fields at mid-month days
meta=np.load('eden_daily/crop_meta.npy'); ey0,ex0,eh,ew=meta; eh,ew=int(eh),int(ew); dis=[edates.get_loc(d) for d in mid]
EF={k:np.asarray(np.load(f'eden_daily/{k}.npy',mmap_mode='r')[dis],dtype='float32') for k in ['cell','cell30','Rcell','5km30','5km90','R5']}
# ---- flows + anchors
AF=pd.read_pickle('all_features.pkl'); AF['date']=pd.to_datetime(AF.date)
FL=['SRS30','SRS90','SRS180','TSC11130','TSC11190','TSC111180','S18C30','S18C90','S19790']
fl=AF.groupby('date')[FL].first().reindex(edates).interpolate(limit=5)
for c in FL: fl[c+'_dep']=fl[c]-fl[c].rolling(365,center=True,min_periods=180).mean()
ANCH=['ENPWP','ENPSR','ENPWW','ENPGB','TAYLORS3','ENPWB']; piv=AF.pivot_table(index='date',columns='station',values='salinity').reindex(edates); anc=pd.DataFrame(index=edates)
for s in ANCH:
    a30=piv[s].rolling(30,min_periods=10).mean(); anc[f'anc_{s}_30']=a30; anc[f'anc_{s}_3']=piv[s].rolling(3,min_periods=1).mean(); anc[f'anc_{s}_dep']=a30-piv[s].rolling(365,center=True,min_periods=120).mean()
flm=[fl.loc[months[m]:months[m]+pd.offsets.MonthEnd(0),['SRS180','TSC111180']].mean() for m in range(NM)]
print('setup',round(time.time()-t0),flush=True)
C=np.lib.format.open_memmap(f'{OUT}/C_int16.npy','r+' if os.path.exists(f'{OUT}/C_int16.npy') else 'w+','int16',(NM,H,W))
donef=f'{OUT}/tiles_done.txt'; done=set(open(donef).read().split()) if os.path.exists(donef) else set()
def bil(A2,fy,fx,iy,ix):   # NaN-aware bilinear on a 2-D grid
    v=[A2[iy,ix],A2[iy,ix+1],A2[iy+1,ix],A2[iy+1,ix+1]]; w=[(1-fy)*(1-fx),(1-fy)*fx,fy*(1-fx),fy*fx]
    num=sum(np.where(np.isfinite(a),a*b,0) for a,b in zip(v,w)); den=sum(np.where(np.isfinite(a),b,0) for a,b in zip(v,w))
    return np.where(den>0,num/np.maximum(den,1e-9),np.nan).astype('float32')
for r0 in range(0,H,TILE):
    r1=min(r0+TILE,H)
    if str(r0) in done: continue
    a,b=np.searchsorted(idx,[r0*W,r1*W]); n=b-a
    if n==0: open(donef,'a').write(f'{r0}\n'); continue
    tt=time.time(); pix=idx[a:b]; rr,cc=pix//W,pix%W
    st={k:np.asarray(STA[i,a:b]) for i,k in enumerate(STK)}
    lon,lat=np.asarray(LL[0,a:b]),np.asarray(LL[1,a:b]); X,Y=np.asarray(UT[0,a:b]),np.asarray(UT[1,a:b])
    # gridMET / EDEN / 250 m map / spectral bilinear weights
    giy=np.clip(np.searchsorted(glat_s,lat)-1,0,len(glat_s)-2); gix=np.clip(np.searchsorted(glon,lon)-1,0,len(glon)-2)
    gfy=np.clip((lat-glat_s[giy])/(glat_s[giy+1]-glat_s[giy]),0,1); gfx=np.clip((lon-glon[gix])/(glon[gix+1]-glon[gix]),0,1)
    ry=(Y-ey0)/400; rx=(X-ex0)/400; eiy=np.clip(np.floor(ry).astype(int),0,eh-2); eix=np.clip(np.floor(rx).astype(int),0,ew-2); efy=np.clip(ry-eiy,0,1); efx=np.clip(rx-eix,0,1)
    fi=(rr-4)/8.0; miy=np.clip(np.floor(fi).astype(int),0,478); mfy=np.clip(fi-miy,0,1); fj=(cc-4)/8.0; mix=np.clip(np.floor(fj).astype(int),0,519); mfx=np.clip(fj-mix,0,1)
    vi=(rr-4)/16.0; va=max(int(np.floor(vi.min()))-1,0); vb=min(int(np.floor(vi.max()))+3,SPG.shape[2])
    viy=np.clip(np.floor(vi).astype(int),va,vb-2)-va; vfy=np.clip(vi-(viy+va),0,1); vj=(cc-4)/8.0; vix=np.clip(np.floor(vj).astype(int),0,SPG.shape[3]-2); vfx=np.clip(vj-vix,0,1)
    SPt=np.asarray(SPG[:,:,va:vb,:],dtype='float32')
    # ---- step 1: monthly site mean at 30 m, smoothed over 12 months
    Bm=np.empty((NM,n),'float32')
    for m in range(NM):
        cols=dict(st)
        for q,k in enumerate(SPK): cols[k]=bil(SPt[q,m],vfy,vfx,viy,vix)
        cols['P365']=bil(GAm['P365'][m],gfy,gfx,giy,gix); cols['P180']=bil(GAm['P180'][m],gfy,gfx,giy,gix)
        cols['stage_5km_cm90']=bil(EF['5km90'][m],efy,efx,eiy,eix)
        cols['SRS180']=np.full(n,flm[m].SRS180,'float32'); cols['TSC111180']=np.full(n,flm[m].TSC111180,'float32')
        Bm[m]=pred('m1',cols,F1)
    cs=np.concatenate([np.zeros((1,n)),np.cumsum(Bm,0,dtype='float64')]); ii=np.arange(NM); hi=np.minimum(ii+6,NM); lo=np.maximum(ii-6,0)
    Bs=((cs[hi]-cs[lo])/(hi-lo)[:,None]).astype('float32')        # centered 12-month window (pandas center=True: i-6 .. i+5)
    # ---- step 2 + blend on the 15th of each month; detail = full30 - bilinear(250 m map)
    for m in range(NM):
        Dt=mid[m]; di=dis[m]; cols=dict(st)
        for k in CLIMV+['P30_dep','P90_dep','P180_dep']: cols[k]=bil(GAm[k][m],gfy,gfx,giy,gix)
        doy=Dt.dayofyear; cols['doy_cos']=np.full(n,np.cos(2*np.pi*doy/365.25),'float32'); cols['doy_sin']=np.full(n,np.sin(2*np.pi*doy/365.25),'float32')
        e={k:bil(EF[k][m],efy,efx,eiy,eix) for k in EF}
        cols['stage_cell_cm']=e['cell']; cols['stage_cell_cm30']=e['cell30']; cols['stage_5km_cm30']=e['5km30']; cols['stage_5km_cm90']=e['5km90']
        cols['stage_cell_cm_dep']=e['cell']-e['Rcell']; cols['stage_cell_cm30_dep']=e['cell30']-e['Rcell']; cols['stage_5km_cm30_dep']=e['5km30']-e['R5']; cols['stage_5km_cm90_dep']=e['5km90']-e['R5']
        for c in fl.columns: cols[c]=np.full(n,fl.at[Dt,c],'float32')
        for c in anc.columns: cols[c]=np.full(n,anc.at[Dt,c],'float32')
        cols['src_flow90_w']=(st['w_SRS']+st['w_BigCypress'])*cols['SRS90']+(st['w_TaylorSlough']+st['w_C111'])*cols['TSC11190']
        full=np.clip(0.7*(Bs[m]+pred('m2',cols,F2))+0.3*pred('m3',cols,F3),0,None)
        g=np.asarray(MED[di],dtype='float32')/10; g[g<0]=np.nan
        b250=bil(g,mfy,mfx,miy,mix)
        d=np.where(np.isfinite(b250),np.clip(np.round((full-b250)*10),-30000,30000),-32768).astype('int16')
        tile=np.full((r1-r0)*W,-32768,'int16'); tile[pix-r0*W]=d; C[m,r0:r1,:]=tile.reshape(r1-r0,W)
    C.flush(); open(donef,'a').write(f'{r0}\n')
    print(f'tile rows {r0}-{r1} px {n} {round(time.time()-tt)}s total {round(time.time()-t0)}s',flush=True)
print('done')
