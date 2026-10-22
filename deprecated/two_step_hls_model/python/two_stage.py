"""Two-step salinity model.
Stage 1: time-varying site mean B(station, month) = centered 12-month mean salinity,
         predicted from the vegetation spectral signature around the site (HLS L30
         NDVI, EVI2, NDMI, SWIR1; centered 12-mo / 3-mo composites, long-term median,
         and the 12-mo departure from the long-term median) + static site features.
Stage 2: daily anomaly (salinity - B) predicted from daily-resolution drivers
         (gridMET windows, EDEN daily stage, structure flows, anchor stations) and
         their departures from their own centered 12-month means.
Evaluation: same 10 spatial blocks and 2025 hold-out as v3."""
import pandas as pd, numpy as np, pickle, rasterio, warnings, glob, os, sys; warnings.filterwarnings('ignore')
from scipy import ndimage as ndi
from sklearn.ensemble import HistGradientBoostingRegressor
D=pickle.load(open('v3_tables.pkl','rb')); A=D['A'].copy(); B_=D['B'].copy(); FEATS_V3=D['FEATS']
A['date']=pd.to_datetime(A.date); B_['date']=pd.to_datetime(B_.date)
def gbm(): return HistGradientBoostingRegressor(max_iter=400,learning_rate=0.05,max_leaf_nodes=31,min_samples_leaf=200,l2_regularization=1.0,random_state=0)
def gbm_small(): return HistGradientBoostingRegressor(max_iter=300,learning_rate=0.05,max_leaf_nodes=15,min_samples_leaf=40,l2_regularization=1.0,random_state=0)
def metrics(df,p):
    y=df.salinity.values; r=dict(RMSE=np.sqrt(((y-p)**2).mean()),R2=1-((y-p)**2).sum()/((y-y.mean())**2).sum())
    g=pd.DataFrame({'s':df.station.values,'y':y,'p':p})
    st=g.groupby('s').apply(lambda x:pd.Series({'mean_err':x.p.mean()-x.y.mean(),'r':np.corrcoef(x.y,x.p)[0,1] if len(x)>30 and x.y.std()>0 else np.nan}))
    r|=dict(median_station_r=st.r.median(),median_abs_bias=st.mean_err.abs().median(),station_mean_R2=1-((st.mean_err)**2).sum()/((g.groupby('s').y.mean()-g.y.mean())**2).sum())
    return r
# ------------------------------------------------------------------ vegetation grid
IDX=['NDVI','EVI2','NDMI','SWIR1']; II=[0,1,2,5]
rows=list(range(4,3839,16)); V=np.load('veg_rows_L30.npz')
M=np.stack([V[f'row_{r:04d}'][:,II,:].astype('float32') for r in rows],2)   # (156,4,nr,521)
for q,(lo_,hi_) in enumerate([(-1,1),(-1,1.5),(-1,1),(0,1)]): M[:,q][(M[:,q]<lo_)|(M[:,q]>hi_)]=np.nan
# robust screen: drop monthly values > 4 MAD from the pixel's own median
med=np.nanmedian(M,0,keepdims=True); mad=np.nanmedian(np.abs(M-med),0,keepdims=True); M[np.abs(M-med)>4*1.4826*np.maximum(mad,0.01)]=np.nan
wf=rasterio.open('SERVER/site_variables/ENP_WaterFrequency_L30.tif'); T=wf.transform
NM=M.shape[0]; months=pd.date_range('2013-01-01',periods=NM,freq='MS')
dy=16*abs(T.e)*111.0; dx=8*T.a*111.0*np.cos(np.deg2rad(25.3))           # km per grid step
def disc(rkm):
    ry=int(np.ceil(rkm/dy)); rx=int(np.ceil(rkm/dx)); yy,xx=np.mgrid[-ry:ry+1,-rx:rx+1]
    return (((yy*dy)**2+(xx*dx)**2)<=rkm**2+1e-9).astype('float32')
RAD={'1km':1.0,'3km':3.0}
NB={}
for rn,rk in RAD.items():
    k=disc(rk); out=np.full(M.shape,np.nan,'float32')
    for m in range(NM):
        for i in range(len(II)):
            g=M[m,i]; v=np.isfinite(g); num=ndi.convolve(np.where(v,g,0),k,mode='constant'); den=ndi.convolve(v.astype('float32'),k,mode='constant')
            out[m,i]=np.where(den>0,num/np.maximum(den,1e-6),np.nan)
    NB[rn]=out
print('neighbourhood spectra done',{k:v.shape for k,v in NB.items()},flush=True)
def temporal(series):   # series (NM,) -> dict of composites (NM,)
    s=pd.Series(series,index=months)
    c12=s.rolling(12,center=True,min_periods=3).median(); c3=s.rolling(3,center=True,min_periods=1).median()
    lt=np.nanmedian(series) if np.isfinite(series).any() else np.nan
    return dict(c12=c12.values,c3=c3.values,lt=np.full(NM,lt),d12=c12.values-lt)
stco=pd.read_csv('sitevars/stations_coast.csv').set_index('station')
SPEC=[]; spec_rows=[]
for s in sorted(set(A.station)|set(B_.station)):
    la,lo=stco.loc[s,'latitude'],stco.loc[s,'longitude']; r_,c_=rasterio.transform.rowcol(T,lo,la)
    i=int(np.clip(round((r_-4)/16),0,len(rows)-1)); j=int(np.clip(c_//8,0,M.shape[3]-1))
    d={'station':[s]*NM,'month':months}
    for rn in RAD:
        for q,nm in enumerate(IDX):
            for kk,v in temporal(NB[rn][:,q,i,j]).items(): d[f'{nm}_{kk}_{rn}']=v
    spec_rows.append(pd.DataFrame(d))
SP=pd.concat(spec_rows); SPEC=[c for c in SP if c not in ('station','month')]
print('station spectral features',SP.shape,'missing frac',round(SP[SPEC].isna().mean().mean(),3),flush=True)
# ------------------------------------------------------------------ baseline B (centered 12-mo mean)
def baseline(df):
    mm=df.assign(month=df.date.dt.to_period('M').dt.to_timestamp()).groupby(['station','month']).salinity.agg(['sum','count']).reset_index()
    out=[]
    for s,g in mm.groupby('station'):
        g=g.set_index('month').reindex(months).fillna(0)
        S_=g['sum'].rolling(12,center=True,min_periods=1).sum(); C_=g['count'].rolling(12,center=True,min_periods=1).sum()
        nmon=(g['count']>0).astype(int).rolling(12,center=True,min_periods=1).sum()
        b=np.where((C_>=90)&(nmon>=4),S_/C_.clip(lower=1),np.nan)
        out.append(pd.DataFrame({'station':s,'month':months,'B':b}))
    return pd.concat(out)
BA=baseline(A)                     # from <=2024 data only
BT=baseline(pd.concat([A,B_]))     # used only to score stage 1 on 2025
for df in (A,B_): df['month']=df.date.dt.to_period('M').dt.to_timestamp()
A=A.merge(BA,on=['station','month'],how='left').merge(SP,on=['station','month'],how='left')
B_=B_.merge(BT,on=['station','month'],how='left').merge(SP,on=['station','month'],how='left')
# ------------------------------------------------------------------ stage-2 daily drivers + departures
AF=pd.read_pickle('all_features.pkl'); AF['date']=pd.to_datetime(AF.date)
ANCH=['ENPWP','ENPSR','ENPWW','ENPGB','TAYLORS3','ENPWB']
piv=AF.pivot_table(index='date',columns='station',values='salinity').reindex(pd.date_range('2013-01-01','2025-12-31'))
for s in ANCH:
    base=piv[s].rolling(365,center=True,min_periods=120).mean(); a30=piv[s].rolling(30,min_periods=10).mean()
    dep=(a30-base)
    for df in (A,B_): df[f'anc_{s}_dep']=dep.reindex(df.date).values
FL=['SRS30','SRS90','SRS180','TSC11130','TSC11190','TSC111180','S18C30','S18C90','S19790']
fl=AF.groupby('date')[FL].first().reindex(pd.date_range('2013-01-01','2025-12-31')).interpolate(limit=5)
for c in FL:
    dep=fl[c]-fl[c].rolling(365,center=True,min_periods=180).mean()
    for df in (A,B_): df[c]=fl[c].reindex(df.date).values; df[c+'_dep']=dep.reindex(df.date).values
LOC=['P30','P90','P180','stage_cell_cm','stage_cell_cm30','stage_5km_cm30','stage_5km_cm90']
both=pd.concat([A.assign(_s='A'),B_.assign(_s='B')]).sort_values(['station','date'])
for c in LOC:
    both[c+'_dep']=both[c]-both.groupby('station',group_keys=False).apply(lambda g:g.set_index('date')[c].rolling('365D',center=True,min_periods=60).mean().reset_index(drop=True).set_axis(g.index))
A=both[both._s=='A'].drop(columns='_s').reset_index(drop=True); B_=both[both._s=='B'].drop(columns='_s').reset_index(drop=True)
STATIC=['distkm','p_Open water_1000m','p_Mangrove_1000m','p_Freshwater marsh/prairie_1000m','p_Swamp scrub/shrubland_1000m','p_Salt marsh_1000m','p_Open water_250m','p_Mangrove_250m','dist_to_eden_cell_km',
        'chan_dist_soft_km','chan_dist_strict_km','chan_sinuosity','chan_density_500m','chan_density_1km','inland_water_1km','src_dist_SRS_km','src_dist_TaylorSlough_km','src_dist_C111_km','src_dist_BigCypressGulf_km',
        'w_SRS','w_TaylorSlough','w_C111','w_BigCypress','basin_cost_SRS_km','basin_cost_TaylorSlough_km','basin_cost_C111_km','basin_cost_BigCypress_km']
SLOW=['P365','P180','stage_5km_cm90','SRS180','TSC111180']   # slow hydro-climate state, optional in stage 1
DAILY=['P30','P90','P180','P365','R30','T90','doy_cos','doy_sin','stage_cell_cm','stage_cell_cm30','stage_5km_cm30','stage_5km_cm90']+FL+[c+'_dep' for c in FL+LOC]+\
      [f'anc_{s}_{n}' for s in ANCH for n in ['30','3','dep']]
S2F=DAILY+STATIC
pickle.dump(dict(A=A,B=B_,SPEC=SPEC,STATIC=STATIC,S2F=S2F),open('two_stage_tables.pkl','wb'))
# ------------------------------------------------------------------ CV
S1V={'static only':STATIC,'static + vegetation spectra':STATIC+SPEC,'vegetation spectra + distance':SPEC+['distkm'],'static + spectra + slow hydro-climate':STATIC+SPEC+SLOW}
def s1table(df):   # station-month table (features averaged over the month's rows)
    cols=list(dict.fromkeys(STATIC+SPEC+SLOW))
    return df.groupby(['station','month']).agg({**{c:'mean' for c in cols},'B':'first','block':'first','salinity':'count'}).rename(columns={'salinity':'n'}).reset_index()
TA=s1table(A); TB=s1table(B_)
res=[]; preds={}
for vn,f1 in S1V.items():
    pA=np.full(len(A),np.nan); pB=np.full(len(B_),np.nan); bA=np.full(len(A),np.nan); bB=np.full(len(B_),np.nan); s1A=np.full(len(TA),np.nan); s1B=np.full(len(TB),np.nan)
    for g in sorted(A.block.unique()):
        tr1=TA[(TA.block!=g)&TA.B.notna()]; m1=gbm_small().fit(tr1[f1],tr1.B,sample_weight=np.sqrt(tr1.n))
        hA=TA.block==g; hB=TB.block==g; s1A[hA.values]=m1.predict(TA.loc[hA,f1]); s1B[hB.values]=m1.predict(TB.loc[hB,f1])
        tr2=A[(A.block!=g)&A.B.notna()]; m2=gbm().fit(tr2[S2F],tr2.salinity-tr2.B)
        ha=(A.block==g).values; hb=(B_.block==g).values
        kA=A.loc[ha,['station','month']].merge(TA.loc[hA,['station','month']].assign(bh=s1A[hA.values]),how='left').bh.values
        kB=B_.loc[hb,['station','month']].merge(TB.loc[hB,['station','month']].assign(bh=s1B[hB.values]),how='left').bh.values
        aA=m2.predict(A.loc[ha,S2F]); aB=m2.predict(B_.loc[hb,S2F]) if hb.any() else []
        pA[ha]=kA+aA; bA[ha]=kA
        if hb.any(): pB[hb]=kB+aB; bB[hb]=kB
    preds[vn]=(pA,pB)
    ok1=TA.B.notna(); ok1b=TB.B.notna()
    s1r=lambda y,p:1-((y-p)**2).sum()/((y-y.mean())**2).sum()
    m=metrics(A,pA); m2=metrics(B_,pB)
    res.append(dict(stage1=vn,stage1_R2_B=s1r(TA.B[ok1].values,s1A[ok1.values]),stage1_RMSE_B=np.sqrt(((TA.B[ok1]-s1A[ok1.values])**2).mean()),stage1_R2_B_2025=s1r(TB.B[ok1b].values,s1B[ok1b.values]),
                    RMSE=m['RMSE'],R2=m['R2'],RMSE_2025=m2['RMSE'],R2_2025=m2['R2'],station_r=m['median_station_r'],station_r_2025=m2['median_station_r'],station_mean_R2=m['station_mean_R2'],station_mean_R2_2025=m2['station_mean_R2'],med_abs_bias=m['median_abs_bias']))
    print(res[-1],flush=True)
# oracle: observed B + predicted anomaly (upper bound for a perfect stage 1)
pA=np.full(len(A),np.nan)
for g in sorted(A.block.unique()):
    tr2=A[(A.block!=g)&A.B.notna()]; m2=gbm().fit(tr2[S2F],tr2.salinity-tr2.B); ha=((A.block==g)&A.B.notna()).values; pA[ha]=A.B[ha]+m2.predict(A.loc[ha,S2F])
ok=np.isfinite(pA); mo=metrics(A[ok],pA[ok]); res.append(dict(stage1='ORACLE observed B (upper bound)',RMSE=mo['RMSE'],R2=mo['R2'],station_r=mo['median_station_r'])); print(res[-1],flush=True)
R=pd.DataFrame(res); R.to_csv('two_stage_cv_results.csv',index=False); pickle.dump(preds,open('two_stage_preds.pkl','wb'))
print(R.round(3).to_string())
