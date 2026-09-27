import pandas as pd, numpy as np, warnings; warnings.filterwarnings('ignore')
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.cluster import KMeans
a=pd.read_pickle('all_features.pkl')
sv=pd.read_csv('sitevars/station_site_variables.csv').set_index('station')
vf=['p_Open water_1000m','p_Mangrove_1000m','p_Freshwater marsh/prairie_1000m','p_Swamp scrub/shrubland_1000m','p_Salt marsh_1000m','p_Open water_250m','p_Mangrove_250m']
for c in vf: a[c]=a.station.map(sv[c]).fillna(0)
unm=a.station.map(sv['mapped_frac_1000m']).fillna(0)<0.5
a.loc[unm,'p_Open water_1000m']=1; a.loc[unm,'p_Open water_250m']=1
a['lat']=a.station.map(sv.latitude); a['lon']=a.station.map(sv.longitude)
SITE=['distkm']+vf+['dist_to_eden_cell_km']
CLIM=['P30','P90','P180','P365','R30','T90','doy_cos','doy_sin']
HYD=['stage_5km_cm30','stage_cell_cm30','stage_5km_cm90','SRS90','TSC11190']
BANDS=['B02','B03','B04','B05','B06','B07','NDVI','NDSI','SRSI']
a=a[a.distkm.notna()]
tr=a[a.year<=2024].copy(); te=a[a.year==2025].copy()
# spatial blocks
co=a.groupby('station')[['lat','lon']].first()
km=KMeans(n_clusters=10,n_init=10,random_state=1).fit(co.values); blk=pd.Series(km.labels_,index=co.index)
a['block']=a.station.map(blk); tr['block']=tr.station.map(blk); te['block']=te.station.map(blk)
print('stations',a.station.nunique(),'train rows',len(tr),'2025 rows',len(te)); print('block sizes',blk.value_counts().sort_index().tolist())
def gbm(): return HistGradientBoostingRegressor(max_iter=400,learning_rate=0.05,max_leaf_nodes=31,min_samples_leaf=200,l2_regularization=1.0,random_state=0)
def hinge(df):
    d=df.distkm.clip(-25,30).values; K=[-5,-0.5,0.5,3,10]
    B=np.column_stack([d]+[np.maximum(d-k,0) for k in K]); s=df[['doy_cos','doy_sin']].fillna(0).values
    X=np.column_stack([B,B*s[:,[0]],B*s[:,[1]],df[vf].values,df[CLIM+HYD].fillna(df[CLIM+HYD].median()).values])
    return np.column_stack([np.ones(len(df)),X])
def metrics(df,p):
    y=df.salinity.values; r=dict(RMSE=np.sqrt(((y-p)**2).mean()),MAE=np.abs(y-p).mean(),R2=1-((y-p)**2).sum()/((y-y.mean())**2).sum())
    g=pd.DataFrame({'s':df.station.values,'y':y,'p':p})
    st=g.groupby('s').apply(lambda x:pd.Series({'mean_err':x.p.mean()-x.y.mean(),'r':np.corrcoef(x.y,x.p)[0,1] if len(x)>30 and x.y.std()>0 else np.nan,'rmse':np.sqrt(((x.y-x.p)**2).mean())}))
    r|=dict(median_station_r=st.r.median(),median_abs_station_bias=st.mean_err.abs().median(),station_mean_R2=1-((st.mean_err)**2).sum()/((g.groupby('s').y.mean()-g.y.mean())**2).sum())
    return r,st
def cv(df,feats,kind,groups,testdf=None):
    pred=np.full(len(df),np.nan); predt=np.full(len(testdf),np.nan) if testdf is not None else None
    for gv in np.unique(groups):
        trm=groups!=gv; tem=~trm
        if kind=='gbm':
            m=gbm().fit(df.loc[trm,feats],df.salinity[trm]); pred[tem]=m.predict(df.loc[tem,feats])
            if testdf is not None:
                tm=(testdf.station.isin(df.station[tem].unique())).values
                if tm.any(): predt[tm]=m.predict(testdf.loc[tm,feats])
        else:
            Xtr=hinge(df[trm]); b=np.linalg.lstsq(Xtr,df.salinity[trm],rcond=None)[0]; pred[tem]=hinge(df[tem])@b
            if testdf is not None:
                tm=(testdf.station.isin(df.station[tem].unique())).values
                if tm.any(): predt[tm]=hinge(testdf[tm])@b
    return pred,predt
res=[]; stations={}
F0=SITE+CLIM+HYD
for name,kind,feats,grp in [('linear (dist×season+veg+clim+EDEN)','lin',None,'block'),('GBM no bands, spatial blocks','gbm',F0,'block'),('GBM no bands, leave-one-station-out','gbm',F0,'station'),('GBM + lat/lon, spatial blocks','gbm',F0+['lat','lon'],'block')]:
    g=tr[grp].values if grp=='block' else pd.factorize(tr.station)[0]
    p,pt=cv(tr,feats,kind,g,te)
    m,st=metrics(tr,p); res.append(dict(model=name,eval='held-out stations 2013-24',**m)); stations[name]=st
    ok=~np.isnan(pt); m2,_=metrics(te[ok],pt[ok]); res.append(dict(model=name,eval='held-out stations, 2025',**m2))
# bands: clear-image days only, compare with and without bands on same rows
cl=tr.dropna(subset=['B03']).copy(); clt=te.dropna(subset=['B03']).copy()
for name,feats in [('GBM no bands (clear days)',F0),('GBM + Landsat bands (clear days)',F0+BANDS),('GBM bands+site only, no clim/EDEN (clear days)',SITE+BANDS+['doy_cos','doy_sin'])]:
    p,pt=cv(cl,feats,'gbm',cl.block.values,clt)
    m,_=metrics(cl,p); res.append(dict(model=name,eval='held-out stations 2013-24',**m))
    ok=~np.isnan(pt); m2,_=metrics(clt[ok],pt[ok]); res.append(dict(model=name,eval='held-out stations, 2025',**m2))
R=pd.DataFrame(res); pd.set_option('display.width',250)
print(R.round(3).to_string(index=False)); R.to_csv('spatial_cv_results.csv',index=False)
st=stations['GBM no bands, spatial blocks'].join(sv[['distCoast_signed_m']])
bins=[-30e3,-5e3,-500,500,3e3,10e3,30e3]; lab=['open bay <-5 km','bay -5 to -0.5','shoreline ±0.5','0.5-3 km inland','3-10 km','>10 km']
st['zone']=pd.cut(st.distCoast_signed_m,bins,labels=lab)
print('\nGBM spatial-block CV by zone (station level)'); print(st.groupby('zone',observed=True).agg(n=('rmse','size'),median_rmse=('rmse','median'),median_bias=('mean_err','median'),median_abs_bias=('mean_err',lambda x:x.abs().median()),median_r=('r','median')).round(2).to_string())
print('\nworst stations:'); print(st.sort_values('rmse',ascending=False).head(10).round(2).to_string())
st.to_csv('spatial_cv_station_gbm.csv')
