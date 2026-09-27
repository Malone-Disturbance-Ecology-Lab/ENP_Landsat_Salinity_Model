import pandas as pd, numpy as np, warnings; warnings.filterwarnings('ignore')
exec(open('spatial_cv.py').read().split("res=[]; stations={}")[0])
F0=SITE+CLIM+HYD
S=pd.read_csv('station_connectivity.csv',index_col=0)
co=a.groupby('station')[['lat','lon','distkm']].first()
def region(s):
    la,lo,d=co.loc[s]
    if d<-0.5: return 'GulfWaters' if (lo<-80.95 and la>25.2) else 'FloridaBay'
    if la>=25.55 and lo<-80.95: return 'BigCypressTTI'
    if lo<=-80.75 or (la>25.3 and lo<=-80.8): return 'SRS'
    if lo<=-80.61: return 'TaylorSlough'
    return 'C111'
S['region']=[region(s) for s in S.index]
print(S.region.value_counts().to_dict())
S['chan_dist_strict_km']=S.chan_dist_strict_km.clip(upper=100); S['chan_sinuosity']=S.chan_sinuosity.clip(upper=10)
for c in [c for c in S if c.startswith('src_dist')]: S[c]=S[c].clip(upper=150)
R=pd.get_dummies(S.region,prefix='reg').astype(float); S=S.join(R)
CONN=['chan_dist_soft_km','chan_dist_strict_km','chan_sinuosity','chan_density_500m','chan_density_1km','inland_water_1km']+[c for c in S if c.startswith('src_dist')]+list(R.columns)
for df in (tr,te):
    for c in CONN: df[c]=df.station.map(S[c])
# regional flow relevant to the station's source
fmap={'SRS':'SRS90','BigCypressTTI':'SRS90','GulfWaters':'SRS90','TaylorSlough':'TSC11190','C111':'TSC11190','FloridaBay':'TSC11190'}
for df in (tr,te):
    reg=df.station.map(S.region); df['own_source_flow90']=[r_[fmap[g]] if pd.notna(g) else np.nan for r_,g in zip(df[['SRS90','TSC11190']].to_dict('records'),reg)]
CONN+=['own_source_flow90']
# anchors
ANCH=['ENPWP','ENPSR','ENPWW','ENPGB','TAYLORS3','ENPWB']
allsal=pd.concat([tr,te])[['station','date','salinity']]
piv=allsal.pivot_table(index='date',columns='station',values='salinity').reindex(pd.date_range('2013-01-01','2025-12-31'))
A30=piv[ANCH].rolling(30,min_periods=10).mean(); A3=piv[ANCH].rolling(3,min_periods=1).mean()
for s in ANCH:
    for nm,fr in [('30',A30),('3',A3)]:
        tr[f'anc_{s}_{nm}']=fr[s].reindex(tr.date).values; te[f'anc_{s}_{nm}']=fr[s].reindex(te.date).values
AF=[f'anc_{s}_{n}' for s in ANCH for n in ['30','3']]
trn=tr[~tr.station.isin(ANCH)].copy(); ten=te[~te.station.isin(ANCH)].copy()
S.to_csv('station_connectivity_with_region.csv')
res=[]; zt={}
for name,feats in [('baseline',F0),('+ connectivity',F0+CONN),('+ anchors',F0+AF),('+ connectivity + anchors',F0+CONN+AF)]:
    p,pt=cv(trn,feats,'gbm',trn.block.values,ten)
    m,st=metrics(trn,p); ok=~np.isnan(pt); m2,st2=metrics(ten[ok],pt[ok])
    res.append(dict(model=name,RMSE=m['RMSE'],R2=m['R2'],RMSE_2025=m2['RMSE'],R2_2025=m2['R2'],station_r=m['median_station_r'],station_mean_R2=m['station_mean_R2'],station_mean_R2_2025=m2['station_mean_R2'],med_abs_bias=m['median_abs_station_bias']))
    st=st.join(sv[['distCoast_signed_m']]); st['zone']=pd.cut(st.distCoast_signed_m,[-30e3,-5e3,-500,500,3e3,10e3,30e3],labels=['open bay','bay','shoreline','0.5-3 km','3-10 km','>10 km'])
    zt[name]=st.groupby('zone',observed=True).rmse.median()
    np.save(f'cvpred_{name.replace(" ","_").replace("+","p")}.npy',np.vstack([p]))
    if name=='+ connectivity + anchors':
        trn['cvpred']=p; ten['cvpred']=pt; st.to_csv('conn_anchor_station_cv.csv')
pd.set_option('display.width',250)
print(pd.DataFrame(res).round(3).to_string(index=False)); print('\nmedian station RMSE by zone'); print(pd.DataFrame(zt).round(2).to_string())
pd.DataFrame(res).to_csv('connectivity_cv_results.csv',index=False)
trn[['station','date','salinity','cvpred','block']].to_pickle('cv_best_train.pkl'); ten[['station','date','salinity','cvpred','block']].to_pickle('cv_best_test.pkl')
