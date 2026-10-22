"""Final models (all non-anchor stations, 2013-2025) + interval calibration from out-of-block residuals."""
exec(open('ts_setup.py').read())
from sklearn.ensemble import HistGradientBoostingRegressor as H
S1CFG=dict(max_iter=700,learning_rate=0.03,max_leaf_nodes=4,min_samples_leaf=400,l2_regularization=3.0,random_state=0)
F1=STATIC+SPEC+SLOW
# observed centered 12-month mean from ALL data (2013-2025)
ALL=pd.concat([A,B_],ignore_index=True)
months=pd.date_range('2013-01-01',periods=156,freq='MS')
mm=ALL.groupby(['station','month']).salinity.agg(['sum','count']).reset_index(); out=[]
for s,g in mm.groupby('station'):
    g=g.set_index('month').reindex(months).fillna(0)
    S_=g['sum'].rolling(12,center=True,min_periods=1).sum(); C_=g['count'].rolling(12,center=True,min_periods=1).sum(); nm=(g['count']>0).astype(int).rolling(12,center=True,min_periods=1).sum()
    out.append(pd.DataFrame({'station':s,'month':months,'Bfull':np.where((C_>=90)&(nm>=4),S_/C_.clip(lower=1),np.nan)}))
ALL=ALL.merge(pd.concat(out),on=['station','month'],how='left')
TT=ALL.groupby(['station','month']).agg({**{c:'mean' for c in F1},'Bfull':'first','salinity':'count'}).rename(columns={'salinity':'n'}).reset_index()
t1=TT[TT.Bfull.notna()]; m1=H(**S1CFG).fit(t1[F1],t1.Bfull,sample_weight=np.sqrt(t1.n))
t2=ALL[ALL.Bfull.notna()]; m2=gbm().fit(t2[S2F],t2.salinity-t2.Bfull)
m3=gbm().fit(ALL[FEATS_V3],ALL.salinity)
print('final fits: stage1 rows',len(t1),'stage2 rows',len(t2),'v3 rows',len(ALL),flush=True)
# ---- interval calibration (split conformal on out-of-block residuals of the blend, 2013-2024)
d=pickle.load(open('best_two_step_oof.pkl','rb')); P=d['P']
W={'K1 two-step':0.7,'v3 one-step':0.3}
kA=d['keyA'].copy(); kB=d['keyB'].copy(); kA['p']=sum(w*P[k][0] for k,w in W.items()); kB['p']=sum(w*P[k][1] for k,w in W.items())
kA['p']=kA.p.clip(lower=0); kB['p']=kB.p.clip(lower=0)
kA['e']=(kA.salinity-kA.p).abs(); kB['e']=(kB.salinity-kB.p).abs()
bins=[-1,5,15,30,1e9]; lab=['0-5','5-15','15-30','>30']
for k_ in (kA,kB): k_['lev']=pd.cut(k_.p,bins,labels=lab)
schemes={'basin':['basin'],'level':['lev'],'basin x level':['basin','lev']}
res=[]; Q={}
for nm,keys in schemes.items():
    q=kA.dropna(subset=['p']).groupby(keys,observed=True).e.quantile(0.9); Q[nm]=q
    t=kB.dropna(subset=['p']); qq=t.set_index(keys).index.map(q.to_dict()) if len(keys)>1 else t[keys[0]].map(q)
    qq=np.asarray(qq,dtype=float); qq=np.where(np.isfinite(qq),qq,kA.e.quantile(0.9))
    lo=np.clip(t.p-qq,0,None); hi=t.p+qq; cov=(t.salinity>=lo)&(t.salinity<=hi); sc=cov.groupby(t.station).mean()
    res.append(dict(scheme=nm,cov_2025=cov.mean(),median_width=float(np.median(hi-lo)),station_p10=sc.quantile(.1),n_below70=int((sc<0.7).sum())))
R=pd.DataFrame(res); print(R.round(3).to_string()); R.to_csv('final_interval_calibration.csv',index=False)
best=R.sort_values(['n_below70','median_width']).scheme.iloc[0]; print('chosen interval scheme:',best)
pickle.dump(dict(m1=m1,m2=m2,m3=m3,F1=F1,S2F=S2F,FEATS_V3=FEATS_V3,W=W,Q=Q[best],Qkeys=schemes[best],Qall=float(kA.e.quantile(0.9)),bins=bins,lab=lab,S1CFG=S1CFG),open('final_models.pkl','wb'))
print('saved')
