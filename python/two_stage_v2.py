"""Two-step model, round 2: stronger / smoother stage-1 mean, and a stacked variant.
Reuses two_stage_tables.pkl written by two_stage.py."""
import pandas as pd, numpy as np, pickle, warnings; warnings.filterwarnings('ignore')
from sklearn.ensemble import HistGradientBoostingRegressor
exec(open('two_stage.py').read().split('# ------------------------------------------------------------------ vegetation grid')[0].split('D=pickle.load')[0])
exec('def gbm'+open('two_stage.py').read().split('def gbm',1)[1].split('# ------------------------------------------------------------------ vegetation grid')[0])
T=pickle.load(open('two_stage_tables.pkl','rb')); A=T['A']; B_=T['B']; SPEC=T['SPEC']; STATIC=T['STATIC']; S2F=T['S2F']
FEATS_V3=pickle.load(open('v3_tables.pkl','rb'))['FEATS']
ANCH=['ENPWP','ENPSR','ENPWW','ENPGB','TAYLORS3','ENPWB']
# anchor 12-month centered means (regional salinity state, observed daily at 6 sites)
AF=pd.read_pickle('all_features.pkl'); AF['date']=pd.to_datetime(AF.date)
piv=AF.pivot_table(index='date',columns='station',values='salinity').reindex(pd.date_range('2013-01-01','2025-12-31'))
ANCB=[f'anc_{s}_B' for s in ANCH]
for s in ANCH:
    b=piv[s].rolling(365,center=True,min_periods=120).mean()
    for df in (A,B_): df[f'anc_{s}_B']=b.reindex(df.date).values
# vegetation seasonal amplitude (12-mo centered 80th-20th pct of the 3-mo composite) and 12-mo trend
for df in (A,B_): pass
SPc=pd.concat([A[['station','month']+SPEC],B_[['station','month']+SPEC]]).drop_duplicates(['station','month']).sort_values(['station','month'])
AMP=[]
for c in [c for c in SPEC if '_c3_' in c]:
    n=c.replace('_c3_','_amp_'); g=SPc.groupby('station')[c]
    SPc[n]=g.transform(lambda x:x.rolling(12,center=True,min_periods=6).quantile(.8)-x.rolling(12,center=True,min_periods=6).quantile(.2)); AMP.append(n)
A=A.merge(SPc[['station','month']+AMP],on=['station','month'],how='left'); B_=B_.merge(SPc[['station','month']+AMP],on=['station','month'],how='left')
SPEC_NOC3=[c for c in SPEC if '_c3_' not in c]
SLOW=['P365','P180','stage_5km_cm90','SRS180','TSC111180']
cols=list(dict.fromkeys(STATIC+SPEC+AMP+SLOW+ANCB))
def s1table(df): return df.groupby(['station','month']).agg({**{c:'mean' for c in cols},'B':'first','block':'first','salinity':'count'}).rename(columns={'salinity':'n'}).reset_index()
TA=s1table(A); TB=s1table(B_)
def smooth(tab,col):   # 12-month centered rolling mean of the predicted mean, per station
    t=tab[['station','month',col]].sort_values(['station','month'])
    t[col+'_s']=t.groupby('station',group_keys=False).apply(lambda g:g.set_index('month')[col].rolling('365D',center=True,min_periods=1).mean().set_axis(g.index))
    return tab.join(t[[col+'_s']])[col+'_s'].values
def fit_s1(tr,f): return gbm_small().fit(tr[f],tr.B,sample_weight=np.sqrt(tr.n))
def s1_oof(f,blocks_tr):   # out-of-block stage-1 predictions for station-months in blocks_tr (nested)
    out=pd.Series(np.nan,index=TA.index)
    for g in blocks_tr:
        tr=TA[TA.block.isin(blocks_tr)&(TA.block!=g)&TA.B.notna()]; h=TA.block==g; out[h]=fit_s1(tr,f).predict(TA.loc[h,f])
    return out.to_numpy(copy=True)
def attach(df,tab,vals):
    return df[['station','month']].merge(tab[['station','month']].assign(v=vals),how='left').v.values
s1r=lambda y,p:1-((y-p)**2).sum()/((y-y.mean())**2).sum()
V={'S1: static+spectra+slow (smoothed)':STATIC+SPEC+SLOW,
   'S1: static+spectra(no 3-mo)+amp+slow+anchor means (smoothed)':STATIC+SPEC_NOC3+AMP+SLOW+ANCB,
   'S1: same, WITHOUT vegetation spectra (smoothed)':STATIC+SLOW+ANCB}
res=[]; keep={}
blocks=sorted(A.block.unique())
import os
for vn,f1 in ([] if os.environ.get('STACK_ONLY') else V.items()):
    pA=np.full(len(A),np.nan); pB=np.full(len(B_),np.nan); s1A=np.full(len(TA),np.nan); s1B=np.full(len(TB),np.nan)
    for g in blocks:
        m1=fit_s1(TA[(TA.block!=g)&TA.B.notna()],f1); hA=(TA.block==g).values; hB=(TB.block==g).values
        s1A[hA]=m1.predict(TA.loc[hA,f1]); s1B[hB]=m1.predict(TB.loc[hB,f1])
    TA['bh']=s1A; TB['bh']=s1B
    both=pd.concat([TA[['station','month','bh']],TB[['station','month','bh']]]).reset_index(drop=True); both['bhs']=smooth(both,'bh')
    sA=both.iloc[:len(TA)].bhs.values; sB=both.iloc[len(TA):].bhs.values
    for g in blocks:
        tr2=A[(A.block!=g)&A.B.notna()]; m2=gbm().fit(tr2[S2F],tr2.salinity-tr2.B)
        ha=(A.block==g).values; hb=(B_.block==g).values
        pA[ha]=attach(A[ha],TA,sA)+m2.predict(A.loc[ha,S2F])
        if hb.any(): pB[hb]=attach(B_[hb],TB,sB)+m2.predict(B_.loc[hb,S2F])
    ok1=TA.B.notna().values; ok1b=TB.B.notna().values; m=metrics(A,pA); mb=metrics(B_,pB)
    res.append(dict(model='two-step | '+vn,stage1_R2_B=s1r(TA.B.values[ok1],sA[ok1]),stage1_R2_B_2025=s1r(TB.B.values[ok1b],sB[ok1b]),RMSE=m['RMSE'],R2=m['R2'],RMSE_2025=mb['RMSE'],R2_2025=mb['R2'],station_r=m['median_station_r'],station_mean_R2=m['station_mean_R2'],station_mean_R2_2025=mb['station_mean_R2'],med_abs_bias=m['median_abs_bias']))
    keep[vn]=(pA,pB); print(res[-1],flush=True)
# ---- stacked: v3 one-step GBM + nested out-of-block spectral stage-1 mean as a feature
f1=V['S1: static+spectra+slow (smoothed)']
pA=np.full(len(A),np.nan); pB=np.full(len(B_),np.nan)
for g in blocks:
    trb=[b for b in blocks if b!=g]
    oof=s1_oof(f1,trb)                                    # training rows: nested OOF
    m1=fit_s1(TA[(TA.block!=g)&TA.B.notna()],f1); hA=(TA.block==g).values; hB=(TB.block==g).values
    oof[hA]=m1.predict(TA.loc[hA,f1]); bB=np.full(len(TB),np.nan); bB[hB]=m1.predict(TB.loc[hB,f1])
    both=pd.concat([TA[['station','month']].assign(bh=oof),TB[['station','month']].assign(bh=bB)]).reset_index(drop=True); both['bhs']=smooth(both,'bh')
    A['Bhat']=attach(A,TA,both.iloc[:len(TA)].bhs.values); B_['Bhat']=attach(B_,TB,both.iloc[len(TA):].bhs.values)
    F=FEATS_V3+['Bhat']; tr=A[A.block!=g]; m=gbm().fit(tr[F],tr.salinity)
    ha=(A.block==g).values; hb=(B_.block==g).values; pA[ha]=m.predict(A.loc[ha,F])
    if hb.any(): pB[hb]=m.predict(B_.loc[hb,F])
    print('stack fold',g,flush=True)
m=metrics(A,pA); mb=metrics(B_,pB)
res.append(dict(model='stacked | v3 features + spectral stage-1 mean',RMSE=m['RMSE'],R2=m['R2'],RMSE_2025=mb['RMSE'],R2_2025=mb['R2'],station_r=m['median_station_r'],station_mean_R2=m['station_mean_R2'],station_mean_R2_2025=mb['station_mean_R2'],med_abs_bias=m['median_abs_bias']))
keep['stacked']=(pA,pB); print(res[-1],flush=True)
v3=pd.read_csv('v3_cv_results.csv').iloc[[-1]].rename(columns={'station_r':'station_r'}); v3['model']='v3 one-step (reference)'
R=pd.concat([v3,pd.DataFrame(res)],ignore_index=True); R.to_csv('two_stage_v2_cv_results'+('_stack' if os.environ.get('STACK_ONLY') else '')+'.csv',index=False); pickle.dump(keep,open('two_stage_v2_preds.pkl','wb'))
print(R[['model','stage1_R2_B','stage1_R2_B_2025','RMSE','R2','RMSE_2025','R2_2025','station_r','station_mean_R2','station_mean_R2_2025']].round(3).to_string())
