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
