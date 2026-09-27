"""Gauge check of the monthly-anomaly trend: seasonal Kendall (Hirsch-Slack) + seasonal Theil-Sen on monthly mean
observed salinity vs the 30 m model on the same days, at stations with >= 10 years of data."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from trend30_seasonal import seasonal_kendall
S=pd.read_csv('/mnt/user-data/outputs/salinity30m/station_daily_model_series_2013_2025_250m_30m.csv.gz',parse_dates=['date'])
S=S.dropna(subset=['gauge_salinity']); S['y']=S.date.dt.year; S['m']=S.date.dt.month
ANCH=['ENPWP','ENPSR','ENPWW','ENPGB','TAYLORS3','ENPWB']
M=S.groupby(['station','y','m']).agg(obs=('gauge_salinity','mean'),mod=('model_salinity_30m','mean'),n=('gauge_salinity','size')).reset_index()
M=M[M.n>=10]
rows=[]
for st,g in M.groupby('station'):
    if g.y.nunique()<10: continue
    X=np.full((13,12,2),np.nan,'float32')
    X[g.y.values-2013,g.m.values-1,0]=g.obs.values; X[g.y.values-2013,g.m.values-1,1]=g['mod'].values
    o=seasonal_kendall(X)
    d=dict(station=st,anchor=st in ANCH,years=g.y.nunique(),months=len(g))
    for k in o:
        d[f'obs_slope_{k}']=o[k]['slope'][0]; d[f'obs_p_{k}']=o[k]['p'][0]; d[f'mod_slope_{k}']=o[k]['slope'][1]; d[f'mod_p_{k}']=o[k]['p'][1]
    rows.append(d)
R=pd.DataFrame(rows); R.round(4).to_csv('change30/station_trend_check_monthly_anomaly_30m.csv',index=False)
for k in ['monthly','drymonths','wetmonths']:
    na=R[~R.anchor]
    print(k,'stations',len(R),'obs sig up/down',int(((R[f'obs_p_{k}']<.05)&(R[f'obs_slope_{k}']>0)).sum()),int(((R[f'obs_p_{k}']<.05)&(R[f'obs_slope_{k}']<0)).sum()),
          'model sig up/down',int(((R[f'mod_p_{k}']<.05)&(R[f'mod_slope_{k}']>0)).sum()),int(((R[f'mod_p_{k}']<.05)&(R[f'mod_slope_{k}']<0)).sum()),
          'median slope obs/mod',round(R[f'obs_slope_{k}'].median(),3),round(R[f'mod_slope_{k}'].median(),3),
          'sign agree',round((np.sign(R[f'obs_slope_{k}'])==np.sign(R[f'mod_slope_{k}'])).mean(),2),
          'rho',round(spearmanr(R[f'obs_slope_{k}'],R[f'mod_slope_{k}'])[0],2),'rho non-anchor',round(spearmanr(na[f'obs_slope_{k}'],na[f'mod_slope_{k}'])[0],2))
print(R[['station','years','obs_slope_monthly','obs_p_monthly','mod_slope_monthly','mod_p_monthly']][R.obs_p_monthly<0.05].round(3).to_string())
