"""Event days per station-year at gauges vs masked maps at gauge cells (same days), and period change 2000-12 -> 2013-25."""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import sys, numpy as np, pandas as pd
U=f'{ROOT}'
D=pd.read_csv(f'{U}/dbhydro_hindcast/results/hindcast_station_daily_2000_2025.csv.gz',parse_dates=['date'])[['station','date','gauge_salinity']]
M=pd.read_csv(sys.argv[1],parse_dates=['date']); M['station']=M.station.astype(str); D['station']=D.station.astype(str)
X=D.merge(M,on=['station','date']).dropna(subset=['gauge_salinity','model_salinity_30m_onestep_masked'])
def evdays(s):
    """s: daily series (may have gaps); event = >=7 consecutive calendar days above own 95th pct"""
    q=s.quantile(0.95); s=s.asfreq('D'); ex=(s>q).astype(int).values; out=np.zeros(len(s),bool); i=0
    while i<len(ex):
        if ex[i]:
            j=i
            while j<len(ex) and ex[j]: j+=1
            if j-i>=7: out[i:j]=True
            i=j
        else: i+=1
    return pd.Series(out,index=s.index)
rows=[]
for st,g in X.groupby('station'):
    g=g.set_index('date').sort_index()
    if len(g)<730: continue
    eo=evdays(g.gauge_salinity); em=evdays(g.model_salinity_30m_onestep_masked)
    days=g.index; ny=pd.Series(1,index=days).groupby(days.year).size()
    for y,n in ny.items():
        if n<300: continue
        sel=days[days.year==y]; rows.append(dict(station=st,year=y,obs=int(eo.reindex(sel).sum()),model=int(em.reindex(sel).sum())))
E=pd.DataFrame(rows); E.to_csv(sys.argv[2],index=False)
X['m']=X.date.dt.to_period('M'); g=X.groupby(['station','m']).agg(n=('gauge_salinity','size'),y=('gauge_salinity','mean'),p=('model_salinity_30m_onestep_masked','mean')).reset_index(); g=g[g.n>=15]
C=[]
for st,gg in g.groupby('station'):
    if len(gg)<200: continue
    gg=gg.assign(per=np.where(gg.m.dt.year<2013,'a','b'),mo=gg.m.dt.month)
    q=gg.groupby(['per','mo'])[['y','p']].mean().groupby('per').mean()
    if len(q)==2: C.append(dict(station=st,obs_change=q.y['b']-q.y['a'],model_change=q.p['b']-q.p['a']))
pd.DataFrame(C).round(2).to_csv(sys.argv[3],index=False)
print('event station-years',len(E),'stations',E.station.nunique(),'r',round(np.corrcoef(E.obs,E.model)[0,1],2),'mean obs/model',round(E.obs.mean(),1),round(E.model.mean(),1))
C=pd.DataFrame(C); print('period change stations',len(C),'median obs',round(C.obs_change.median(),2),'model',round(C.model_change.median(),2))
