import pandas as pd, numpy as np
T=pd.read_parquet('hindcast_station_predictions.parquet'); o=T[T.salinity.notna()].copy()
def pooled(y,p): return pd.Series(dict(n=len(y),stations=None,R2=1-((y-p)**2).sum()/((y-y.mean())**2).sum(),RMSE=np.sqrt(((y-p)**2).mean()),bias=(p-y).mean()))
def met(g,col):
    r=pooled(g.salinity,g[col]); r['stations']=g.station.nunique()
    s=[]
    for _,x in g.groupby('station'):
        if len(x)<50 or x.salinity.std()==0: continue
        y,p=x.salinity.values,x[col].values; rr=np.corrcoef(y,p)[0,1]
        kge=1-np.sqrt((rr-1)**2+(p.std()/y.std()-1)**2+(p.mean()/y.mean()-1)**2) if y.mean()>0 else np.nan
        s.append((rr,kge))
    s=np.array(s); r['station_r_med']=np.nanmedian(s[:,0]) if len(s) else np.nan; r['KGE_med']=np.nanmedian(s[:,1]) if len(s) else np.nan; r['n_st_50d']=len(s)
    return r
o['era']=pd.cut(o.date,pd.to_datetime(['2000-01-01','2004-01-01','2008-01-01','2013-04-01','2025-01-01','2026-01-01']),labels=['2000-03','2004-07','2008-13Q1','2013-24','2025'],right=False)
print('C) reproduction of Table 2 one-step (CV trained 2013-24):')
print(' 2013-24 held-out space:',met(o[o.era=='2013-24'],'p_cv24')[['n','R2','RMSE','station_r_med']].round(3).to_dict())
print(' 2025 held-out space+time:',met(o[o.era=='2025'],'p_cv24')[['n','R2','RMSE']].round(3).to_dict(), '| table 2: 0.786/6.23 and 0.701/6.74')
rows=[]
for col,lab in [('p_final','A released model (held out in time)'),('p_cv','B block CV (held out in space; + time before 2013)')]:
    for era,g in o.groupby('era',observed=True):
        r=met(g,col); r['model']=lab; r['era']=era; rows.append(r)
R=pd.DataFrame(rows)[['model','era','n','stations','R2','RMSE','bias','station_r_med','KGE_med','n_st_50d']]
print(R.round(3).to_string(index=False)); R.round(4).to_csv('hindcast_skill_by_era.csv',index=False)
# sondes vs grab before 2013
o['kind']=np.where(o.groupby('station').salinity.transform('size')>2000,'sonde','grab')
pre=o[o.date<'2013-04-01']
for k,g in pre.groupby('kind'): print('pre-2013',k,met(g,'p_cv')[['n','stations','R2','RMSE','bias']].round(3).to_dict())
# by basin (block CV, pre-2013 vs 2013-25)
B=pd.read_csv('/mnt/user-data/uploads/MaloneLab--Research--ENP_Landsat_Salinity_Model/best_two_step/timeseries/station_timeseries_metrics.csv').set_index('station').basin
o['basin']=o.station.map(B)
rows=[]
for (b,per),g in o.groupby(['basin',np.where(o.date<'2013-04-01','2000-13','2013-25')]):
    r=met(g,'p_cv'); rows.append(dict(basin=b,period=per,n=int(r.n),stations=r.stations,R2=r.R2,RMSE=r.RMSE,bias=r.bias))
print(pd.DataFrame(rows).round(2).to_string(index=False))
