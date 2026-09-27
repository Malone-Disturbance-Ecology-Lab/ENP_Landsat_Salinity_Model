import pandas as pd, numpy as np
d=pd.read_csv('SERVER/cleaned_data/DBSAL_daily_clean.csv',parse_dates=['formatted_date'])
d=d[d.formatted_date<='2025-12-31'].sort_values(['station','formatted_date'])
rows=[]
for st,g in d.groupby('station'):
    g=g.set_index('formatted_date').asfreq('D'); g['station']=st
    p=g.precip
    for w in [30,90,180,365]: g[f'P{w}']=p.rolling(w,min_periods=int(w*.8)).sum()
    g['R30']=g.srad.rolling(30,min_periods=24).mean(); g['T90']=g.tavg.rolling(90,min_periods=72).mean()
    g['doy_cos']=np.cos(2*np.pi*g.index.dayofyear/365.25); g['doy_sin']=np.sin(2*np.pi*g.index.dayofyear/365.25)
    rows.append(g.reset_index())
a=pd.concat(rows); a=a[a.salinity.notna()].rename(columns={'formatted_date':'date'})
# flows
F=pd.read_csv('SERVER/dbhydro_flows/DBHydro_structure_flows_daily_2013_2026.csv',parse_dates=['date'],index_col='date')
for c,nm in zip(['SRS_inflow_cfs','TS_C111_detention_inflow_cfs','C111_S18C_cfs','C111_S197_cfs'],['SRS','TSC111','S18C','S197']):
    for w in [30,90,180]: F[f'{nm}{w}']=F[c].rolling(w,min_periods=int(w*.8)).mean()
FL=[f'{n}{w}' for n in ['SRS','TSC111','S18C','S197'] for w in [30,90,180]]
a=a.merge(F[FL],left_on='date',right_index=True,how='left')
# gauge stage (paired gauges)
W=pd.read_csv('SERVER/eden_stage/EDEN_gage_stage_daily_2013_2026_wide.csv',parse_dates=['date'],index_col='date')
P=pd.read_csv('SERVER/eden_stage/station_eden_gage_pairs.csv',index_col=0)
W30=W.rolling(30,min_periods=24).mean(); W90=W.rolling(90,min_periods=72).mean()
out=[]
for s,g in a.groupby('station'):
    g=g.copy()
    if s in P.index:
        p=P.loc[s]; dt=g.date
        for nm,fr in [('',W),('30',W30),('90',W90)]:
            g['H_marsh'+nm]=fr[p.marsh_gage].reindex(dt).values; g['H_tidal'+nm]=fr[p.tidal_gage].reindex(dt).values
        g['dH']=g.H_marsh-g.H_tidal; g['dH90']=g.H_marsh90-g.H_tidal90
    out.append(g)
a=pd.concat(out)
# EDEN grid at station cell
G=pd.read_csv('downloads/EDEN_grid_stage_at_stations_2013_2025.csv.gz',parse_dates=['date']).sort_values(['station','date'])
for c in ['stage_cell_cm','stage_5km_cm']:
    G[c+'30']=G.groupby('station')[c].transform(lambda s:s.rolling(30,min_periods=24).mean())
    G[c+'90']=G.groupby('station')[c].transform(lambda s:s.rolling(90,min_periods=72).mean())
a=a.merge(G,on=['station','date'],how='left')
cells=pd.read_csv('downloads/station_cells.csv').set_index('station')
st=pd.read_csv('sitevars/stations_coast.csv').set_index('station')
a['distkm']=a.station.map(st.distCoast_signed_m)/1000
a['dist_to_eden_cell_km']=a.station.map(cells.dist_to_cell_m)/1000
a['ym']=a.date.dt.to_period('M').astype(str); a['year']=a.date.dt.year; a['month']=a.date.dt.month
a.to_pickle('all_features.pkl'); print(a.shape, a.date.min(), a.date.max(), a.station.nunique())
print(cells.dist_to_cell_m.describe().round(0).to_dict(), (cells.dist_to_cell_m<=400).sum(),'stations inside a valid EDEN cell')
