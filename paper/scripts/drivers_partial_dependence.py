"""Partial dependence of the final model on main daily drivers, by station salinity zone (training rows; 3000 sampled station-days per zone)."""
import pickle, numpy as np, pandas as pd
from sklearn.inspection import partial_dependence
E = 'work/exp'; M = pickle.load(open('work/final/final_model_onestep.pkl', 'rb')); m = M['m3']; FE = M['F']
T = pd.read_parquet(f'{E}/T.parquet'); T = T[T.salinity.notna()]
z = T.groupby('station').salinity.transform('mean'); T['zone'] = pd.cut(z, [-1, 5, 25, 99], labels=['< 5 psu', '5–25 psu', '> 25 psu'])
D = {'stage_5km_cm90': 'EDEN water level, 90-day mean within 5 km (cm)', 'SRS90': 'Shark River Slough inflow, 90-day mean',
     'P180': 'precipitation, previous 180 days (mm)', 'anc_ENPGB_30': 'Florida Bay anchor (ENPGB) salinity, 30-day mean (psu)', 'msl_12m': 'sea level, 12-month mean anomaly (m)'}
rows = []
for zn, g in T.groupby('zone', observed=True):
    X = g.sample(3000, random_state=0)[FE].astype('float64').values
    for f in D:
        lo, hi = np.nanpercentile(T[f], [5, 95]); grid = np.linspace(lo, hi, 25)
        r = partial_dependence(m, X, [FE.index(f)], kind='average', method='brute', grid_resolution=25, percentiles=(0.05, 0.95))
        for v, p in zip(r['grid_values'][0], r['average'][0]): rows.append(dict(zone=zn, driver=f, value=v, salinity=p))
P = pd.DataFrame(rows); P.to_csv('partial_dependence.csv', index=False)
print(P.groupby(['driver', 'zone']).salinity.agg(lambda s: s.max() - s.min()).unstack().round(2))
for f in D: print(f, 'units range', np.nanpercentile(T[f], [5, 50, 95]).round(3))
