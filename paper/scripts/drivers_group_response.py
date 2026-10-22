"""Group response curves of the final model: for each driver group, set all its predictors at each station-day to that
station's own q-th percentile (q = 5..95) and average the predictions over sampled station-days by salinity zone."""
import pickle, numpy as np, pandas as pd
E = 'work/exp'; M = pickle.load(open('work/final/final_model_onestep.pkl', 'rb')); m = M['m3']; FE = M['F']
T = pd.read_parquet(f'{E}/T.parquet'); T = T[T.salinity.notna()].reset_index(drop=True)
z = T.groupby('station').salinity.transform('mean'); T['zone'] = pd.cut(z, [-1, 5, 25, 99], labels=['< 5 psu', '5–25 psu', '> 25 psu'])
MS = ['msl_1m', 'msl_3m', 'msl_12m']
G = {'water level': ['stage_cell_cm30', 'stage_5km_cm30', 'stage_5km_cm90'], 'flows': ['SRS90', 'TSC11190', 'src_flow90_w'],
     'climate': ['P30', 'P90', 'P180', 'P365', 'R30', 'T90'], 'anchor stations': [c for c in FE if c.startswith('anc_')], 'sea level': MS}
Q = np.arange(5, 96, 5); rows = []
S = T.groupby('zone', observed=True).sample(2500, random_state=0)
for g, cols in G.items():
    pct = T.groupby('station')[cols].quantile(Q / 100)             # (station, q) x cols
    for q in Q:
        X = S[FE].copy()
        v = pct.xs(q / 100, level=1).reindex(S.station).values; X.loc[:, cols] = np.where(np.isfinite(v), v, X[cols].values)
        p = np.clip(m.predict(X.astype('float32').values), 0, None)
        for zn in S.zone.cat.categories: rows.append(dict(group=g, q=q, zone=zn, salinity=float(p[(S.zone == zn).values].mean())))
    print(g, flush=True)
R = pd.DataFrame(rows); R.to_csv('group_response.csv', index=False)
W = R.pivot_table(index=['group', 'zone'], columns='q', values='salinity'); print((W[95] - W[5]).round(2).unstack())
