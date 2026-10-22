"""Diagnostics quoted in Sects. 3.3, 4.1, 4.2 and 4.4 (review checks of the held-out predictions).

- share of station-days from continuous sondes vs grab-sample stations
- days on which at least one anchor-station input is missing (left as NaN; HistGradientBoosting routes them)
- held-out skill vs distance to the nearest anchor station
- residual lag-1 autocorrelation at sondes; fresh-site error vs the lowest-bin interval half-width
- held-out bias by era, pooled and per station (which stations drive the shift)
Writes anchor_distance_skill.csv and era_bias_by_station.csv; prints the numbers used in the text.
"""
import os, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr

ROOT = os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model')
WORK = os.environ.get('ENP_WORK', 'work')
AN = ['ENPWP', 'ENPSR', 'ENPWW', 'ENPGB', 'TAYLORS3', 'ENPWB']

C = pd.read_parquet(f'{WORK}/final/cv_predictions_final.parquet'); C = C[C.salinity.notna()].copy()
S = pd.read_csv(f'{ROOT}/site_variables/station_site_variables.csv').set_index('station')
C['e'] = C.pred_cv - C.salinity; C['ae'] = C.e.abs()
M = C[~C.fresh_site].copy(); M['grab'] = M.station.map(S.grab)

# 1. training composition
print('share of station-days, sondes / grab:', M.groupby('grab').size().div(len(C)).round(3).to_dict())

# 2. missing anchor inputs on map days
D = pd.read_csv(f'{WORK}/final/hpc_inputs/daily_drivers_1998_2025.csv', parse_dates=['date']).set_index('date').loc['2000-04-01':'2025-12-31']
miss = D[[c for c in D.columns if c.startswith('anc_')]].isna().any(axis=1)
print(f'map days with an anchor input missing: {int(miss.sum())} of {len(D)} ({int(miss[miss.index.year <= 2012].sum())} before 2013)')

# 3. skill vs distance to nearest anchor
def km(la1, lo1, la2, lo2): return 6371 * np.hypot(np.radians(la2 - la1), np.radians(lo2 - lo1) * np.cos(np.radians(la1)))
def kge(o, p):
    r = np.corrcoef(o, p)[0, 1]; return 1 - np.sqrt((r - 1) ** 2 + (p.std() / o.std() - 1) ** 2 + (p.mean() / o.mean() - 1) ** 2)
A = S.loc[AN, ['latitude', 'longitude']]; rows = []
for s, g in M.groupby('station'):
    if len(g) < 30 or g.salinity.std() == 0: continue
    d = km(S.loc[s, 'latitude'], S.loc[s, 'longitude'], A.latitude.values, A.longitude.values)
    rows.append(dict(station=s, dist_km=d.min(), nearest=AN[d.argmin()], r=np.corrcoef(g.salinity, g.pred_cv)[0, 1],
                     KGE=kge(g.salinity.values, g.pred_cv.values), RMSE=np.sqrt((g.e ** 2).mean())))
T = pd.DataFrame(rows); T.round(3).to_csv('anchor_distance_skill.csv', index=False)
for k in ['r', 'KGE', 'RMSE']: print(k, 'Spearman rho vs distance', round(spearmanr(T.dist_km, T[k])[0], 2))
T['band'] = pd.cut(T.dist_km, [0, 5, 15, 100], labels=['<5 km', '5-15 km', '>15 km'])
print(T.groupby('band', observed=True).agg(n=('r', 'size'), r=('r', 'median'), KGE=('KGE', 'median')).round(2))

# 4. autocorrelation and fresh-site interval width
ac = []
for s, g in M[M.grab == 0].groupby('station'):
    g = g.sort_values('date'); one = g.date.diff().dt.days.values[1:] == 1; e = g.e.values
    if one.sum() > 100: ac.append(np.corrcoef(e[:-1][one], e[1:][one])[0, 1])
q0 = json.load(open(f'{WORK}/final/q_levels.json'))['q'][0]
F = C[C.fresh_site]
print(f'median lag-1 residual autocorrelation at sondes {np.median(ac):.2f}; fresh sites mean obs {F.salinity.mean():.2f}, '
      f'median |error| {F.ae.median():.2f}, lowest-bin half-width {q0:.2f} psu')

# 5. held-out bias by era
M['era'] = np.where(M.date.dt.year <= 2012, '2000-2012', '2013-2025')
print('pooled bias by era', M.groupby('era').e.mean().round(2).to_dict())
sb = M.groupby(['station', 'era']).e.mean().unstack(); n = M.groupby(['station', 'era']).size().unstack(fill_value=0)
sb = sb[n.min(axis=1) >= 30]; sb['change'] = sb['2013-2025'] - sb['2000-2012']
sb['days'] = M.groupby('station').size(); sb['contribution'] = sb.change * sb.days / sb.days.sum()
sb.sort_values('contribution').round(3).to_csv('era_bias_by_station.csv')
print('median station change', round(sb.change.median(), 2)); print(sb.sort_values('contribution').head(6).round(2))
