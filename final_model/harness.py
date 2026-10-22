"""Model experiments (one-step, station level). 10-block spatial CV (k-means on station coordinates), training on chosen years,
predictions for every held-out station-day 2000-2025. Metrics: skill by era, fresh-station bias, long-gauge trend agreement."""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import sys, time, pickle, json, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.ensemble import HistGradientBoostingRegressor as HGB
from sklearn.cluster import KMeans
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
U = f'{ROOT}'
ANCH = ['ENPWP', 'ENPSR', 'ENPWW', 'ENPGB', 'TAYLORS3', 'ENPWB']
def load():
    T = pd.read_parquet(f'{U}/dbhydro_hindcast/features/station_features_2000_2025_v2.parquet')
    T = T[(T.date >= '2000-04-01') & ~T.station.isin(ANCH)].reset_index(drop=True)
    # sea level: mean anomaly of Key West + Vaca Key + Virginia Key monthly MSL (each relative to its 2000-2025 mean), metres
    L = []
    for s in (8724580, 8723970, 8723214):
        d = pd.read_csv(f'{U}/sea_level/noaa_monthly_mean_{s}.csv'); d.columns = [c.strip() for c in d.columns]
        d['m'] = pd.to_datetime(dict(year=d.Year, month=d.Month, day=1)); d = d.set_index('m').MSL.astype(float)
        L.append(d - d['2000':'2025'].mean())
    msl = pd.concat(L, axis=1).mean(axis=1).sort_index()
    full = msl.reindex(pd.date_range('1998-01-01', '2025-12-01', freq='MS')).interpolate()
    m12 = full.rolling(12, min_periods=9).mean(); m3 = full.rolling(3, min_periods=2).mean()
    mo = T.date.dt.to_period('M').dt.to_timestamp()
    T['msl_1m'] = mo.map(full).values; T['msl_3m'] = mo.map(m3).values; T['msl_12m'] = mo.map(m12).values
    T['year_frac'] = T.date.dt.year + T.date.dt.dayofyear / 365.25
    S = pd.read_csv(f'{U}/site_variables/station_site_variables.csv').set_index('station'); st = sorted(set(T.station))
    km = KMeans(10, n_init=10, random_state=0).fit(S.loc[st, ['longitude', 'latitude']].values); T['block'] = T.station.map(dict(zip(st, km.labels_)))
    return T
from model_config import F0, MS, FE, PRM   # predictor list and hyperparameters of the final model
def run_cv(T, feats, train_years=(2000, 2025), target='raw', weight=None, prm=None, mono=None, tag=''):
    prm = dict(PRM if prm is None else prm)
    if mono: prm['monotonic_cst'] = mono
    obs = T.salinity.notna(); yr = T.date.dt.year; trn_base = obs & (yr >= train_years[0]) & (yr <= train_years[1])
    X = T[feats].astype('float32'); y = T.salinity.values
    yt = np.log1p(np.clip(y, 0, None)) if target == 'log' else (np.sqrt(np.clip(y, 0, None)) if target == 'sqrt' else y)
    p = np.full(len(T), np.nan); t0 = time.time()
    for b in range(10):
        trn = trn_base & (T.block != b); h = (T.block == b).values
        w = None if weight is None else weight[trn.values]
        m = HGB(**prm).fit(X[trn], yt[trn.values], sample_weight=w); q = m.predict(X[h])
        p[h] = np.expm1(q) if target == 'log' else (q ** 2 if target == 'sqrt' else q)
    print(f'{tag}: CV {round(time.time() - t0)} s', flush=True)
    return np.clip(p, 0, None)
def metrics(T, p, tag):
    o = T.salinity.notna().values & np.isfinite(p); y = T.salinity.values; d = T.date
    out = dict(exp=tag)
    for nm, sel in [('2000-12', d < '2013-01-01'), ('2013-25', d >= '2013-01-01')]:
        s = o & sel.values; e = p[s] - y[s]; out[f'R2_{nm}'] = 1 - np.sum(e ** 2) / np.sum((y[s] - y[s].mean()) ** 2); out[f'RMSE_{nm}'] = np.sqrt(np.mean(e ** 2)); out[f'bias_{nm}'] = e.mean()
    sm = pd.DataFrame({'st': T.station[o], 'y': y[o], 'p': p[o]}).groupby('st').mean()
    fresh = sm[sm.y < 2]; out['fresh_bias'] = (fresh.p - fresh.y).mean(); out['fresh_pred_mean'] = fresh.p.mean()
    low = sm[sm.y < 10]; out['lt10_bias'] = (low.p - low.y).mean()
    out['station_mean_R2'] = 1 - np.sum((sm.p - sm.y) ** 2) / np.sum((sm.y - sm.y.mean()) ** 2)
    # long-gauge trends (monthly means with >=15 obs days, stations with >=200 months), same months
    from gauge_trend_long_fn import sk
    D = pd.DataFrame({'station': T.station[o].values, 'm': d[o].dt.to_period('M').values, 'y': y[o], 'p': p[o]})
    g = D.groupby(['station', 'm']).agg(n=('y', 'size'), y=('y', 'mean'), p=('p', 'mean')).reset_index(); g = g[g.n >= 15]
    rows = []
    for s, gg in g.groupby('station'):
        if len(gg) < 200: continue
        gg = gg.set_index('m'); rows.append((s, sk(gg.y)[0], sk(gg.p)[0]))
    R = pd.DataFrame(rows, columns=['st', 'obs', 'mod'])
    out['n_long'] = len(R); out['trend_obs_med'] = R.obs.median(); out['trend_mod_med'] = R['mod'].median()
    out['trend_sign_agree'] = np.mean(np.sign(R.obs) == np.sign(R['mod'])); out['trend_rho'] = R.obs.corr(R['mod'], method='spearman')
    pc = D[D.station.isin(R.st)].assign(per=lambda x: np.where(x.m.dt.year < 2013, 'a', 'b'), mo=lambda x: x.m.dt.month)
    q = pc.groupby(['station', 'per', 'mo'])[['y', 'p']].mean().groupby(['station', 'per']).mean().unstack('per')
    out['change_obs'] = (q[('y', 'b')] - q[('y', 'a')]).median(); out['change_mod'] = (q[('p', 'b')] - q[('p', 'a')]).median()
    return out
