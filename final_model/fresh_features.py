"""Model inputs (one-step feature set, 53 predictors) at the 13 new DBHydro freshwater stations on their sample dates,
built with the same code as the maps (onestep_maps.py features())."""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import sys, numpy as np, pandas as pd, time
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); os.chdir(f'{WORK}/map')
import onestep_maps as om
d = pd.read_parquet(f'{WORK}/exp/cond_daily.parquet'); S = pd.read_csv(f'{WORK}/exp/cond_station_summary.csv', index_col=0, parse_dates=['first', 'last'])
new = S[S.inside & (S.n >= 24) & (S.km_to_model_station > 0.3)]
d = d[d.station.isin(new.index) & (d.date >= '2000-04-01') & (d.date <= '2025-12-31')].copy()
H, W = om.H, om.W; r = np.floor((om.TOP - new.lat.values) / om.D30).astype(int); c = np.floor((new.lon.values - om.LEFT) / om.D30).astype(int)
pos = np.searchsorted(om.IDX, r * W + c); assert (om.IDX[pos] == r * W + c).all()
st = np.ascontiguousarray(np.asarray(om.STA[:, pos]).T); wt = om.weights(np.asarray(om.LL[0][pos]), np.asarray(om.LL[1][pos]), np.asarray(om.UT[0][pos]), np.asarray(om.UT[1][pos]))
smap = dict(zip(new.index, range(len(new)))); rows = []; t0 = time.time()
for y in sorted(d.date.dt.year.unique()):
    E, ed = om.eden_year(y); eidx = pd.Index(ed); dy = d[d.date.dt.year == y]
    for day, g in dy.groupby('date'):
        X = om.features(st, wt, day, E, eidx.get_loc(day), len(pos))
        for _, s in g.iterrows():
            rows.append(dict(date=day, station=s.station, salinity=s.sal, **dict(zip(om.F3, X[smap[s.station]]))))
    print(y, len(rows), round(time.time() - t0), 's', flush=True)
F = pd.DataFrame(rows); F['src'] = 'dbhydro_spc'; F.to_parquet(f'{WORK}/exp/fresh_station_features.parquet'); print(F.shape, F.isna().mean().round(3).sort_values().tail(5))
