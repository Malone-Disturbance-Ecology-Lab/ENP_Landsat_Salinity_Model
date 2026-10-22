"""Score the daily wet/dry rule at EDEN gauges: EDEN 400 m grid stage (bilinear, native cells only) vs EDEN DEM
lower/upper bounds (30 m mean and 10 m point) against the gauge's own D (dry) flag, 2013-2025."""
import numpy as np, pandas as pd, glob, scipy.io as sio, rasterio
from pyproj import Transformer
U = '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'; M = 'work/map'
D10 = f'{U}/EDEN_DEMs_With_Uncertainty_Treatment_v2/EDEN_DEMs_With_Uncertainty_Treatment_v2/10m_DEMs'
raw = pd.concat([pd.read_csv(f) for f in glob.glob(f'{U}/eden_stage/raw_downloads/eden2_*.csv') if 'coords' not in f])
raw = raw.drop_duplicates(['site', 'date']); raw['date'] = pd.to_datetime(raw.date)
inv = pd.read_csv(f'{U}/eden_stage/EDEN_gage_inventory.csv')
G = inv[['site', 'lat', 'lon', 'type', 'pct_dry_flag']].dropna(subset=['lat'])
tr = Transformer.from_crs('EPSG:4326', 'EPSG:26917', always_xy=True); G['x'], G['y'] = tr.transform(G.lon.values, G.lat.values)
anyv = np.load(f'{M}/in/eden_anyv.npy'); EH, EW = anyv.shape
fy = (G.y.values - 2790200) / 400; fx = (G.x.values - 463400) / 400
iy = np.floor(fy).astype(int); ix = np.floor(fx).astype(int); wy = fy - iy; wx = fx - ix
G['in_grid'] = (iy >= 0) & (iy < EH - 1) & (ix >= 0) & (ix < EW - 1)
G['eden_native'] = G.in_grid & anyv[np.clip(np.round(fy).astype(int), 0, EH - 1), np.clip(np.round(fx).astype(int), 0, EW - 1)]
# DEM: 30 m mean (from the static product) and 10 m point
IDXr = None
H, W = 3839, 4171; D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873
r30 = np.floor((TOP - G.lat.values) / D30).astype(int); c30 = np.floor((G.lon.values - LEFT) / D30).astype(int)
ok30 = (r30 >= 0) & (r30 < H) & (c30 >= 0) & (c30 < W)
for k in ('best', 'lower', 'upper'):
    A = np.load(f'out/dem30_{k}_cm.npy'); G[f'dem30_{k}'] = np.where(ok30, A[np.clip(r30, 0, H - 1), np.clip(c30, 0, W - 1)], np.nan)
    with rasterio.open(f'{D10}/EDEN_10m_{k.capitalize()}_Elevation_v2.tif') as r:
        v = np.array([x[0] for x in r.sample(list(zip(G.x, G.y)))], 'float32'); v[(v == r.nodata) | (np.abs(v) > 9000)] = np.nan; G[f'dem10_{k}'] = v
S = np.load('out/static_class30.npy'); G['static'] = np.where(ok30, S[np.clip(r30, 0, H - 1), np.clip(c30, 0, W - 1)], 255)
# EDEN daily grid stage at gauges, 2013-2025
rows = []
for f in sorted(glob.glob(f'{M}/eden/eden_20[12][0-9]_q[1-4].nc')):
    if int(f.split('_')[-2]) < 2013: continue
    nc = sio.netcdf_file(f, 'r', mmap=False); s = nc.variables['stage'][:].astype('float32'); s[(s < -900) | (s > 9000)] = np.nan
    u = nc.variables['time'].units; u = u.decode() if isinstance(u, bytes) else u
    dates = pd.Timestamp(u.split('since')[1].strip()[:10]) + pd.to_timedelta(nc.variables['time'][:].astype(float), unit='D')
    iyc, ixc = np.clip(iy, 0, EH - 2), np.clip(ix, 0, EW - 2)
    v = [s[:, iyc, ixc], s[:, iyc, ixc + 1], s[:, iyc + 1, ixc], s[:, iyc + 1, ixc + 1]]
    w = [(1 - wy) * (1 - wx), (1 - wy) * wx, wy * (1 - wx), wy * wx]
    num = sum(np.where(np.isfinite(a), a * b, 0) for a, b in zip(v, w)); den = sum(np.where(np.isfinite(a), b, 0) for a, b in zip(v, w))
    st = np.where(den > 0, num / np.maximum(den, 1e-9), np.nan)
    rows.append(pd.DataFrame(st, index=pd.DatetimeIndex(dates).normalize(), columns=G.site.values))
ST = pd.concat(rows).stack().rename('grid_cm').reset_index().rename(columns={'level_0': 'date', 'level_1': 'site'})
D = raw.merge(ST, on=['site', 'date']).merge(G, on='site')
D = D[D.flag.isin(['O', 'E', 'D'])].copy(); D['dry_obs'] = D.flag.eq('D'); D['gauge_cm'] = D.stage_ft_NAVD88 * 30.48
def cls(stage, lo, up): return np.where(stage < lo, 'dry', np.where(stage > up, 'wet', 'uncertain'))
D['c_grid30'] = cls(D.grid_cm, D.dem30_lower, D.dem30_upper); D['c_grid10'] = cls(D.grid_cm, D.dem10_lower, D.dem10_upper)
D['c_gauge10'] = cls(D.gauge_cm, D.dem10_lower, D.dem10_upper)
D.to_parquet('out/gauge_wetdry_validation.parquet')
use = D[D.eden_native & (D.type == 'marsh') & np.isfinite(D.dem30_lower)]
print('gauges used', use.site.nunique(), 'gauge-days', len(use), 'dry share', round(use.dry_obs.mean(), 3))
for c in ['c_grid30', 'c_grid10', 'c_gauge10']:
    t = pd.crosstab(use[c], use.dry_obs.map({True: 'flag D', False: 'flag wet'}), normalize='columns').round(3)
    print('\n==', c); print(t.to_string())
    dd = use[use[c] != 'uncertain']; acc = ((dd[c] == 'dry') == dd.dry_obs).mean(); print('accuracy excl. uncertain', round(acc, 3), 'share uncertain', round((use[c] == 'uncertain').mean(), 3))
per = use.groupby('site').apply(lambda g: pd.Series(dict(dry_obs=g.dry_obs.mean(), dry_pred=(g.c_grid30 == 'dry').mean(), unc=(g.c_grid30 == 'uncertain').mean(),
      bias_cm=(g.grid_cm - g.gauge_cm).median(), dem30_minus_dem10=float(g.dem30_best.iloc[0] - g.dem10_best.iloc[0]))))
print(per[per.dry_obs > 0.05].round(2).sort_values('dry_obs', ascending=False).to_string())
per.to_csv('out/gauge_wetdry_by_site.csv')
