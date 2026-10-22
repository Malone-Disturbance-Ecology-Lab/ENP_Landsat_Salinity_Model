"""Wet-day monthly means from the final one-step daily 30 m maps, one month per task (0 = 2000-04 ... 308 = 2025-12).
Writes <root>/monthly30m_wetdays_onestep_final/ENP_salinity30m_monthly_wetdays_YYYY_MM.tif:
band 1 mean salinity over days kept (codes 0 dry and 9 never wet dropped; psu x 10; -1 none), band 2 number of days kept.
Usage: python wetdays_final.py <month_index> <root>"""
import sys, os, numpy as np, pandas as pd, rasterio
from rasterio.transform import Affine
m = int(sys.argv[1]); R = sys.argv[2]
H, W = 3839, 4171; D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873
M0 = pd.date_range('2000-04-01', '2025-12-01', freq='MS')[m]; tag = f'{M0.year}_{M0.month:02d}'
O = f'{R}/monthly30m_wetdays_onestep_final'; os.makedirs(O, exist_ok=True)
with rasterio.open(f'{R}/daily30m_onestep_final/{M0.year}/ENP_salinity30m_onestep_daily_{tag}.tif') as s: S = s.read()
with rasterio.open(f'{R}/wetdry30m/{M0.year}/ENP_wetdry30m_daily_{tag}.tif') as c: C = c.read()
assert C.shape == S.shape
ok = S >= 0; n = ok.sum(0); keep = ok & (C != 0) & (C != 9); nk = keep.sum(0); tk = np.where(keep, S, 0).sum(0, dtype='int64')
mean = np.where(nk > 0, np.round(tk / np.maximum(nk, 1)), -1).astype('int16')
base = dict(driver='GTiff', crs='EPSG:4326', compress='deflate', predictor=2, tiled=True, blockxsize=256, blockysize=256)
with rasterio.open(f'{O}/ENP_salinity30m_monthly_wetdays_{tag}.tif', 'w', width=W, height=H, count=2, dtype='int16', nodata=-1, transform=Affine(D30, 0, LEFT, 0, -D30, TOP), **base) as o:
    o.write(mean, 1); o.write(np.where(n > 0, nk, -1).astype('int16'), 2)
    o.set_band_description(1, 'mean salinity over days with surface water (psu x 10)'); o.set_band_description(2, 'days with surface water')
    o.update_tags(source=f'ENP_salinity30m_onestep_daily_{tag}.tif (final one-step model)', dropped_codes='[0, 9]')
print('done', tag)
