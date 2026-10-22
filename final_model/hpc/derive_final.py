"""Products derived from the final one-step daily 30 m maps, one month per task (0 = 2000-04 ... 308 = 2025-12):
  onestep_final_products/monthly30/ENP_salinity30m_onestep_monthly_YYYY_MM.tif   monthly mean of all days (int16 psu x 10; -1 no data)
  onestep_final_products/daily250/YYYY/ENP_salinity_onestep_daily_YYYY_MM.tif    daily values at the 250 m nodes (30 m cells 8i+4, 8j+4)
  onestep_final_products/parts250/YYYY_MM.npz                                     monthly 250 m mean / lower90 / upper90 (assembled by reduce)
  monthly30m_wetdays_onestep_final/ENP_salinity30m_monthly_wetdays_YYYY_MM.tif   months before 2013 only: band 1 mean over days kept
      (codes 0 dry and 9 never wet dropped; psu x 10; -1 none), band 2 number of days kept
Usage: python derive_final.py <month_index> <root>   (root = ~/enp_salinity)"""
import sys, os, json, warnings, numpy as np, pandas as pd, rasterio; warnings.filterwarnings("ignore")
from rasterio.transform import Affine
m = int(sys.argv[1]); R = sys.argv[2]
H, W = 3839, 4171; H2, W2 = 480, 521; D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873
MONTHS = pd.date_range('2000-04-01', '2025-12-01', freq='MS'); M0 = MONTHS[m]; tag = f'{M0.year}_{M0.month:02d}'
O = f'{R}/onestep_final_products'
for d in ('monthly30', f'daily250/{M0.year}', 'parts250'): os.makedirs(f'{O}/{d}', exist_ok=True)
MODEL = 'final one-step model (no vegetation spectra; sea level), trained 2000-2025 with 13 fresh sites (2026-10-13)'
Q = json.load(open(f'{R}/hpc_inputs/q_levels.json')); bins = np.array(Q['bins']); q = np.array(Q['q'])
base = dict(driver='GTiff', crs='EPSG:4326', compress='deflate', predictor=2, tiled=True, blockxsize=256, blockysize=256)
with rasterio.open(f'{R}/daily30m_onestep_final/{M0.year}/ENP_salinity30m_onestep_daily_{tag}.tif') as s:
    S = s.read(); desc = s.descriptions
nd = S.shape[0]
# monthly 30 m mean over all days
ok = S >= 0; n = ok.sum(0); tot = np.where(ok, S, 0).sum(0, dtype='int64')
mon = np.where(n > 0, np.round(tot / np.maximum(n, 1)), -1).astype('int16')
with rasterio.open(f'{O}/monthly30/ENP_salinity30m_onestep_monthly_{tag}.tif', 'w', width=W, height=H, count=1, dtype='int16', nodata=-1, transform=Affine(D30, 0, LEFT, 0, -D30, TOP), **base) as o:
    o.write(mon, 1); o.set_band_description(1, f'{M0.year}-{M0.month:02d}')
    o.update_tags(units='monthly mean salinity, psu x 10 (divide by 10); -1 = no data', grid='HLS/Landsat 30 m (EPSG:4326, 0.000270923 deg)', model=MODEL + '; mean of the daily 30 m maps')
# daily 250 m nodes
D2 = np.ascontiguousarray(S[:, 4::8, 4::8][:, :H2, :W2])
with rasterio.open(f'{O}/daily250/{M0.year}/ENP_salinity_onestep_daily_{tag}.tif', 'w', width=W2, height=H2, count=nd, dtype='int16', nodata=-1, transform=Affine(D30 * 8, 0, LEFT, 0, -D30 * 8, TOP), **base) as o:
    o.write(D2); [o.set_band_description(i + 1, desc[i]) for i in range(nd)]
    o.update_tags(units='salinity psu x 10 (divide by 10); -1 = no data', model=MODEL + '; 250 m nodes of the daily 30 m maps')
p = D2.astype('float32') / 10; p[D2 < 0] = np.nan; hw = q[np.clip(np.digitize(np.nan_to_num(p), bins[1:-1]), 0, 3)]
with np.errstate(all='ignore'):
    np.savez(f'{O}/parts250/{tag}.npz', mean=np.nanmean(p, 0), lower90=np.nanmean(np.clip(p - hw, 0, None), 0), upper90=np.nanmean(p + hw, 0))
# wet-day monthly means before 2013 (from 2013 the primary HLS product is used, unchanged)
if M0 < pd.Timestamp('2013-01-01'):
    os.makedirs(f'{R}/monthly30m_wetdays_onestep_final', exist_ok=True)
    with rasterio.open(f'{R}/wetdry30m/{M0.year}/ENP_wetdry30m_daily_{tag}.tif') as c: C = c.read()
    assert C.shape == S.shape
    keep = ok & (C != 0) & (C != 9); nk = keep.sum(0); tk = np.where(keep, S, 0).sum(0, dtype='int64')
    mean = np.where(nk > 0, np.round(tk / np.maximum(nk, 1)), -1).astype('int16')
    with rasterio.open(f'{R}/monthly30m_wetdays_onestep_final/ENP_salinity30m_monthly_wetdays_{tag}.tif', 'w', width=W, height=H, count=2, dtype='int16', nodata=-1, transform=Affine(D30, 0, LEFT, 0, -D30, TOP), **base) as o:
        o.write(mean, 1); o.write(np.where(n > 0, nk, -1).astype('int16'), 2)
        o.set_band_description(1, 'mean salinity over days with surface water (psu x 10)'); o.set_band_description(2, 'days with surface water')
        o.update_tags(source=f'ENP_salinity30m_onestep_daily_{tag}.tif (final one-step model)', dropped_codes='[0, 9]')
print('done', tag, 'domain cells', int((n > 0).sum()), 'mean psu', round(float(tot.sum() / n.sum() / 10), 2))
