"""Climatology of the final one-step record from the wet-day monthly files (days without surface water excluded), 2000-04..2025-12.
Writes <root>/climatology_final/:
  ENP_salinity30m_climatology_2000_2025.tif  (int16 psu x 10 unless noted; -1 = no data)
     1 mean over all wet days, 2000-04..2025-12 (day-weighted)
     2 dry season (March-May) mean            3 wet season (August-October) mean
     4 May mean (end of dry season)           5 October mean (end of wet season)
     6 interannual SD of annual wet-day means, 2001-2025 (psu x 10)
     7 share of days with surface water (percent x 10)
     8 95th percentile of monthly wet-day means (psu x 10)
  ENP_salinity30m_annual_mean_2000_2025.tif  (26 bands; annual wet-day means; 2000 = Apr-Dec)
Usage: python climatology_final.py <root>"""
import sys, os, numpy as np, pandas as pd, rasterio
R = sys.argv[1]; I = f'{R}/monthly30m_wetdays_onestep_final'; O = f'{R}/climatology_final'; os.makedirs(O, exist_ok=True)
MONTHS = pd.date_range('2000-04-01', '2025-12-01', freq='MS'); YRS = list(range(2000, 2026))
with rasterio.open(f'{I}/ENP_salinity30m_monthly_wetdays_2000_04.tif') as s: prof = s.profile; H, W = s.height, s.width
N = H * W
acc = {k: (np.zeros(N), np.zeros(N)) for k in ('all', 'dry', 'wet', 'may', 'oct')}
yr = {y: (np.zeros(N, 'float32'), np.zeros(N, 'float32')) for y in YRS}
tot_days = np.zeros(N, 'float32'); wet_days = np.zeros(N, 'float32')
mon = None; DOM = None
for i, M0 in enumerate(MONTHS):
    with rasterio.open(f'{I}/ENP_salinity30m_monthly_wetdays_{M0.year}_{M0.month:02d}.tif') as s: b = s.read()
    v = b[0].ravel(); n = b[1].ravel().astype('float32'); dom = b[1].ravel() >= 0; ok = v >= 0
    if mon is None: DOM = np.flatnonzero(dom); mon = np.full((len(MONTHS), len(DOM)), np.nan, 'float32')
    x = np.where(ok, v / 10.0, 0); w = np.where(ok, n, 0); mon[i] = np.where(ok, v / 10.0, np.nan)[DOM]
    tot_days += np.where(dom, M0.days_in_month, 0); wet_days += w
    keys = ['all'] + (['dry'] if M0.month in (3, 4, 5) else []) + (['wet'] if M0.month in (8, 9, 10) else []) + (['may'] if M0.month == 5 else []) + (['oct'] if M0.month == 10 else [])
    for k in keys: acc[k][0][:] += x * w; acc[k][1][:] += w
    yr[M0.year][0][:] += (x * w).astype('float32'); yr[M0.year][1][:] += w
    if i % 24 == 0: print('read', M0.date(), flush=True)
def fin(s, w): return np.where(w > 0, s / np.maximum(w, 1), np.nan)
A = {k: fin(*v) for k, v in acc.items()}
Y = np.stack([fin(*yr[y]) for y in YRS])
with np.errstate(all='ignore'):
    sd = np.nanstd(Y[1:], 0, ddof=1); sd[np.isfinite(Y[1:]).sum(0) < 20] = np.nan
    p95 = np.full(N, np.nan, 'float32')
    for c0 in range(0, len(DOM), 400000): p95[DOM[c0:c0 + 400000]] = np.nanpercentile(mon[:, c0:c0 + 400000], 95, axis=0)
pw = np.where(tot_days > 0, 100 * wet_days / np.maximum(tot_days, 1), np.nan)
def i16(a): return np.where(np.isfinite(a), np.round(a * 10), -1).astype('int16').reshape(H, W)
bands = [('mean salinity, all wet days 2000-04..2025-12', A['all']), ('dry season (Mar-May) mean', A['dry']), ('wet season (Aug-Oct) mean', A['wet']),
         ('May mean', A['may']), ('October mean', A['oct']), ('interannual SD of annual means 2001-2025', sd), ('percent of days with surface water (x10)', pw),
         ('95th percentile of monthly wet-day means', p95)]
p = dict(prof); p.update(count=len(bands), dtype='int16', nodata=-1)
with rasterio.open(f'{O}/ENP_salinity30m_climatology_2000_2025.tif', 'w', **p) as o:
    for j, (d, a) in enumerate(bands): o.write(i16(a), j + 1); o.set_band_description(j + 1, d)
    o.update_tags(units='psu x 10 (divide by 10); band 7 percent x 10; -1 = no data', model='final one-step model; days without surface water excluded')
p.update(count=len(YRS))
with rasterio.open(f'{O}/ENP_salinity30m_annual_mean_2000_2025.tif', 'w', **p) as o:
    for j, y in enumerate(YRS): o.write(i16(Y[j]), j + 1); o.set_band_description(j + 1, str(y) + (' (Apr-Dec)' if y == 2000 else ''))
    o.update_tags(units='annual mean salinity over wet days, psu x 10; -1 = no data')
print('climatology done')
