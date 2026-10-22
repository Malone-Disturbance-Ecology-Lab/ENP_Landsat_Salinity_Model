"""Masked daily 30 m salinity, final one-step model, April 2000 - December 2025.
For each month: salinity from daily30m_onestep_final/ is set to no data (-1) on days when the surface-water mask
(wetdry30m/) is 0 (dry: EDEN stage below the lower DEM) or 9 (never wet: upland or built). Uncertain days (code 2) are kept,
as in the paper. Output: daily30m_masked_final/YYYY/ENP_salinity30m_masked_daily_YYYY_MM.tif, one int16 band per day
(psu x 10, nodata -1), deflate + predictor 2, 256 x 256 tiles; band descriptions are the dates.
Usage: python mask_daily30.py <root> <month_index...>  (0 = 2000-04 ... 308 = 2025-12)"""
import sys, os, numpy as np, pandas as pd, rasterio
from rasterio.windows import Window
R = sys.argv[1]; MON = pd.date_range('2000-04-01', '2025-12-01', freq='MS')
for m in map(int, sys.argv[2:]):
    M0 = MON[m]; tag = f'{M0.year}_{M0.month:02d}'; od = f'{R}/daily30m_masked_final/{M0.year}'; os.makedirs(od, exist_ok=True)
    out = f'{od}/ENP_salinity30m_masked_daily_{tag}.tif'
    if os.path.exists(out): print('skip', tag); continue
    with rasterio.open(f'{R}/daily30m_onestep_final/{M0.year}/ENP_salinity30m_onestep_daily_{tag}.tif') as s, \
         rasterio.open(f'{R}/wetdry30m/{M0.year}/ENP_wetdry30m_daily_{tag}.tif') as c:
        assert s.count == c.count == M0.days_in_month and s.shape == c.shape
        prof = s.profile.copy(); prof.update(dtype='int16', nodata=-1, compress='deflate', predictor=2, tiled=True, blockxsize=256, blockysize=256, BIGTIFF='IF_SAFER')
        tmp = out + f'.{os.getpid()}.part'; days = pd.date_range(M0, periods=s.count)
        with rasterio.open(tmp, 'w', **prof) as o:
            for r0 in range(0, s.height, 512):
                w = Window(0, r0, s.width, min(512, s.height - r0)); S = s.read(window=w); C = c.read(window=w)
                S[(C == 0) | (C == 9)] = -1; o.write(S, window=w)
            for i, d in enumerate(days, 1): o.set_band_description(i, d.strftime('%Y-%m-%d'))
            o.update_tags(product='Daily 30 m surface-water salinity, Everglades National Park (final model, masked)', units='psu x 10',
                          scale='0.1', nodata='-1 = outside domain, or no surface water that day (mask code 0 dry or 9 never wet)',
                          mask=f'wetdry30m/{M0.year}/ENP_wetdry30m_daily_{tag}.tif; codes kept 1,2,3,4,5 (2 = uncertain, see companion mask)',
                          source=f'daily30m_onestep_final/{M0.year}/ENP_salinity30m_onestep_daily_{tag}.tif')
    os.replace(tmp, out); print('done', tag, flush=True)
