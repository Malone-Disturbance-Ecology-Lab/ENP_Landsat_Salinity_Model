"""Assemble yearly monthly-250 m files (mean, lower90, upper90) from derive_final.py parts. Usage: python reduce_monthly250.py <root>"""
import sys, json, glob, numpy as np, pandas as pd, rasterio
from rasterio.transform import Affine
R = sys.argv[1]; O = f'{R}/onestep_final_products'; import os; os.makedirs(f'{O}/monthly250', exist_ok=True)
D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873
Q = json.load(open(f'{R}/hpc_inputs/q_levels.json')); MONTHS = pd.date_range('2000-04-01', '2025-12-01', freq='MS')
assert len(glob.glob(f'{O}/parts250/*.npz')) == len(MONTHS)
base = dict(driver='GTiff', crs='EPSG:4326', compress='deflate', predictor=3, tiled=True, blockxsize=256, blockysize=256)
for y in sorted(set(MONTHS.year)):
    ms = MONTHS[MONTHS.year == y]; P = [np.load(f'{O}/parts250/{m.year}_{m.month:02d}.npz') for m in ms]
    for k in ('mean', 'lower90', 'upper90'):
        A = np.stack([p[k] for p in P]).astype('float32')
        with rasterio.open(f'{O}/monthly250/ENP_salinity_onestep_monthly_{k}_{y}.tif', 'w', width=521, height=480, count=len(ms), dtype='float32', nodata=np.nan, transform=Affine(D30 * 8, 0, LEFT, 0, -D30 * 8, TOP), **base) as o:
            o.write(A); [o.set_band_description(i + 1, f'{m.year}-{m.month:02d}') for i, m in enumerate(ms)]
            o.update_tags(units='psu', model='final one-step model (no vegetation spectra; sea level), trained 2000-2025 with 13 fresh sites; 250 m',
                          interval='90% split-conformal half-widths by predicted level ' + str(np.round(Q['q'], 2).tolist()) + ' from ' + Q['source'])
print('monthly250 done')
