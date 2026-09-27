"""Pack per-file outputs of the extraction scripts into the single files the model scripts read.
Run from the working directory (see README): after extract_eden_grid_at_stations.py,
eden_monthly_grids.py and veg_spectra_rows_L30.py."""
import numpy as np, pandas as pd, glob, os
# 1. EDEN grid at stations (one CSV per quarter) -> one table
f=sorted(glob.glob('downloads/20[12][0-9]_q[1-4].csv'))
if f:
    G=pd.concat([pd.read_csv(x) for x in f]); G=G[(G.date>='2013-01-01')&(G.date<='2025-12-31')]
    G.to_csv('downloads/EDEN_grid_stage_at_stations_2013_2025.csv.gz',index=False); print('EDEN station table',G.shape)
# 2. EDEN 30/90-day mean grids at mid-month -> one npz (keys YYYY-MM-15_m30 / _m90)
f=sorted(glob.glob('downloads/eden_monthly/*.npz'))
if f:
    out={}
    for x in f:
        d=np.load(x); k=os.path.basename(x)[:-4]
        for w in ('m30','m90'): out[f'{k}_{w}']=d[w]
    np.savez_compressed('downloads/eden_monthly_mean_grids_2013_2025.npz',**out); print('EDEN monthly grids',len(f))
# 3. HLS L30 vegetation rows -> veg_rows_L30.npz (keys row_RRRR / cnt_RRRR)
f=sorted(glob.glob('veg_rows/row_*.npz'))
if f:
    out={}
    for x in f:
        r=os.path.basename(x)[4:8]; z=np.load(x); out[f'row_{r}']=z['M']; out[f'cnt_{r}']=z['C']
    np.savez_compressed('veg_rows_L30.npz',**out); print('veg rows',len(f))
