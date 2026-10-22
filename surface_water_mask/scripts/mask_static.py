"""Static 30 m layers for the surface-water mask (HLS/Landsat 30 m grid, EPSG:4326, 3839 x 4171).
- dem30_{best,lower,upper}_cm.npy: EDEN DEMs v2 (Simons et al. 2024/2025, doi:10.5066/P1XTHZ7D), 10 m -> 30 m MINIMUM (lowest 10 m point in the 30 m cell; 'any surface water in the cell'), cm NAVD88; 30 m means kept as dem30_*_cm.npy
- static class (uint8): 3 permanent water | 9 never wet (upland / built) | 10 EDEN-assessed (daily wet/dry) |
  4 tidal/intermittent coastal wetland (outside EDEN) | 5 not assessed (other wetland outside EDEN) | 255 outside domain
- eden_rc30.npy: EDEN 400 m (row, col) of each domain cell, for daily stage lookup"""
import numpy as np, pandas as pd, rasterio, os, sys
from rasterio.warp import reproject, Resampling
from rasterio.transform import Affine
U = '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'; M = 'work/map/in'; O = 'out'
DEMD = f'{U}/EDEN_DEMs_With_Uncertainty_Treatment_v2/EDEN_DEMs_With_Uncertainty_Treatment_v2/10m_DEMs'
H, W = 3839, 4171; D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873
T30 = Affine(D30, 0, LEFT, 0, -D30, TOP)
IDX = np.load(f'{M}/px30_index.npy'); UT = np.load(f'{M}/utm30.npy')
dom = np.zeros(H * W, bool); dom[IDX] = True
for k in ['Best', 'Lower', 'Upper']:
    if os.path.exists(f'{O}/dem30_{k.lower()}_cm.npy'): continue
    with rasterio.open(f'{DEMD}/EDEN_10m_{k}_Elevation_v2.tif') as r:
        a = r.read(1).astype('float32'); nd = r.nodata; a[(a == nd) | (a <= -9000) | (a >= 30000)] = np.nan
        dst = np.full((H, W), np.nan, 'float32')
        reproject(a, dst, src_transform=r.transform, src_crs=r.crs.to_epsg() and 'EPSG:6346' or r.crs, src_nodata=np.nan,
                  dst_transform=T30, dst_crs='EPSG:4326', dst_nodata=np.nan, resampling=Resampling.average)
    np.save(f'{O}/dem30_{k.lower()}_cm.npy', dst); print(k, np.nanmin(dst), np.nanmedian(dst), np.nanmax(dst), int(np.isfinite(dst).sum()), flush=True); del a
B, L, Up = [np.load(f'{O}/dem30min_{k}_cm.npy') for k in ('best', 'lower', 'upper')]
# EDEN 400 m cell of each domain cell (cell centres at 463400 + 400 c, 2790200 + 400 r; rows increase northward)
anyv = np.load(f'{M}/eden_anyv.npy'); EH, EW = anyv.shape
rr = np.round((UT[1] - 2790200) / 400).astype(int); cc = np.round((UT[0] - 463400) / 400).astype(int)
inside = (rr >= 0) & (rr < EH) & (cc >= 0) & (cc < EW); rc = np.where(inside, rr, -1).astype('int16'), np.where(inside, cc, -1).astype('int16')
np.save(f'{O}/eden_rc30.npy', np.stack(rc))
eden_ok = np.zeros(H * W, bool); eden_ok[IDX] = inside & anyv[np.clip(rr, 0, EH - 1), np.clip(cc, 0, EW - 1)]
veg = rasterio.open(f'{U}/site_variables/ENP_VegMap_30m.tif').read(1).ravel()
codes = pd.read_csv(f'{U}/site_variables/vegmap_codes.csv'); cmap = dict(zip(codes.id, codes.cls)); vcode = dict(zip(codes.id, codes.VegCode))
ids = np.arange(veg.max() + 1)
is_water_v = np.isin(veg, [i for i in ids if vcode.get(i) in ('ONW', 'OAC')])
is_never = np.isin(veg, [i for i in ids if cmap.get(i) == 'Upland' or vcode.get(i) in ('OA', 'OAA', 'OAG', 'OAL', 'OAR', 'OAS', 'ONO')])
is_coastal = np.isin(veg, [i for i in ids if cmap.get(i) in ('Mangrove', 'Salt marsh', 'Tidal/salt flat, beach')])
wf = rasterio.open(f'{U}/site_variables/ENP_WaterFrequency_L30.tif').read(1).ravel()
perm = (wf >= 0.5) | is_water_v
demok = np.isfinite(B.ravel()) & np.isfinite(L.ravel()) & np.isfinite(Up.ravel())
S = np.full(H * W, 255, 'uint8'); S[dom] = 5
S[dom & eden_ok & demok] = 10
S[dom & is_coastal] = 4          # mangrove / salt marsh / tidal flat: tidal flooding, not captured by daily-mean EDEN stage
S[dom & is_never] = 9
S[dom & perm] = 3
np.save(f'{O}/static_class30.npy', S.reshape(H, W))
lab = {3: 'permanent water', 9: 'never wet (upland/built)', 10: 'EDEN-assessed daily', 4: 'tidal/intermittent coastal wetland', 5: 'not assessed'}
n = dom.sum()
for k, v in lab.items(): print(f'{v:38s} {100 * (S == k).sum() / n:5.1f}% ({(S == k).sum():,} cells)')
print('EDEN-native cells with DEM:', int((dom & eden_ok & demok).sum()), ' EDEN-native without DEM:', int((dom & eden_ok & ~demok).sum()))
