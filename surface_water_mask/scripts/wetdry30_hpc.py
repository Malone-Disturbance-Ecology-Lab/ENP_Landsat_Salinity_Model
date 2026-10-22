"""Daily surface-water (wet/dry) mask on the 30 m salinity grid, plus wet-day monthly salinity, one month per task.
Usage: python wetdry30_hpc.py <month_index> <mask_input_dir> <root>
  month_index 0 = 2000-04 ... 308 = 2025-12.
  mask_input_dir: static_class30.npy, dem30_lower_cm_int16.npy, dem30_upper_cm_int16.npy, eden_anyv.npy, lonlat30.npy/utm30.npy/px30_index.npy
  root: has data/EDEN_GRIDS_Q/YYYY_qN.nc, daily30m/ (2013-2025) and daily30m_onestep/ (2000-2012)
Writes <root>/wetdry30m/YYYY/ENP_wetdry30m_daily_YYYY_MM.tif (uint8, one band per day) and
       <root>/monthly30m_wetdays/ENP_salinity30m_monthly_wetdays_YYYY_MM.tif (band 1: mean salinity over days kept, psu x 10, -1 none;
       band 2: number of days kept).
Codes: 0 dry | 1 wet | 2 uncertain (EDEN water level within the DEM lower-upper band) | 3 permanent water |
       4 tidal/intermittent coastal wetland (not assessed daily) | 5 not assessed | 9 never wet (upland/built) | 255 outside.
Daily rule for EDEN cells: dry if stage < lower DEM, wet if stage > upper DEM, else uncertain. EDEN stage is bilinear
on native EDEN cells only; DEMs are the EDEN v2 10 m DEMs (Simons et al. 2025, doi:10.5066/P1XTHZ7D), lowest 10 m point in each 30 m cell (any surface water in the cell).
Salinity is dropped (NA) on codes 0 and 9; uncertain days are kept (set DROP_UNCERTAIN=1 to drop them too)."""
import sys, os, numpy as np, pandas as pd, scipy.io as sio, rasterio
from rasterio.transform import Affine
m_job = int(sys.argv[1]); IN = sys.argv[2]; R = sys.argv[3]
H, W = 3839, 4171; D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873
MONTHS = pd.date_range('2000-04-01', '2025-12-01', freq='MS'); M0 = MONTHS[m_job]; days = pd.date_range(M0, M0 + pd.offsets.MonthEnd(0))
S = np.load(f'{IN}/static_class30.npy').ravel(); A = np.flatnonzero(S == 10)
LO = np.load(f'{IN}/dem30_lower_cm_int16.npy').ravel()[A].astype('float32'); UPv = np.load(f'{IN}/dem30_upper_cm_int16.npy').ravel()[A].astype('float32')
IDX = np.load(f'{IN}/px30_index.npy'); UT = np.load(f'{IN}/utm30.npy', mmap_mode='r'); pos = np.searchsorted(IDX, A)
X, Y = np.asarray(UT[0][pos]), np.asarray(UT[1][pos]); anyv = np.load(f'{IN}/eden_anyv.npy'); EH, EW = anyv.shape
fy = (Y - 2790200) / 400; fx = (X - 463400) / 400; iy = np.clip(np.floor(fy).astype(int), 0, EH - 2); ix = np.clip(np.floor(fx).astype(int), 0, EW - 2)
wy = np.clip(fy - iy, 0, 1); wx = np.clip(fx - ix, 0, 1)
def stage_days():
    out = {}
    for q in sorted({(d.year, (d.month - 1) // 3 + 1) for d in days}):
        nc = sio.netcdf_file(f'{R}/data/EDEN_GRIDS_Q/{q[0]}_q{q[1]}.nc', 'r', mmap=False); s = nc.variables['stage'][:].astype('float32')
        s[(s < -900) | (s > 9000)] = np.nan; s[:, ~anyv] = np.nan
        u = nc.variables['time'].units; u = u.decode() if isinstance(u, bytes) else u
        dts = pd.DatetimeIndex(pd.Timestamp(u.split('since')[1].strip()[:10]) + pd.to_timedelta(nc.variables['time'][:].astype(float), unit='D')).normalize()
        for t, d in enumerate(dts):
            if d in days: out[d] = s[t]
    return out
SD = stage_days(); assert len(SD) == len(days), (len(SD), len(days))
CODE = np.empty((len(days), H * W), 'uint8')
for t, d in enumerate(days):
    g = SD[d]; v = [g[iy, ix], g[iy, ix + 1], g[iy + 1, ix], g[iy + 1, ix + 1]]; w = [(1 - wy) * (1 - wx), (1 - wy) * wx, wy * (1 - wx), wy * wx]
    num = sum(np.where(np.isfinite(a), a * b, 0) for a, b in zip(v, w)); den = sum(np.where(np.isfinite(a), b, 0) for a, b in zip(v, w))
    st = np.where(den > 0, num / np.maximum(den, 1e-9), np.nan)
    c = np.where(st < LO, 0, np.where(st > UPv, 1, 2)).astype('uint8'); c[~np.isfinite(st)] = 5
    row = S.copy(); row[A] = c; CODE[t] = row
prof = dict(driver='GTiff', width=W, height=H, crs='EPSG:4326', transform=Affine(D30, 0, LEFT, 0, -D30, TOP), compress='deflate', predictor=1, tiled=True, blockxsize=256, blockysize=256)
os.makedirs(f'{R}/wetdry30m/{M0.year}', exist_ok=True)
with rasterio.open(f'{R}/wetdry30m/{M0.year}/ENP_wetdry30m_daily_{M0.year}_{M0.month:02d}.tif', 'w', count=len(days), dtype='uint8', nodata=255, **prof) as o:
    o.write(CODE.reshape(len(days), H, W)); [o.set_band_description(i + 1, str(d.date())) for i, d in enumerate(days)]
    o.update_tags(codes='0 dry; 1 wet; 2 uncertain; 3 permanent water; 4 tidal/intermittent coastal wetland; 5 not assessed; 9 never wet (upland/built); 255 outside',
                  rule='EDEN cells: dry if EDEN stage < EDEN v2 lower DEM, wet if > upper DEM (lowest 10 m point of the EDEN v2 10 m DEMs in each 30 m cell)')
# wet-day monthly salinity
f = (f'{R}/daily30m/{M0.year}/ENP_salinity30m_daily_{M0.year}_{M0.month:02d}.tif' if M0 >= pd.Timestamp('2013-01-01')
     else f'{R}/daily30m_onestep/{M0.year}/ENP_salinity30m_onestep_daily_{M0.year}_{M0.month:02d}.tif')
if not os.path.exists(f): print(M0.date(), 'no salinity file; mask only'); sys.exit(0)
drop = [0, 9] + ([2] if os.environ.get('DROP_UNCERTAIN') == '1' else [])
with rasterio.open(f) as r: SAL = r.read().reshape(len(days), -1)
keep = (SAL >= 0) & ~np.isin(CODE, drop)
n = keep.sum(0); mean = np.where(n > 0, np.where(keep, SAL, 0).sum(0, dtype='int64') / np.maximum(n, 1), -1)
os.makedirs(f'{R}/monthly30m_wetdays', exist_ok=True)
with rasterio.open(f'{R}/monthly30m_wetdays/ENP_salinity30m_monthly_wetdays_{M0.year}_{M0.month:02d}.tif', 'w', count=2, dtype='int16', nodata=-1, **{**prof, 'predictor': 2}) as o:
    o.write(np.round(mean).astype('int16').reshape(H, W), 1); o.write(np.where(SAL[0] >= 0, n, -1).astype('int16').reshape(H, W), 2)
    o.set_band_description(1, 'mean salinity over days with surface water (psu x 10)'); o.set_band_description(2, 'days with surface water')
    o.update_tags(source=os.path.basename(f), dropped_codes=str(drop))
dry = (CODE[:, A] == 0).mean(); unc = (CODE[:, A] == 2).mean()
print(M0.date(), f'EDEN cells dry {dry:.3f} uncertain {unc:.3f}; cells with all days dropped {int(((n == 0) & (SAL[0] >= 0)).sum())}')
