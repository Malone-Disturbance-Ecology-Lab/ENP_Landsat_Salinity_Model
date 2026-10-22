"""Full-model DAILY salinity at 30 m from the final ONE-STEP (no-spectra) model, for HPC array jobs.
Model: released one-step predictors + sea level (msl_1m, msl_3m, msl_12m), trained 2000-04..2025-12 at all non-anchor stations
plus 13 fresh DBHydro conductivity sites (final_model_onestep.pkl / m3_final.so, 2026-10-13).
Month index 0 = 2000-04 ... 308 = 2025-12.
Usage:  python daily30_onestep_hpc.py <month_index> <input_dir> <output_dir>
Writes <output_dir>/YYYY/ENP_salinity30m_onestep_daily_YYYY_MM.tif (one band per day; int16 psu*10; -1 = no data; EPSG:4326).

Same model, inputs and feature construction as best_two_step/hindcast_onestep (onestep_maps.py), run at every 30 m cell
for every day instead of only on the 15th. The monthly 30 m hindcast maps are bilinear(250 m daily) + the time-interpolated
15th-of-month detail layer; this job predicts every day directly.

<input_dir> needs (from setup_bouchet.slurm plus the hindcast extras listed in README_HPC.md):
  final_model_onestep.pkl, m3_final.so (optional, ~12x faster), msl_by_date.csv, static30.npy, static30_keys.json, px30_index.npy, lonlat30.npy, utm30.npy,
  gridmet_filled_1998_2013.npz (before 2013), gridmet_filled.npz (2013 on), daily_drivers_1998_2025.csv, eden_anyv.npy
EDEN v3 quarterly NetCDFs: $EDEN_DIR (default <input_dir>/../data/EDEN_GRIDS_Q), named YYYY_qN.nc ($EDEN_PREFIX prepended if set).
fillers.py must sit next to this script."""
import sys, os, json, time, pickle, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp, rasterio
from scipy import ndimage as ndi
from scipy.signal import fftconvolve
from rasterio.transform import Affine
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fillers import laplace_solver, nearest_inside
m_job = int(sys.argv[1]); IN = sys.argv[2]; OUTD = sys.argv[3]; t0 = time.time()
EDIR = os.environ.get('EDEN_DIR', f'{IN}/../data/EDEN_GRIDS_Q'); EPRE = os.environ.get('EDEN_PREFIX', '')
H, W = 3839, 4171; D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873
MONTHS = pd.date_range('2000-04-01', '2025-12-01', freq='MS'); M0 = MONTHS[m_job]
days = pd.date_range(M0, M0 + pd.offsets.MonthEnd(0))
FM = pickle.load(open(f'{IN}/final_model_onestep.pkl', 'rb')); F3 = list(FM['m3'].feature_names_in_); assert F3 == FM['F'] and 'msl_12m' in F3
try:
    import tl2cgen; PRED = tl2cgen.Predictor(f'{IN}/m3_final.so', nthread=int(os.environ.get('NTHREAD', '4')))
    def pred(X): return PRED.predict(tl2cgen.DMatrix(X)).ravel()
except Exception as e:
    print('tl2cgen unavailable, scikit-learn predict:', e, flush=True); pred = lambda X: FM['m3'].predict(X)
STK = json.load(open(f'{IN}/static30_keys.json')); STA = np.load(f'{IN}/static30.npy', mmap_mode='r')
IDX = np.load(f'{IN}/px30_index.npy'); LL = np.load(f'{IN}/lonlat30.npy', mmap_mode='r'); UT = np.load(f'{IN}/utm30.npy', mmap_mode='r'); NP = len(IDX)
DRV = pd.read_csv(f'{IN}/daily_drivers_1998_2025.csv', parse_dates=['date']).set_index('date')
DRV = DRV.join(pd.read_csv(f'{IN}/msl_by_date.csv', parse_dates=['date']).set_index('date'))   # sea level, NOAA monthly MSL anomaly (m)
DCOLS = ['SRS90', 'TSC11190'] + [c for c in DRV.columns if c.startswith('anc_')] + ['msl_1m', 'msl_3m', 'msl_12m']
CLIM = ['P30', 'P90', 'P180', 'P365', 'R30', 'T90']
G = np.load(f'{IN}/gridmet_filled_1998_2013.npz' if M0 < pd.Timestamp('2013-01-01') else f'{IN}/gridmet_filled.npz')
glat_s, glon = G['glat_s'], G['glon']; GD = pd.DatetimeIndex(pd.to_datetime(G['dates']))
gi = [GD.get_loc(d) for d in days]; GA = {k: G[k][gi] for k in CLIM}               # only this month's days
# ---- EDEN geometry and harmonic fill (as onestep_maps.py / eden_daily.py)
anyv = np.load(f'{IN}/eden_anyv.npy')
r_ = 12; yy, xx = np.mgrid[-r_:r_ + 1, -r_:r_ + 1]; KK = ((yy ** 2 + xx ** 2) <= 12.5 ** 2).astype('float32')
reg5 = ndi.convolve(anyv.astype('float32'), KK, mode='constant') > 0; den5 = ndi.convolve(anyv.astype('float32'), KK, mode='constant')
PS, PW, PE, PN = 140, 70, 10, 10; y0p = 2790200 - PS * 400; x0p = 463400 - PW * 400
X30, Y30 = UT[0], UT[1]
er0 = max(int(np.floor((np.min(Y30) - y0p) / 400)) - 2, 0); er1 = int(np.ceil((np.max(Y30) - y0p) / 400)) + 3
ec0 = max(int(np.floor((np.min(X30) - x0p) / 400)) - 2, 0); ec1 = int(np.ceil((np.max(X30) - x0p) / 400)) + 3
FILL = {'c': laplace_solver(np.pad(anyv, ((PS, PN), (PW, PE)))), '5': laplace_solver(np.pad(reg5, ((PS, PN), (PW, PE))))}
KNOWN = {'c': np.pad(anyv, ((PS, PN), (PW, PE))), '5': np.pad(reg5, ((PS, PN), (PW, PE)))}
EK = {'stage_cell_cm30': 'c', 'stage_5km_cm30': '5', 'stage_5km_cm90': '5'}

def eden_month():
    """Filled, cropped trailing-window EDEN fields for every day of this month: dict k -> (ndays, ch, cw).
    Reads the quarters covering the 90 days before the month through its end."""
    qs = pd.period_range((M0 - pd.Timedelta(days=100)).to_period('Q'), days[-1].to_period('Q'), freq='Q')
    files = [f'{EDIR}/{EPRE}{q.year}_q{q.quarter}.nc' for q in qs]; files = [f for f in files if os.path.exists(f)]
    C, F5, dates = [], [], []
    for f in files:
        nc = sio.netcdf_file(f, 'r', mmap=False); s = nc.variables['stage'][:].astype('float32'); s[(s < -900) | (s > 9000)] = np.nan
        u = nc.variables['time'].units; u = u.decode() if isinstance(u, bytes) else u
        dates += list(pd.Timestamp(u.split('since')[1].strip()[:10]) + pd.to_timedelta(nc.variables['time'][:].astype(float), unit='D'))
        for k in range(s.shape[0]): s[k] = nearest_inside(s[k], anyv)
        num = fftconvolve(np.nan_to_num(s), KK[None], mode='same', axes=(1, 2))
        C.append(s); F5.append(np.where(reg5[None], num / np.maximum(den5, 1)[None], np.nan).astype('float32'))
    C = np.concatenate(C); F5 = np.concatenate(F5); dates = pd.DatetimeIndex(dates).normalize()
    assert (np.diff(dates.values).astype('timedelta64[D]').astype(int) == 1).all()
    def trailing(x, w, minp):
        o = np.empty(x.shape, 'float32')
        for r0 in range(0, x.shape[1], 45):
            xb = x[:, r0:r0 + 45]; v = np.isfinite(xb); c = np.cumsum(np.where(v, xb, 0), 0, dtype='float64'); n = np.cumsum(v, 0)
            cs = np.concatenate([np.zeros((1,) + xb.shape[1:]), c]); ns = np.concatenate([np.zeros((1,) + xb.shape[1:]), n])
            i = np.arange(1, len(xb) + 1); j = np.maximum(i - w, 0); S = cs[i] - cs[j]; N = ns[i] - ns[j]
            o[:, r0:r0 + 45] = np.where(N >= minp, S / np.maximum(N, 1), np.nan)
        return o
    NAT = {'stage_cell_cm30': trailing(C, 30, 10), 'stage_5km_cm30': trailing(F5, 30, 10), 'stage_5km_cm90': trailing(F5, 90, 30)}
    keep = [dates.get_loc(d) for d in days]; E = {}
    for k, A in NAT.items():
        g = np.pad(A[keep], ((0, 0), (PS, PN), (PW, PE)), constant_values=np.nan); ok = np.isfinite(g).any(axis=(1, 2))
        res = FILL[EK[k]](np.where(np.isfinite(g), g, 0)).astype('float32'); res[~ok] = np.nan
        res[~np.isfinite(g) & np.isfinite(res) & KNOWN[EK[k]][None]] = np.nan
        E[k] = res[:, er0:er1, ec0:ec1]
    return E

def _interp_matrix(fy, fx, iy, ix, ncol_grid, ncells):
    n = len(fy); rows = np.repeat(np.arange(n), 4)
    cols = np.stack([iy * ncol_grid + ix, iy * ncol_grid + ix + 1, (iy + 1) * ncol_grid + ix, (iy + 1) * ncol_grid + ix + 1], 1).ravel()
    vals = np.stack([(1 - fy) * (1 - fx), (1 - fy) * fx, fy * (1 - fx), fy * fx], 1).ravel()
    return sp.csr_matrix((vals.astype('float32'), (rows, cols)), shape=(n, ncells))
def weights(lon, lat, X, Y):
    giy = np.clip(np.searchsorted(glat_s, lat) - 1, 0, len(glat_s) - 2); gix = np.clip(np.searchsorted(glon, lon) - 1, 0, len(glon) - 2)
    gfy = np.clip((lat - glat_s[giy]) / (glat_s[giy + 1] - glat_s[giy]), 0, 1); gfx = np.clip((lon - glon[gix]) / (glon[gix + 1] - glon[gix]), 0, 1)
    ry = (Y - y0p) / 400 - er0; rx = (X - x0p) / 400 - ec0; eh, ew = er1 - er0, ec1 - ec0
    eiy = np.clip(np.floor(ry).astype(int), 0, eh - 2); eix = np.clip(np.floor(rx).astype(int), 0, ew - 2)
    return dict(g=_interp_matrix(gfy, gfx, giy, gix, len(glon), len(glat_s) * len(glon)),
                e=_interp_matrix(np.clip(ry - eiy, 0, 1), np.clip(rx - eix, 0, 1), eiy, eix, ew, eh * ew))
def sbil(M, A2):
    a = np.asarray(A2, 'float32').ravel(); ok = np.isfinite(a); num = M @ np.where(ok, a, 0); den = M @ ok.astype('float32')
    return np.where(den > 0, num / np.maximum(den, 1e-9), np.nan).astype('float32')
COLI = {c: i for i, c in enumerate(F3)}; SI = [COLI[k] for k in STK]
def features(st, wt, t, E, n):
    d = days[t]; X = np.empty((n, len(F3)), 'float32'); X[:, SI] = st
    for k in CLIM: X[:, COLI[k]] = sbil(wt['g'], GA[k][t])
    for k in EK: X[:, COLI[k]] = sbil(wt['e'], E[k][t])
    r = DRV.loc[d]
    for c in DCOLS: X[:, COLI[c]] = r[c]
    doy = d.dayofyear; X[:, COLI['doy_cos']] = np.cos(2 * np.pi * doy / 365.25); X[:, COLI['doy_sin']] = np.sin(2 * np.pi * doy / 365.25)
    w = st[:, STK.index('w_SRS')] + st[:, STK.index('w_BigCypress')]; v = st[:, STK.index('w_TaylorSlough')] + st[:, STK.index('w_C111')]
    X[:, COLI['src_flow90_w']] = w * X[:, COLI['SRS90']] + v * X[:, COLI['TSC11190']]
    return X

E = eden_month(); print(M0.date(), 'EDEN', round(time.time() - t0), 's', flush=True)
OUT = np.full((len(days), H * W), -1, 'int16')
RR = [int(v) for v in os.environ.get('ROWS', '0,%d' % H).split(',')]            # optional row range for testing
a0, b0 = np.searchsorted(IDX, [RR[0] * W, RR[1] * W]); TILE = int(os.environ.get('TILE_CELLS', '1000000'))
for a in range(a0, b0, TILE):
    b = min(a + TILE, b0); n = b - a
    st = np.ascontiguousarray(np.asarray(STA[:, a:b]).T)
    wt = weights(np.asarray(LL[0, a:b]), np.asarray(LL[1, a:b]), np.asarray(UT[0, a:b]), np.asarray(UT[1, a:b]))
    for t in range(len(days)):
        OUT[t, IDX[a:b]] = np.round(np.clip(pred(features(st, wt, t, E, n)), 0, None) * 10).astype('int16')
    print(f'cells {a}-{b} {round(time.time() - t0)}s', flush=True)
os.makedirs(f'{OUTD}/{M0.year}', exist_ok=True)
prof = dict(driver='GTiff', width=W, height=H, count=len(days), dtype='int16', nodata=-1, crs='EPSG:4326',
            transform=Affine(D30, 0, LEFT, 0, -D30, TOP), compress='deflate', predictor=2, tiled=True, blockxsize=256, blockysize=256, BIGTIFF='IF_SAFER')
with rasterio.open(f'{OUTD}/{M0.year}/ENP_salinity30m_onestep_daily_{M0.year}_{M0.month:02d}.tif', 'w', **prof) as o:
    o.write(OUT.reshape(len(days), H, W)); [o.set_band_description(i + 1, str(d.date())) for i, d in enumerate(days)]
    o.update_tags(units='daily salinity, psu x 10 (divide by 10); -1 = no data',
                  model='final one-step model (no vegetation spectra; sea level), trained 2000-2025 with 13 fresh sites; 30 m full model every day')
print('done', M0.date(), round(time.time() - t0), 's')
