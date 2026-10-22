"""Daily one-step-model (v3, no spectra) features at the 90 station 30 m cells, 2000-01-01 .. 2025-12-31.

Reproduces the mapping pipeline exactly (eden_daily.py + daily30_hpc.py), but samples only the station cells:
  static    : static30 (27 layers) at the station's 30 m cell
  gridMET   : filled P30/P90/P180/P365/R30/T90, bilinear at the cell centre (1998-2013 rebuilt; 2012-2025 original)
  EDEN      : daily stage -> nearest_inside -> 5 km disc mean -> trailing 30/90-day windows -> harmonic extrapolation
              (same anyv mask from the 2013-2014 files, same padding) -> NaN-aware bilinear at the cell centre (UTM 17N)
  drivers   : SRS90, TSC11190 (90-day trailing means, >= 72 days), anchors (30-day >= 10 days, 3-day), doy, src_flow90_w
Steps are cached so the script can be re-run in pieces:  python3 station_features_hindcast.py <step> ...
  step static | eden <year> | assemble"""
import sys, os, glob, gzip, json, numpy as np, pandas as pd
SERVER = os.environ.get('ENP_DATA_ROOT', os.path.expanduser('~/mnt/MaloneLab--Research--ENP_Landsat_Salinity_Model'))
WORK = os.environ.get('HC_WORK', os.path.expanduser('~/work'))
HC = f'{SERVER}/dbhydro_hindcast'; HPC = f'{SERVER}/best_two_step/salinity30m/hpc_inputs'
CACHE = f'{WORK}/hc_cache'; os.makedirs(CACHE, exist_ok=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
H30, W30 = 3839, 4171; D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873


def stations():
    s = pd.read_csv(f'{SERVER}/site_variables/station_site_variables.csv')[['station', 'longitude', 'latitude', 'grab']]
    s['row'] = np.floor((TOP - s.latitude) / D30).astype(int); s['col'] = np.floor((s.longitude - LEFT) / D30).astype(int)
    s['lon30'] = LEFT + (s.col + 0.5) * D30; s['lat30'] = TOP - (s.row + 0.5) * D30            # cell centre, as in the maps
    from pyproj import Transformer
    s['X'], s['Y'] = Transformer.from_crs(4326, 26917, always_xy=True).transform(s.lon30.values, s.lat30.values)
    return s


def step_static():
    import rasterio
    s = stations(); wf = rasterio.open(f'{SERVER}/site_variables/ENP_WaterFrequency_L30.tif').read(1)
    idx = np.flatnonzero((wf > -1).ravel()); flat = s.row.values * W30 + s.col.values
    pos = np.searchsorted(idx, flat); inside = (pos < len(idx)) & (idx[np.minimum(pos, len(idx) - 1)] == flat)
    print('stations inside the 30 m domain:', int(inside.sum()), 'of', len(s))
    keys = json.load(open(f'{HPC}/static30_keys.json')); NP = len(idx)
    parts = sorted(glob.glob(f'{HPC}/static30.npy.gz.part*'))
    class Cat:                                                   # concatenated parts as one stream
        def __init__(s_, fs): s_.fs = iter(fs); s_.f = open(next(s_.fs), 'rb')
        def read(s_, n=-1):
            b = s_.f.read(n)
            while (n < 0 or len(b) < n):
                try: s_.f = open(next(s_.fs), 'rb')
                except StopIteration: break
                b += s_.f.read(-1 if n < 0 else n - len(b))
            return b
    g = gzip.GzipFile(fileobj=Cat(parts)); ver = np.lib.format.read_magic(g); shp, fo, dt = np.lib.format._read_array_header(g, ver)
    assert shp == (len(keys), NP) and not fo, (shp, fo)
    want = np.where(inside, pos, -1); out = np.full((len(keys), len(s)), np.nan, 'float32'); chunk = 4_000_000
    for k in range(len(keys)):
        done = 0
        while done < NP:
            n = min(chunk, NP - done); a = np.frombuffer(g.read(n * 4), dtype=dt, count=n)
            m = (want >= done) & (want < done + n); out[k, m] = a[want[m] - done]; done += n
    T = pd.DataFrame(out.T, columns=keys); T.insert(0, 'station', s.station.values)
    T.to_csv(f'{CACHE}/static_at_stations.csv', index=False); print(T.describe().T[['count', 'mean']].round(2).head(30))


# ---------------------------------------------------------------- EDEN (exactly as eden_daily.py)
def eden_setup():
    import scipy.io as sio
    from scipy import ndimage as ndi
    from fillers import laplace_solver
    fn = f'{CACHE}/eden_anyv.npy'
    if os.path.exists(fn): anyv = np.load(fn)
    else:
        files = sorted(glob.glob(f'{SERVER}/data/EDEN_GRIDS_Q/20[12][0-9]_q[1-4].nc')); files = [f for f in files if 2013 <= int(os.path.basename(f)[:4]) <= 2025]
        anyv = np.zeros((405, 287), bool)
        for f in files[:8]:
            s_ = sio.netcdf_file(f, 'r', mmap=False).variables['stage'][:].astype('float32'); s_[(s_ < -900) | (s_ > 9000)] = np.nan; anyv |= np.isfinite(s_).any(0)
        np.save(fn, anyv)
    r_ = 12; yy, xx = np.mgrid[-r_:r_ + 1, -r_:r_ + 1]; kk = ((yy ** 2 + xx ** 2) <= 12.5 ** 2).astype('float32')
    reg5 = ndi.convolve(anyv.astype('float32'), kk, mode='constant') > 0; den5 = ndi.convolve(anyv.astype('float32'), kk, mode='constant')
    PS, PW, PE, PN = 140, 70, 10, 10
    return anyv, kk, reg5, den5, None, None, (PS, PN, PW, PE)   # fills are applied through eden_station_weights.npz


def eden_files(y0, y1):
    """Quarterly files: 2000-2012 downloaded v3 (~/work/eden/YYYY_qN_v3.nc), 2013-2025 the server copies."""
    out = []
    for y in range(y0, y1 + 1):
        for q in range(1, 5):
            f = f'{WORK}/eden/{y}_q{q}_v3.nc' if y <= 2012 else f'{SERVER}/data/EDEN_GRIDS_Q/{y}_q{q}.nc'
            out.append(f)
    return out


def step_eden(year):
    """Native daily cell and 5 km fields for one year (+90 days before for the windows) -> station series."""
    import scipy.io as sio
    from scipy.signal import fftconvolve
    from fillers import nearest_inside
    fn = f'{CACHE}/eden_station_{year}.csv'
    if os.path.exists(fn): print('cached', fn); return
    anyv, kk, reg5, den5, fill_c, fill_5, (PS, PN, PW, PE) = eden_setup()
    files = eden_files(year - 1, year)[3:]                        # Q4 of the previous year for the 30/90-day windows
    files = [f for f in files if os.path.exists(f)]
    cells, fives, dates = [], [], []
    for f in files:
        nc = sio.netcdf_file(f, 'r', mmap=False); s_ = nc.variables['stage'][:].astype('float32'); s_[(s_ < -900) | (s_ > 9000)] = np.nan
        tv = nc.variables['time']; units = tv.units.decode() if isinstance(tv.units, bytes) else tv.units
        t0 = pd.Timestamp(units.split('since')[1].strip()[:10]); dates += list(t0 + pd.to_timedelta(tv[:].astype(float), unit='D'))
        for d in range(s_.shape[0]): s_[d] = nearest_inside(s_[d], anyv)
        num = fftconvolve(np.nan_to_num(s_), kk[None], mode='same', axes=(1, 2))
        cells.append(s_); fives.append(np.where(reg5[None], num / np.maximum(den5, 1)[None], np.nan).astype('float32'))
    C = np.concatenate(cells); F5 = np.concatenate(fives); dates = pd.DatetimeIndex(dates).normalize()
    assert dates.is_monotonic_increasing and (np.diff(dates.values).astype('timedelta64[D]').astype(int) == 1).all(), 'EDEN dates not daily/continuous'
    def trailing(x, w, minp):
        v = np.isfinite(x); c = np.cumsum(np.where(v, x, 0), 0, dtype='float64'); n = np.cumsum(v, 0)
        cs = np.concatenate([np.zeros((1,) + x.shape[1:]), c]); ns = np.concatenate([np.zeros((1,) + x.shape[1:]), n])
        i = np.arange(1, len(x) + 1); j = np.maximum(i - w, 0); S = cs[i] - cs[j]; N = ns[i] - ns[j]
        return np.where(N >= minp, S / np.maximum(N, 1), np.nan).astype('float32')
    def blockwise(x, w, minp):                                   # row blocks of 45, as in eden_daily.py pass A (memory)
        o = np.empty(x.shape, 'float32')
        for r0 in range(0, x.shape[1], 45): o[:, r0:r0 + 45] = trailing(x[:, r0:r0 + 45], w, minp)
        return o
    FLD = {'cell30': (blockwise(C, 30, 10), None, anyv), '5km30': (blockwise(F5, 30, 10), None, reg5), '5km90': (blockwise(F5, 90, 30), None, reg5)}
    keep = dates.year == year; s = stations()
    y0p = 2790200 - PS * 400; x0p = 463400 - PW * 400
    ry = (s.Y.values - y0p) / 400; rx = (s.X.values - x0p) / 400
    iy = np.floor(ry).astype(int); ix = np.floor(rx).astype(int); fy = ry - iy; fx = rx - ix
    out = {'date': np.repeat(dates[keep], len(s)), 'station': np.tile(s.station.values, int(keep.sum()))}
    del C, F5, cells, fives
    # harmonic extrapolation evaluated only where needed: exact station weights from eden_weights.py
    import scipy.sparse as sp
    Wt = np.load(f'{HC}/eden_station_weights.npz', allow_pickle=True)
    assert (Wt['station'] == s.station.values).all()
    kp = np.flatnonzero(keep)
    for k, (A, fill, known) in FLD.items():
        nm = 'cell' if k == 'cell30' else '5km'
        ky, kx = Wt[f'{nm}_ky'], Wt[f'{nm}_kx']; R = sp.csr_matrix((Wt[f'{nm}_R_data'], Wt[f'{nm}_R_indices'], Wt[f'{nm}_R_indptr']), shape=tuple(Wt[f'{nm}_R_shape']))
        isk, ck, cu = Wt[f'{nm}_isk'], Wt[f'{nm}_ck'], Wt[f'{nm}_cu']
        Ak = A[kp]; dayok = np.isfinite(Ak).any(axis=(1, 2))
        kyn, kxn = ky - PS, kx - PW; inside = (kyn >= 0) & (kyn < Ak.shape[1]) & (kxn >= 0) & (kxn < Ak.shape[2])   # padded cells are NaN
        kv = np.full((len(kp), len(ky)), np.nan, 'float32'); kv[:, inside] = Ak[:, kyn[inside], kxn[inside]]; U = (R @ np.where(np.isfinite(kv), kv, 0).T).T          # (days, n_unknown_corners)
        U[~dayok] = np.nan
        vals = np.where(isk[None], kv[:, np.where(ck >= 0, ck, 0)], U[:, np.where(cu >= 0, cu, 0)])   # (days, stations, 4)
        w = np.stack([(1 - fy) * (1 - fx), (1 - fy) * fx, fy * (1 - fx), fy * fx], 1)[None]
        num = np.where(np.isfinite(vals), vals * w, 0).sum(2); den = np.where(np.isfinite(vals), w, 0).sum(2)
        out[{'cell30': 'stage_cell_cm30', '5km30': 'stage_5km_cm30', '5km90': 'stage_5km_cm90'}[k]] = np.where(den > 0, num / np.maximum(den, 1e-9), np.nan).ravel()
        del Ak, kv, U
    pd.DataFrame(out).to_csv(fn, index=False); print('eden', year, int(keep.sum()), 'days')


# ---------------------------------------------------------------- assemble
def step_assemble():
    s = stations(); ST = pd.read_csv(f'{CACHE}/static_at_stations.csv')
    days = pd.date_range('2000-01-01', '2025-12-31')
    # gridMET: rebuilt 1998-2013 for dates < 2013-01-01, original filled grids from 2013 on
    import zipfile
    GN = np.load(f'{WORK}/gridmet_filled_1998_2013.npz'); GO = np.load(f'{WORK}/gridmet_filled.npz')
    def gm(Z, k, la, lo, dsel):
        glat_s, glon = Z['glat_s'], Z['glon']; gd = pd.to_datetime(Z['dates'])
        iy = np.clip(np.searchsorted(glat_s, la) - 1, 0, len(glat_s) - 2); ix = np.clip(np.searchsorted(glon, lo) - 1, 0, len(glon) - 2)
        fy = np.clip((la - glat_s[iy]) / (glat_s[iy + 1] - glat_s[iy]), 0, 1); fx = np.clip((lo - glon[ix]) / (glon[ix + 1] - glon[ix]), 0, 1)
        ti = gd.get_indexer(dsel); A = Z[k][ti]
        v = [A[:, iy, ix], A[:, iy, ix + 1], A[:, iy + 1, ix], A[:, iy + 1, ix + 1]]; w = [(1 - fy) * (1 - fx), (1 - fy) * fx, fy * (1 - fx), fy * fx]
        num = sum(np.where(np.isfinite(a), a * b, 0) for a, b in zip(v, w)); den = sum(np.where(np.isfinite(a), b, 0) for a, b in zip(v, w))
        return np.where(den > 0, num / np.maximum(den, 1e-9), np.nan)
    d1 = days[days < '2013-01-01']; d2 = days[days >= '2013-01-01']
    T = pd.DataFrame({'date': np.repeat(days, len(s)), 'station': np.tile(s.station.values, len(days))})
    for k in ['P30', 'P90', 'P180', 'P365', 'R30', 'T90']:
        T[k] = np.concatenate([gm(GN, k, s.lat30.values, s.lon30.values, d1), gm(GO, k, s.lat30.values, s.lon30.values, d2)]).ravel()
    E = pd.concat([pd.read_csv(f, parse_dates=['date']) for f in sorted(glob.glob(f'{CACHE}/eden_station_*.csv'))])
    T = T.merge(E, on=['date', 'station'], how='left')
    # flows: rebuilt 1998-2013 + original 2013-2026 file
    fo = pd.read_csv(f'{SERVER}/dbhydro_flows/DBHydro_structure_flows_daily_2013_2026.csv', parse_dates=['date']).set_index('date')
    fn_ = pd.read_csv(f'{HC}/DBHydro_structure_flows_daily_1998_2013.csv', parse_dates=['date']).set_index('date')
    fl = pd.concat([fn_[fn_.index < '2013-01-01'], fo])[['SRS_inflow_cfs', 'TS_C111_detention_inflow_cfs']].reindex(pd.date_range('1998-01-01', '2025-12-31'))
    DRV = pd.DataFrame(index=fl.index)
    DRV['SRS90'] = fl.SRS_inflow_cfs.rolling(90, min_periods=72).mean(); DRV['TSC11190'] = fl.TS_C111_detention_inflow_cfs.rolling(90, min_periods=72).mean()
    # anchors: observed daily salinity (pre-2013 download + current training data, same processing)
    old = pd.read_csv(f'{HC}/salinity_daily_pre2013.csv', parse_dates=['date'])          # all pre-2013 observations (sonde + grab)
    new = pd.read_csv(f'{SERVER}/ENP_DBHydro_salinity_df.csv', dtype={'station': str}); new['date'] = pd.to_datetime(new.collectDate)
    new = new[(new.value >= 0) & (new.value <= 80) & (new.date >= '2013-04-01')].groupby(['station', 'date']).value.mean().rename('salinity').reset_index()
    sal = pd.concat([old.loc[old.grab == 0, ['station', 'date', 'salinity']], new]); piv = sal.pivot_table(index='date', columns='station', values='salinity').reindex(DRV.index)
    for a in ['ENPWP', 'ENPSR', 'ENPWW', 'ENPGB', 'TAYLORS3', 'ENPWB']:
        DRV[f'anc_{a}_30'] = piv[a].rolling(30, min_periods=10).mean(); DRV[f'anc_{a}_3'] = piv[a].rolling(3, min_periods=1).mean()
    DRV.index.name = 'date'; DRV.to_csv(f'{CACHE}/daily_drivers_1998_2025.csv')
    T = T.merge(DRV.reset_index(), on='date', how='left').merge(ST, on='station', how='left')
    doy = T.date.dt.dayofyear; T['doy_cos'] = np.cos(2 * np.pi * doy / 365.25); T['doy_sin'] = np.sin(2 * np.pi * doy / 365.25)
    T['src_flow90_w'] = (T.w_SRS + T.w_BigCypress) * T.SRS90 + (T.w_TaylorSlough + T.w_C111) * T.TSC11190
    obs = pd.concat([old.groupby(['station', 'date']).salinity.mean().reset_index().assign(src='pre2013'), new.assign(src='training')])
    T = T.merge(obs[['station', 'date', 'salinity', 'src']], on=['station', 'date'], how='left')
    T.to_parquet(f'{CACHE}/station_features_2000_2025.parquet', index=False) if _has_parquet() else T.to_csv(f'{CACHE}/station_features_2000_2025.csv.gz', index=False)
    print(T.shape, 'observed station-days:', int(T.salinity.notna().sum()))
    # check the rebuilt drivers against the training-period drivers file
    G = pd.read_csv(f'{HPC}/daily_global_drivers.csv', parse_dates=['date']).set_index('date')
    for c in ['SRS90', 'TSC11190', 'anc_ENPWP_30', 'anc_TAYLORS3_3']:
        j = pd.concat([DRV[c].rename('new'), G[c].rename('old')], axis=1).dropna(); print('check', c, len(j), 'max abs diff', round(float((j.new - j.old).abs().max()), 4))


def _has_parquet():
    try: import pyarrow; return True
    except ImportError: return False


if __name__ == '__main__':
    st = sys.argv[1]
    if st == 'static': step_static()
    elif st == 'eden': [step_eden(int(y)) for y in sys.argv[2:]]
    elif st == 'assemble': step_assemble()
