"""One-step (no-spectra) salinity maps, 2000-04 .. 2025-12, made the same way as the published 30 m product:
  1. daily 250 m maps (119,019 nodes = 30 m cells 8i+4, 8j+4) from the one-step model m3;
  2. the full one-step model at every 30 m cell on the 15th of each month -> detail C = full30 - bilinear(250 m daily);
  3. monthly 30 m = bilinear(250 m monthly mean) + monthly mean of C interpolated in time (as export30.py).
Inputs (in/): static30.npy, px30_index/lonlat30/utm30.npy, gridmet_filled*.npz, daily_drivers_1998_2025.csv,
eden_anyv.npy, q_levels.json, final_models.pkl; EDEN v3 quarterly NetCDFs (eden/); compiled model m3.so.
Run:  python onestep_maps.py run      (resumable by year)   then   python onestep_maps.py export"""
import sys, os, glob, json, time, pickle, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scipy.io as sio, rasterio
from scipy import ndimage as ndi
from scipy.signal import fftconvolve
from rasterio.transform import Affine
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fillers import laplace_solver, nearest_inside
B = os.path.dirname(os.path.abspath(__file__)); IN = f'{B}/in'; WK = f'{B}/work'; OUT = os.environ.get('MAP_OUT', f'{B}/out'); os.makedirs(WK, exist_ok=True)
H, W = 3839, 4171; D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873; H2, W2 = 480, 521
DAYS = pd.date_range('2000-04-01', '2025-12-31'); MONTHS = pd.date_range('2000-04-01', '2025-12-01', freq='MS'); MID = MONTHS + pd.Timedelta(days=14)
NM, ND = len(MONTHS), len(DAYS)
from model_config import F0; F3 = list(F0)   # inputs only (used by fresh_features.py); predictions use onestep_maps_final.py
def pred(X): raise RuntimeError('onestep_maps.py builds model inputs only; use onestep_maps_final.py for predictions')
STK = json.load(open(f'{IN}/static30_keys.json')); STA = np.load(f'{IN}/static30.npy', mmap_mode='r')
IDX = np.load(f'{IN}/px30_index.npy'); LL = np.load(f'{IN}/lonlat30.npy', mmap_mode='r'); UT = np.load(f'{IN}/utm30.npy', mmap_mode='r'); NP = len(IDX)
DRV = pd.read_csv(f'{IN}/daily_drivers_1998_2025.csv', parse_dates=['date']).set_index('date')
DCOLS = ['SRS90', 'TSC11190'] + [c for c in DRV.columns if c.startswith('anc_')]
GN = np.load(f'{IN}/gridmet_filled_1998_2013.npz'); GO = np.load(f'{IN}/gridmet_filled.npz')
glat_s, glon = GN['glat_s'], GN['glon']; assert np.allclose(glat_s, GO['glat_s']) and np.allclose(glon, GO['glon'])
GD = {'N': pd.DatetimeIndex(pd.to_datetime(GN['dates'])), 'O': pd.DatetimeIndex(pd.to_datetime(GO['dates']))}
CLIM = ['P30', 'P90', 'P180', 'P365', 'R30', 'T90']
GA = {'N': {k: GN[k] for k in CLIM}, 'O': {k: GO[k] for k in CLIM}}          # load once (npz members decompress on every access)
def gfield(k, d):
    s = 'N' if d < pd.Timestamp('2013-01-01') else 'O'; return GA[s][k][GD[s].get_loc(d)]
# ---- EDEN geometry (as eden_daily.py)
anyv = np.load(f'{IN}/eden_anyv.npy')
r_ = 12; yy, xx = np.mgrid[-r_:r_ + 1, -r_:r_ + 1]; KK = ((yy ** 2 + xx ** 2) <= 12.5 ** 2).astype('float32')
reg5 = ndi.convolve(anyv.astype('float32'), KK, mode='constant') > 0; den5 = ndi.convolve(anyv.astype('float32'), KK, mode='constant')
PS, PW, PE, PN = 140, 70, 10, 10; y0p = 2790200 - PS * 400; x0p = 463400 - PW * 400
X30, Y30 = UT[0], UT[1]
er0 = max(int(np.floor((Y30.min() - y0p) / 400)) - 2, 0); er1 = int(np.ceil((Y30.max() - y0p) / 400)) + 3
ec0 = max(int(np.floor((X30.min() - x0p) / 400)) - 2, 0); ec1 = int(np.ceil((X30.max() - x0p) / 400)) + 3
FILL = {'c': laplace_solver(np.pad(anyv, ((PS, PN), (PW, PE)))), '5': laplace_solver(np.pad(reg5, ((PS, PN), (PW, PE))))}
KNOWN = {'c': np.pad(anyv, ((PS, PN), (PW, PE))), '5': np.pad(reg5, ((PS, PN), (PW, PE)))}
EK = {'stage_cell_cm30': 'c', 'stage_5km_cm30': '5', 'stage_5km_cm90': '5'}

def bil(A2, fy, fx, iy, ix):
    v = [A2[iy, ix], A2[iy, ix + 1], A2[iy + 1, ix], A2[iy + 1, ix + 1]]; w = [(1 - fy) * (1 - fx), (1 - fy) * fx, fy * (1 - fx), fy * fx]
    num = sum(np.where(np.isfinite(a), a * b, 0) for a, b in zip(v, w)); den = sum(np.where(np.isfinite(a), b, 0) for a, b in zip(v, w))
    return np.where(den > 0, num / np.maximum(den, 1e-9), np.nan).astype('float32')

import scipy.sparse as sp
def _interp_matrix(fy, fx, iy, ix, ncol_grid, ncells):
    n = len(fy); rows = np.repeat(np.arange(n), 4)
    cols = np.stack([iy * ncol_grid + ix, iy * ncol_grid + ix + 1, (iy + 1) * ncol_grid + ix, (iy + 1) * ncol_grid + ix + 1], 1).ravel()
    vals = np.stack([(1 - fy) * (1 - fx), (1 - fy) * fx, fy * (1 - fx), fy * fx], 1).ravel()
    return sp.csr_matrix((vals.astype('float32'), (rows, cols)), shape=(n, ncells))
def weights(lon, lat, X, Y):
    """Sparse bilinear interpolation matrices (NaN-aware use: (M @ v0) / (M @ ok)) for gridMET and the cropped EDEN grid."""
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
def features(st, wt, d, E, di, n):
    """st: (n, 27) static block in STK order."""
    X = np.empty((n, len(F3)), 'float32'); X[:, SI] = st
    for k in CLIM: X[:, COLI[k]] = sbil(wt['g'], gfield(k, d))
    for k in EK: X[:, COLI[k]] = sbil(wt['e'], E[k][di])
    r = DRV.loc[d]
    for c in DCOLS: X[:, COLI[c]] = r[c]
    doy = d.dayofyear; X[:, COLI['doy_cos']] = np.cos(2 * np.pi * doy / 365.25); X[:, COLI['doy_sin']] = np.sin(2 * np.pi * doy / 365.25)
    w = st[:, STK.index('w_SRS')] + st[:, STK.index('w_BigCypress')]; v = st[:, STK.index('w_TaylorSlough')] + st[:, STK.index('w_C111')]
    X[:, COLI['src_flow90_w']] = w * X[:, COLI['SRS90']] + v * X[:, COLI['TSC11190']]
    return X

def eden_year(y):
    """Filled, cropped daily EDEN window fields for every day of year y: dict k -> (ndays, ch, cw)."""
    files = [f'{B}/eden/eden_{yy}_q{q}.nc' for yy, q in [(y - 1, 4)] + [(y, q) for q in range(1, 5)]]
    files = [f for f in files if os.path.exists(f)]
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
    del C, F5
    keep = np.flatnonzero(dates.year == y); E = {}
    for k, A in NAT.items():
        out = np.full((len(keep), er1 - er0, ec1 - ec0), np.nan, 'float32')
        for c0 in range(0, len(keep), 30):
            g = np.pad(A[keep[c0:c0 + 30]], ((0, 0), (PS, PN), (PW, PE)), constant_values=np.nan); ok = np.isfinite(g).any(axis=(1, 2))
            res = FILL[EK[k]](np.where(np.isfinite(g), g, 0)).astype('float32'); res[~ok] = np.nan
            res[~np.isfinite(g) & np.isfinite(res) & KNOWN[EK[k]][None]] = np.nan
            out[c0:c0 + 30] = res[:, er0:er1, ec0:ec1]
        E[k] = out
    return E, dates[keep]

def run():
    t0 = time.time()
    # 250 m nodes
    dom = np.zeros(H * W, bool); dom[IDX] = True; dom = dom.reshape(H, W)
    ii, jj = np.nonzero(dom[4::8, 4::8][:H2, :W2]); flat = (ii * 8 + 4) * W + (jj * 8 + 4); pos = np.searchsorted(IDX, flat)
    stN = np.ascontiguousarray(np.asarray(STA[:, pos]).T); wtN = weights(LL[0][pos], LL[1][pos], UT[0][pos], UT[1][pos]); nN = len(pos)
    print('nodes', nN, flush=True)
    D250 = np.lib.format.open_memmap(f'{WK}/daily250_int16.npy', 'r+' if os.path.exists(f'{WK}/daily250_int16.npy') else 'w+', 'int16', (ND, H2, W2))
    C30 = np.lib.format.open_memmap(f'{WK}/C30_int16.npy', 'r+' if os.path.exists(f'{WK}/C30_int16.npy') else 'w+', 'int16', (NM, NP))
    global TILE
    TILE = {a: (np.ascontiguousarray(np.asarray(STA[:, a:min(a + 1_000_000, NP)]).T),
                weights(np.asarray(LL[0, a:min(a + 1_000_000, NP)]), np.asarray(LL[1, a:min(a + 1_000_000, NP)]), np.asarray(UT[0, a:min(a + 1_000_000, NP)]), np.asarray(UT[1, a:min(a + 1_000_000, NP)])))
            for a in range(0, NP, 1_000_000)}
    print('30 m tiles cached', round(time.time() - t0), 's', flush=True)
    donef = f'{WK}/years_done.txt'; done = set(open(donef).read().split()) if os.path.exists(donef) else set()
    # 30 m bilinear weights from the 250 m grid (as run30.py)
    for y in sorted(set(DAYS.year)):
        if str(y) in done: continue
        ty = time.time(); E, ed = eden_year(y); eidx = pd.Index(ed); print(y, 'eden', round(time.time() - ty), 's', flush=True)
        dsel = np.flatnonzero(DAYS.year == y)
        for b0 in range(0, len(dsel), 8):                                        # 8 days per prediction batch
            ds = dsel[b0:b0 + 8]; X = np.concatenate([features(stN, wtN, DAYS[d], E, eidx.get_loc(DAYS[d]), nN) for d in ds])
            p = np.clip(pred(X), 0, None).reshape(len(ds), nN)
            for q, d in enumerate(ds):
                g = np.full((H2, W2), -1, 'int16'); g[ii, jj] = np.round(p[q] * 10).astype('int16'); D250[d] = g
        print(y, '250 m daily done', round(time.time() - ty), 's', flush=True)
        for m in np.flatnonzero(MONTHS.year == y):
            tm = time.time(); d = MID[m]; di = eidx.get_loc(d); g = np.asarray(D250[DAYS.get_loc(d)], 'float32') / 10; g[g < 0] = np.nan
            for a in range(0, NP, 1_000_000):
                b = min(a + 1_000_000, NP); n = b - a; pix = IDX[a:b]; rr, cc = pix // W, pix % W
                st, wt = TILE[a]
                full = np.clip(pred(features(st, wt, d, E, di, n)), 0, None)
                fi = (rr - 4) / 8.0; miy = np.clip(np.floor(fi).astype(int), 0, H2 - 2); mfy = np.clip(fi - miy, 0, 1)
                fj = (cc - 4) / 8.0; mix = np.clip(np.floor(fj).astype(int), 0, W2 - 2); mfx = np.clip(fj - mix, 0, 1)
                b250 = bil(g, mfy, mfx, miy, mix)
                C30[m, a:b] = np.where(np.isfinite(b250), np.clip(np.round((full - b250) * 10), -30000, 30000), -32768).astype('int16')
            print(f'  {d.date()} 30 m {round(time.time() - tm)} s', flush=True)
        D250.flush(); C30.flush(); open(donef, 'a').write(f'{y}\n'); print(y, 'done', round(time.time() - ty), 's; total', round(time.time() - t0), 's', flush=True)

def export():
    """GeoTIFFs: daily 250 m (one file per month, one band per day), monthly 250 m mean / lower90 / upper90, monthly 30 m mean."""
    D250 = np.load(f'{WK}/daily250_int16.npy', mmap_mode='r'); C30 = np.load(f'{WK}/C30_int16.npy', mmap_mode='r')
    Q = json.load(open(f'{IN}/q_levels.json')); bins = np.array(Q['bins']); q = np.array(Q['q'])
    T250 = Affine(D30 * 8, 0, LEFT, 0, -D30 * 8, TOP); T30 = Affine(D30, 0, LEFT, 0, -D30, TOP)
    base = dict(driver='GTiff', crs='EPSG:4326', compress='deflate', predictor=2, tiled=True, blockxsize=256, blockysize=256)
    model = 'one-step model (v3; no vegetation spectra), trained 2013-2025; hindcast before 2013-04 (no-spectra hindcast, 2026-10-09)'
    inside = np.zeros(H * W, bool); inside[IDX] = True; inside = inside.reshape(H, W)
    midn = MID.values.astype('datetime64[D]').astype(float); dn = DAYS.values.astype('datetime64[D]').astype(float)
    J = np.clip(np.searchsorted(midn, dn) - 1, 0, NM - 2); Fw = np.clip((dn - midn[J]) / (midn[J + 1] - midn[J]), 0, 1)
    fi = (np.arange(H) - 4) / 8.0; iy = np.clip(np.floor(fi).astype(int), 0, H2 - 2); fy = np.clip(fi - iy, 0, 1)[:, None]
    fj = (np.arange(W) - 4) / 8.0; ix = np.clip(np.floor(fj).astype(int), 0, W2 - 2); fx = np.clip(fj - ix, 0, 1)[None, :]
    def up(g):
        v = [g[iy][:, ix], g[iy][:, ix + 1], g[iy + 1][:, ix], g[iy + 1][:, ix + 1]]; w = [(1 - fy) * (1 - fx), (1 - fy) * fx, fy * (1 - fx), fy * fx]
        num = sum(np.where(np.isfinite(a), a * b, 0) for a, b in zip(v, w)); den = sum(np.where(np.isfinite(a), b, 0) for a, b in zip(v, w))
        return np.where(den > 0, num / np.maximum(den, 1e-9), np.nan).astype('float32')
    for d_ in ['daily250', 'monthly250', 'monthly30']: os.makedirs(f'{OUT}/{d_}', exist_ok=True)
    MM = {k: [] for k in ('mean', 'lower90', 'upper90')}
    UNTIL = pd.Timestamp(os.environ.get('EXPORT_UNTIL', '2100-01-01'))
    for m, M0 in enumerate(MONTHS):
        if M0 > UNTIL: break
        days = np.flatnonzero((DAYS >= M0) & (DAYS < M0 + pd.offsets.MonthBegin(1)))
        raw = np.asarray(D250[days[0]:days[-1] + 1]); os.makedirs(f'{OUT}/daily250/{M0.year}', exist_ok=True)
        with rasterio.open(f'{OUT}/daily250/{M0.year}/ENP_salinity_onestep_daily_{M0.year}_{M0.month:02d}.tif', 'w', width=W2, height=H2, count=len(days), dtype='int16', nodata=-1, transform=T250, **base) as o:
            o.write(raw); [o.set_band_description(i + 1, str(DAYS[d].date())) for i, d in enumerate(days)]
            o.update_tags(units='salinity psu x 10 (divide by 10); -1 = no data', model=model + '; 250 m')
        p = raw.astype('float32') / 10; p[raw < 0] = np.nan; hw = q[np.clip(np.digitize(p, bins[1:-1]), 0, 3)]
        with np.errstate(all='ignore'):
            MM['mean'].append(np.nanmean(p, 0)); MM['lower90'].append(np.nanmean(np.clip(p - hw, 0, None), 0)); MM['upper90'].append(np.nanmean(p + hw, 0))
        # monthly 30 m (export30.py): bilinear(250 m monthly mean) + time-weighted mean of the mid-month detail layers
        wts = {}
        for d in days: wts[J[d]] = wts.get(J[d], 0) + (1 - Fw[d]); wts[J[d] + 1] = wts.get(J[d] + 1, 0) + Fw[d]
        cm = np.zeros(NP, 'float32'); ws = np.zeros(NP, 'float32')
        for k, wk in wts.items():
            if wk <= 0: continue
            c = np.asarray(C30[k], 'float32'); ok = c > -32768; cm += np.where(ok, c / 10 * wk, 0); ws += np.where(ok, wk, 0)
        cmg = np.zeros(H * W, 'float32'); cmg[IDX] = np.where(ws > 0, cm / np.maximum(ws, 1e-9), 0); cmg = cmg.reshape(H, W)
        s = np.clip(up(MM['mean'][-1]) + cmg, 0, None); o30 = np.where(inside & np.isfinite(s), np.round(s * 10), -1).astype('int16')
        with rasterio.open(f'{OUT}/monthly30/ENP_salinity30m_onestep_monthly_{M0.year}_{M0.month:02d}.tif', 'w', width=W, height=H, count=1, dtype='int16', nodata=-1, transform=T30, **base) as o:
            o.write(o30, 1); o.set_band_description(1, f'{M0.year}-{M0.month:02d}')
            o.update_tags(units='monthly mean salinity, psu x 10 (divide by 10); -1 = no data', grid='HLS/Landsat 30 m (EPSG:4326, 0.000270923 deg)', model=model + '; 30 m')
        if m % 12 == 0: print('export', M0.date(), flush=True)
    for k, L in MM.items():
        A = np.stack(L).astype('float32')
        for y in sorted(set(MONTHS.year)):
            sel = np.flatnonzero(MONTHS.year == y); sel = sel[sel < len(A)]
            if not len(sel): continue
            with rasterio.open(f'{OUT}/monthly250/ENP_salinity_onestep_monthly_{k}_{y}.tif', 'w', width=W2, height=H2, count=len(sel), dtype='float32', nodata=np.nan, transform=T250, **{**base, 'predictor': 3}) as o:
                o.write(A[sel]); [o.set_band_description(i + 1, f'{MONTHS[s].year}-{MONTHS[s].month:02d}') for i, s in enumerate(sel)]
                o.update_tags(units='psu', model=model + '; 250 m', interval='90% split-conformal half-widths by predicted level ' + str(np.round(q, 2).tolist()) + ' from ' + Q['source'])
    print('export done')

if __name__ == '__main__':
    if 'test' in sys.argv:                                    # quick timing test: one month of 250 m days + one 30 m tile
        E, ed = eden_year(2005); eidx = pd.Index(ed)
        dom = np.zeros(H * W, bool); dom[IDX] = True; dom = dom.reshape(H, W)
        ii, jj = np.nonzero(dom[4::8, 4::8][:H2, :W2]); pos = np.searchsorted(IDX, (ii * 8 + 4) * W + (jj * 8 + 4))
        st = np.ascontiguousarray(np.asarray(STA[:, pos]).T); wt = weights(LL[0][pos], LL[1][pos], UT[0][pos], UT[1][pos])
        t = time.time(); X = features(st, wt, pd.Timestamp('2005-05-15'), E, eidx.get_loc(pd.Timestamp('2005-05-15')), len(pos)); t1 = time.time() - t
        t = time.time(); p = pred(X); t2 = time.time() - t
        a=2_000_000; b=a+1_000_000; st30=np.ascontiguousarray(np.asarray(STA[:, a:b]).T); wt30=weights(np.asarray(LL[0,a:b]),np.asarray(LL[1,a:b]),np.asarray(UT[0,a:b]),np.asarray(UT[1,a:b]))
        t=time.time(); X30=features(st30, wt30, pd.Timestamp('2005-05-15'), E, eidx.get_loc(pd.Timestamp('2005-05-15')), 1_000_000); t3=time.time()-t; t=time.time(); p30=pred(X30); t4=time.time()-t
        print('30 m tile 1e6: features', round(t3,2), 's predict', round(t4,2), 's; mean', round(float(np.clip(p30,0,None).mean()),2))
        print('nodes', len(pos), 'features', round(t1, 2), 's predict', round(t2, 2), 's', 'NaN share by feature', {c: round(float(np.isnan(X[:, i]).mean()), 3) for i, c in enumerate(F3) if np.isnan(X[:, i]).mean() > 0})
        print('2005-05-15 node salinity: mean', round(float(p.mean()), 2), 'p5/p95', np.round(np.percentile(p, [5, 95]), 1))
    elif 'run' in sys.argv: run()
    elif 'export' in sys.argv: export()
