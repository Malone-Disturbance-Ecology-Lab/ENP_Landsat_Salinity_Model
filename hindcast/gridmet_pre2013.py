"""gridMET climate windows for 1998-2013, built like gridmet_prep.py + make_hpc_inputs.py (same 22 x 32 crop, 7-row
southern pad, harmonic fill of ocean cells). Input: gm_{pr,srad,tmmn,tmmx}.nc from the gridMET THREDDS NCSS subset.
Output: gridmet_filled_1998_2013.npz (glat_s, glon, dates, P30, P90, P180, P365, R30, T90)."""
import numpy as np, pandas as pd, sys
from scipy.io import netcdf_file
from fillers import laplace_solver
def load(v):
    nc = netcdf_file(f'gm_{v}.nc', 'r', mmap=False); k = [k for k in nc.variables if len(nc.variables[k].shape) == 3][0]; a = nc.variables[k]
    raw = a.data; x = raw.astype('float32')
    fv = getattr(a, '_FillValue', None); mv = getattr(a, 'missing_value', None)
    bad = np.zeros(raw.shape, bool)
    for f in (fv, mv, 32767):
        if f is not None: bad |= raw == np.asarray(f).astype(raw.dtype)
    x[bad] = np.nan; x = x * getattr(a, 'scale_factor', 1.0) + getattr(a, 'add_offset', 0.0)
    lat = nc.variables['lat'][:].copy(); lon = nc.variables['lon'][:].copy(); t = pd.Timestamp('1900-01-01') + pd.to_timedelta(nc.variables['day'][:], unit='D')
    return x.astype('float32'), lat, lon, t
pr, lat, lon, t = load('pr'); srad, _, _, _ = load('srad'); tn, _, _, _ = load('tmmn'); tx, _, _, _ = load('tmmx'); tavg = (tn + tx) / 2
def roll(x, w, how):
    c = np.nancumsum(np.nan_to_num(x), 0); n = np.cumsum(np.isfinite(x), 0)
    s = c.copy(); s[w:] = c[w:] - c[:-w]; nn = n.copy(); nn[w:] = n[w:] - n[:-w]
    out = np.where(nn >= int(0.8 * w), (s if how == 'sum' else s / np.maximum(nn, 1)), np.nan); out[:w - 1] = np.nan; return out.astype('float32')
F = {'P30': roll(pr, 30, 'sum'), 'P90': roll(pr, 90, 'sum'), 'P180': roll(pr, 180, 'sum'), 'P365': roll(pr, 365, 'sum'), 'R30': roll(srad, 30, 'mean'), 'T90': roll(tavg, 90, 'mean')}
flip = lat[0] > lat[-1]; glat0 = lat[::-1] if flip else lat; step = np.diff(glat0).mean(); npad = 7; glat_s = np.r_[glat0[0] - step * np.arange(npad, 0, -1), glat0]
GA = {}
for k, A_ in F.items():
    A_ = A_[:, ::-1, :] if flip else A_; A_ = np.concatenate([np.full((A_.shape[0], npad, A_.shape[2]), np.nan, 'float32'), A_], 1)
    kn = np.isfinite(A_[-1]); GA[k] = laplace_solver(kn)(np.where(kn[None], A_, 0)).astype('float32'); GA[k][~np.isfinite(A_).any(axis=(1, 2))] = np.nan
np.savez_compressed('gridmet_filled_1998_2013.npz', glat_s=glat_s, glon=lon, dates=np.array([str(d.date()) for d in t]), **GA)
# overlap check against the training-period file (2013 dates)
Z = np.load('gridmet_filled.npz'); zd = pd.to_datetime(Z['dates']); td = pd.DatetimeIndex(t)
print('grid identical:', np.allclose(Z['glat_s'], glat_s), np.allclose(Z['glon'], lon))
for k in GA:
    i = td.get_indexer(zd[(zd >= '2013-01-01') & (zd <= '2013-12-31')]); j = zd.get_indexer(td[i])
    a = GA[k][i]; b = Z[k][j]; ok = np.isfinite(a) & np.isfinite(b)
    print(k, 'days', len(i), 'max abs diff', float(np.abs(a - b)[ok].max()), 'mean', float(np.abs(a - b)[ok].mean()))
