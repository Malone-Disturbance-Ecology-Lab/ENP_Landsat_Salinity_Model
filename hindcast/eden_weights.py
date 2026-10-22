"""Station weights for the EDEN harmonic (Laplace) extrapolation, so station values equal the full-grid fill of eden_daily.py.
fill(g)[unknown] = A^-1 Bk g[known]  ->  value at an unknown cell u = r_u . g[known], r_u = Bk^T A^-1 e_u (A symmetric).
Only known cells next to the unknown region get weight. Output: eden_station_weights.npz"""
import numpy as np, pandas as pd, scipy.sparse as sp, time
from scipy.sparse.linalg import splu
from scipy import ndimage as ndi
from pyproj import Transformer
anyv = np.load('eden_anyv.npy')
r_ = 12; yy, xx = np.mgrid[-r_:r_ + 1, -r_:r_ + 1]; kk = ((yy ** 2 + xx ** 2) <= 12.5 ** 2).astype('float32')
reg5 = ndi.convolve(anyv.astype('float32'), kk, mode='constant') > 0
PS, PW, PE, PN = 140, 70, 10, 10; y0p = 2790200 - PS * 400; x0p = 463400 - PW * 400
H30, W30 = 3839, 4171; D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873
s = pd.read_csv('station_site_variables.csv')[['station', 'longitude', 'latitude']]
row = np.floor((TOP - s.latitude) / D30).astype(int); col = np.floor((s.longitude - LEFT) / D30).astype(int)
lon30 = LEFT + (col + 0.5) * D30; lat30 = TOP - (row + 0.5) * D30
X, Y = Transformer.from_crs(4326, 26917, always_xy=True).transform(lon30.values, lat30.values)
ry = (Y - y0p) / 400; rx = (X - x0p) / 400; iy = np.floor(ry).astype(int); ix = np.floor(rx).astype(int); fy = ry - iy; fx = rx - ix
out = dict(station=s.station.values, iy=iy, ix=ix, fy=fy, fx=fx)
for name, mask in (('cell', anyv), ('5km', reg5)):
    known = np.pad(mask, ((PS, PN), (PW, PE))); H, W = known.shape; unk = ~known
    idx = -np.ones((H, W), int); uy, ux = np.nonzero(unk); n = len(uy); idx[uy, ux] = np.arange(n)
    ky, kx = np.nonzero(known); kidx = -np.ones((H, W), int); kidx[ky, kx] = np.arange(len(ky))
    rows = []; cols = []; vals = []; krow = []; kcol = []; diag = np.zeros(n)
    for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):            # identical to fillers.laplace_solver
        ny, nx = uy + dy, ux + dx; ok = (ny >= 0) & (ny < H) & (nx >= 0) & (nx < W); diag[ok] += 1
        i = np.nonzero(ok)[0]; nyy, nxx = ny[i], nx[i]; u = unk[nyy, nxx]
        rows += list(i[u]); cols += list(idx[nyy[u], nxx[u]]); vals += [-1.0] * int(u.sum())
        krow += list(i[~u]); kcol += list(kidx[nyy[~u], nxx[~u]])
    A = sp.csc_matrix((np.r_[vals, diag], (np.r_[rows, np.arange(n)], np.r_[cols, np.arange(n)])), shape=(n, n))
    Bk = sp.csr_matrix((np.ones(len(krow)), (krow, kcol)), shape=(n, len(ky)))
    t0 = time.time(); lu = splu(A); print(name, 'unknown', n, 'known', len(ky), 'splu', round(time.time() - t0), 's', flush=True)
    corners = [(iy, ix), (iy, ix + 1), (iy + 1, ix), (iy + 1, ix + 1)]
    cy = np.concatenate([c[0] for c in corners]); cx = np.concatenate([c[1] for c in corners])
    uc = np.unique(idx[cy, cx][unk[cy, cx]])
    E_ = np.zeros((n, len(uc))); E_[uc, np.arange(len(uc))] = 1
    Z = lu.solve(E_, trans='T')                                    # A^-T e_u
    R = (Bk.T @ Z).T                                               # (n_uc, n_known)
    R[np.abs(R) < 1e-9] = 0; Rs = sp.csr_matrix(R)
    print(name, 'unknown corner cells', len(uc), 'nonzero weights', Rs.nnz, 'row sums', np.round(np.asarray(Rs.sum(1)).ravel()[:5], 6), flush=True)
    # per station corner: known -> ('k', known index); unknown -> ('u', row of R)
    ck = np.stack([kidx[c[0], c[1]] for c in corners], 1); cu = np.stack([np.searchsorted(uc, idx[c[0], c[1]]) if True else 0 for c in corners], 1)
    isk = np.stack([known[c[0], c[1]] for c in corners], 1)
    out[f'{name}_ky'] = ky; out[f'{name}_kx'] = kx; out[f'{name}_isk'] = isk; out[f'{name}_ck'] = ck; out[f'{name}_cu'] = np.where(isk, -1, cu)
    out[f'{name}_R_data'] = Rs.data.astype('float64'); out[f'{name}_R_indices'] = Rs.indices; out[f'{name}_R_indptr'] = Rs.indptr; out[f'{name}_R_shape'] = np.array(Rs.shape)
    # self-check: random field -> full fill vs weights at the unknown corners
    import fillers
    g = np.random.default_rng(0).normal(size=(H, W)); full = fillers.laplace_solver(known)(np.where(known, g, 0))
    approx = Rs @ g[ky, kx]; exact = full[uy[uc], ux[uc]]; print(name, 'self-check max abs diff', float(np.abs(approx - exact).max()), flush=True)
np.savez_compressed('eden_station_weights.npz', **out); print('saved')
