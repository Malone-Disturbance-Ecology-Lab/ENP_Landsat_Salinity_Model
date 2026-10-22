"""Consistent 2000-2025 record from the FINAL one-step model (no vegetation spectra; sea level; 2026-10-13) with the surface-water mask, one tile per task:
 - monthly means over wet days -> seasonal Kendall / seasonal Theil-Sen trend of monthly anomalies, 2000-2025 (26 years; Jan-Mar 2000 missing)
 - high-salinity events: >= 7 consecutive wet days above the cell's own 95th percentile (2000-04..2025-12); a dry day ends an event
 - group accumulators (region / ecosystem) and daily series at salinity stations.
Days with mask code 0 (dry) or 9 (never wet) are missing; uncertain days (2) are kept.
Usage: python long_masked_hpc.py <tile> <root>  (tile 0..74 = 15 row blocks of 256 rows x 5 column blocks of 1024 cols)
Inputs: <root>/daily30m_onestep_final/YYYY/*.tif (2000-04..2025-12), <root>/wetdry30m/YYYY/*.tif, <root>/events_masked/in/{px30_index,eco30,region30}.npy,
        station_site_variables.csv.  Writes <root>/long_masked_final/parts/tile_XX.npz (+ stations_tile_XX.csv.gz)."""
import sys, os, time, numpy as np, pandas as pd, rasterio
from rasterio.windows import Window
from scipy.stats import norm, rankdata
tile = int(sys.argv[1]); R = sys.argv[2]; IN = f'{R}/events_masked/in'; OUT = f'{R}/long_masked_final/parts'; os.makedirs(OUT, exist_ok=True)
H, W = 3839, 4171; t0 = time.time()
r0 = (tile // 5) * 256; c0 = (tile % 5) * 1024; r1 = min(r0 + 256, H); c1 = min(c0 + 1024, W)
idx = np.load(f'{IN}/px30_index.npy'); eco_all = np.load(f'{IN}/eco30.npy'); reg_all = np.load(f'{IN}/region30.npy')
a, b = np.searchsorted(idx, [r0 * W, r1 * W]); sel = np.arange(a, b); p = idx[sel]; rr, cc = p // W, p % W
k = (cc >= c0) & (cc < c1); sel, rr, cc = sel[k], rr[k], cc[k]; n = len(sel)
if n == 0: np.savez(f'{OUT}/tile_{tile:02d}.npz', n=0); sys.exit(0)
loc = (rr - r0) * (c1 - c0) + (cc - c0)
ALLM = pd.date_range('2000-04-01', '2025-12-01', freq='MS')
MS = int(os.environ.get('M_START', '0')); MC = int(os.environ.get('M_COUNT', str(len(ALLM))))
MLIST = ALLM[MS:MS + MC]; NM = len(MLIST)
edates = pd.date_range(MLIST[0], MLIST[-1] + pd.offsets.MonthEnd(0)); T = len(edates)
Y0 = 2000; NYR = 26; yrs = edates.year.values - Y0; mon = edates.month.values
X = np.empty((T, n), 'float32'); MM = np.full((NM, n), np.nan, 'float32'); t = 0
win = Window(c0, r0, c1 - c0, r1 - r0)
for m, M0 in enumerate(MLIST):
    with rasterio.open(f'{R}/daily30m_onestep_final/{M0.year}/ENP_salinity30m_onestep_daily_{M0.year}_{M0.month:02d}.tif') as s: S = s.read(window=win)
    with rasterio.open(f'{R}/wetdry30m/{M0.year}/ENP_wetdry30m_daily_{M0.year}_{M0.month:02d}.tif') as q: C = q.read(window=win)
    nd = S.shape[0]; S = S.reshape(nd, -1)[:, loc]; C = C.reshape(nd, -1)[:, loc]
    x = S.astype('float32') / 10; x[(S < 0) | (C == 0) | (C == 9)] = np.nan; X[t:t + nd] = x; t += nd
    with np.errstate(all='ignore'): MM[m] = np.nanmean(x, 0)
assert t == T
print('read', round(time.time() - t0), 's', flush=True)
# ---------- seasonal Kendall on monthly wet-day means, 26 calendar years (missing months NaN)
NY = NYR; iu, ju = np.triu_indices(NY, 1); DT = (ju - iu).astype('float32')
SETS = {'monthly': list(range(12)), 'drymonths': [11, 0, 1, 2, 3, 4], 'wetmonths': [5, 6, 7, 8, 9, 10]}
def seasonal_kendall(Xs):
    n_ = Xs.shape[2]; d = Xs[ju] - Xs[iu]; sg = np.sign(d); sg[np.isnan(sg)] = 0; sg = sg.astype(np.int8); Sg = sg.sum(0, dtype=np.int32)
    ng = np.isfinite(Xs).sum(0); Xf = np.where(np.isfinite(Xs), Xs, np.inf); srt = np.sort(Xf, 0); eq = (srt[1:] == srt[:-1]) & np.isfinite(srt[1:])
    run = np.zeros(Xs.shape[1:], 'float32'); tie = np.zeros(Xs.shape[1:], 'float32')
    for kk in range(eq.shape[0]):
        e = eq[kk]; tt = run + 1; tie += np.where((~e) & (run > 0), tt * (tt - 1) * (2 * tt + 5), 0); run = np.where(e, run + 1, 0)
    tt = run + 1; tie += np.where(run > 0, tt * (tt - 1) * (2 * tt + 5), 0); varg = (ng * (ng - 1) * (2 * ng + 5) - tie) / 18.0
    r = rankdata(Xf, axis=0).astype('float32'); Rk = np.where(np.isfinite(Xs), r, ((ng + 1) / 2.0)[None])
    K = np.empty((12, 12, n_), 'float32')
    for g in range(12): K[g] = (sg[:, g][:, None, :] * sg).sum(0, dtype=np.int32)
    RR = np.einsum('ign,ihn->ghn', Rk, Rk); ng1 = (ng + 1).astype('float32'); cov = (K + 4 * RR - NY * ng1[:, None, :] * ng1[None, :, :]) / 3.0
    for g in range(12): cov[g, g] = varg[g]
    slopes = d / DT[:, None, None]; out = {}
    for kname, ms in SETS.items():
        Sv = Sg[ms].sum(0).astype('float32'); V = np.maximum(cov[np.ix_(ms, ms)].sum((0, 1)), 1e-6)
        z = np.where(Sv > 0, (Sv - 1) / np.sqrt(V), np.where(Sv < 0, (Sv + 1) / np.sqrt(V), 0))
        with np.errstate(all='ignore'): med = np.nanmedian(slopes[:, ms].reshape(-1, n_), 0)
        out[kname] = dict(slope=med.astype('float32'), z=z.astype('float32'), p=(2 * norm.sf(np.abs(z))).astype('float32'))
    return out
full = np.full((NY * 12, n), np.nan, 'float32'); off = (MLIST[0].year - Y0) * 12 + MLIST[0].month - 1; full[off:off + NM] = MM
TR = {f'{q}_{kname}': np.full(n, np.nan, 'float32') for kname in SETS for q in ('slope', 'z', 'p')}
for s0 in range(0, n, 8000):
    s1 = min(s0 + 8000, n); o = seasonal_kendall(full[:, s0:s1].reshape(NY, 12, -1))
    for kname, v in o.items():
        for q in ('slope', 'z', 'p'): TR[f'{q}_{kname}'][s0:s1] = v[q]
print('trend', round(time.time() - t0), 's', flush=True)
# ---------- events
st = pd.read_csv(f'{IN}/station_site_variables.csv'); D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873
srow = np.floor((TOP - st.latitude.values) / D30).astype(int); scol = np.floor((st.longitude.values - LEFT) / D30).astype(int)
gpos = rr * W + cc; spos = dict(zip(srow * W + scol, st.station)); SSER = {spos[g]: X[:, q].copy() for q, g in enumerate(gpos) if g in spos}
with np.errstate(all='ignore'):
    q95 = np.nanpercentile(X, 95, axis=0).astype('float32'); q50 = np.nanmedian(X, axis=0).astype('float32')
nwet = np.isfinite(X).sum(0); hyp = np.stack([(X[yrs == y] > 40).sum(0) for y in range(NYR)])
ex = (X > q95[None]).astype(np.int8)
d = np.diff(np.vstack([np.zeros((1, n), np.int8), ex, np.zeros((1, n), np.int8)]), axis=0)
sp, stt = np.nonzero(d.T == 1); _, ent = np.nonzero(d.T == -1); dur = ent - stt; kk = dur >= 7; sp, stt, ent, dur = sp[kk], stt[kk], ent[kk], dur[kk]
eco = eco_all[sel]; reg = reg_all[sel]; nE, nR = 7, 6
acc = {}
for g, cat, nn in (('reg', reg, nR), ('eco', eco, nE)):
    for kname in ('ev_year', 'evd_year', 'px_event_year', 'hyp_year'): acc[f'{g}_{kname}'] = np.zeros((nn, NYR))
    acc[f'{g}_npx'] = np.bincount(cat[cat >= 0], minlength=nn).astype(float)
    acc[f'{g}_msum'] = np.zeros((nn, NM)); acc[f'{g}_mcnt'] = np.zeros((nn, NM)); acc[f'{g}_wetcells'] = np.zeros((nn, NM))
    for c_ in range(nn):
        mm = cat == c_
        if mm.any():
            v = MM[:, mm]; ok = np.isfinite(v); acc[f'{g}_msum'][c_] = np.where(ok, v, 0).sum(1); acc[f'{g}_mcnt'][c_] = ok.sum(1)
            acc[f'{g}_hyp_year'][c_] += hyp[:, mm].sum(1)
EVY = np.zeros(n, 'float32')
if len(sp):
    y0 = yrs[stt]; EVY = (np.bincount(sp, minlength=n) / (T / 365.25)).astype('float32')
    for g, cat, nn in (('reg', reg, nR), ('eco', eco, nE)):
        cg = cat[sp]; okg = cg >= 0
        np.add.at(acc[f'{g}_ev_year'], (cg[okg], y0[okg]), 1); np.add.at(acc[f'{g}_evd_year'], (cg[okg], y0[okg]), dur[okg])
        pe = np.zeros((n, NYR), bool); pe[sp, y0] = True
        for c_ in range(nn):
            mm = cat == c_
            if mm.any(): acc[f'{g}_px_event_year'][c_] += pe[mm].sum(0)
np.savez_compressed(f'{OUT}/tile_{tile:02d}.npz', n=n, sel=sel, p95=q95, p50=q50, events_per_yr=EVY, wet_days=nwet.astype('int32'),
                    hyper_days_per_yr=(hyp.sum(0) / (T / 365.25)).astype('float32'), months=np.array([str(m.date()) for m in MLIST]), **TR, **acc)
if SSER: pd.DataFrame(SSER, index=edates).astype('float32').to_csv(f'{OUT}/stations_tile_{tile:02d}.csv.gz')
print('done tile', tile, 'cells', n, 'events', len(sp), round(time.time() - t0), 's', flush=True)
