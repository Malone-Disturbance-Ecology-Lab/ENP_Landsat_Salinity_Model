"""Spatial attribution of day-to-day salinity variation to driver groups (final model), at the 250 m nodes.
Days: the 15th of every month, 2000-04..2025-12 (309 days). For each group G the group's predictors are replaced at each node by
that node's mean over the 309 days, and the drop in the node's temporal variance of predicted salinity is G's contribution."""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import os, sys, pickle, time, numpy as np, pandas as pd
os.environ['MODEL_PKL'] = os.path.abspath(f'{WORK}/final/final_model_onestep.pkl'); os.environ['MODEL_SO'] = os.path.abspath(f'{WORK}/final/m3_final.so')
OUTD = os.path.abspath(f'{WORK}/drivers'); os.makedirs(OUTD, exist_ok=True)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'final_model'))
import onestep_maps_final as om   # inputs in final_model/in/ (copy of hpc_inputs)
F = om.F3; C = {c: i for i, c in enumerate(F)}
G = {'water level': ['stage_cell_cm30', 'stage_5km_cm30', 'stage_5km_cm90'],
     'flows': ['SRS90', 'TSC11190', 'src_flow90_w'],
     'climate': ['P30', 'P90', 'P180', 'P365', 'R30', 'T90'],
     'anchor stations': [c for c in F if c.startswith('anc_')],
     'sea level': ['msl_1m', 'msl_3m', 'msl_12m'],
     'season': ['doy_cos', 'doy_sin']}
assert sum(len(v) for v in G.values()) == 29, sum(len(v) for v in G.values())
dom = np.zeros(om.H * om.W, bool); dom[om.IDX] = True; dom = dom.reshape(om.H, om.W)
ii, jj = np.nonzero(dom[4::8, 4::8][:om.H2, :om.W2]); keep = (ii % 2 == 0) & (jj % 2 == 0); ii, jj = ii[keep], jj[keep]; pos = np.searchsorted(om.IDX, (ii * 8 + 4) * om.W + (jj * 8 + 4))
st = np.ascontiguousarray(np.asarray(om.STA[:, pos]).T); wt = om.weights(om.LL[0][pos], om.LL[1][pos], om.UT[0][pos], om.UT[1][pos]); n = len(pos)
# 15th of each month, May 2000 - Dec 2025 (308 days; freq='MS' starts at 1 May, so April 2000 is not included).
# These are the days used for Fig. 10 / Sect. 3.7 -- keep as is to reproduce the published attribution.
days = pd.date_range('2000-04-15', '2025-12-15', freq='MS') + pd.Timedelta(days=14)
X = np.empty((len(days), n, len(F)), 'float32'); t0 = time.time()
for y in sorted(set(days.year)):
    E, ed = om.eden_year(y); ei = pd.Index(ed)
    for k in np.flatnonzero(days.year == y): X[k] = om.features(st, wt, days[k], E, ei.get_loc(days[k]), n)
    print(y, round(time.time() - t0), 's', flush=True)
np.save(f'{OUTD}/node_rc.npy', np.stack([ii, jj]))
m = pickle.load(open(os.environ['MODEL_PKL'], 'rb'))['m3']
def pred(A): return np.clip(m.predict(A.reshape(-1, len(F))), 0, None).reshape(len(days), n).astype('float32')
def predc(A): return np.clip(m.predict(A.reshape(-1, len(F))), 0, None).reshape(A.shape[0], A.shape[1]).astype('float32')
out = {'full': np.zeros(n, 'float32'), 'var_full': np.zeros(n, 'float32'), **{f'var_drop_{g}': np.zeros(n, 'float32') for g in G}}
for a in range(0, n, 4000):
    b = min(a + 4000, n); Xc = X[:, a:b]; P0 = predc(Xc); V0 = P0.var(0); out['full'][a:b] = P0.mean(0); out['var_full'][a:b] = V0
    for g, cols in G.items():
        A = Xc.copy(); idx = [C[c] for c in cols]; A[:, :, idx] = Xc[:, :, idx].mean(0, keepdims=True); out[f'var_drop_{g}'][a:b] = V0 - predc(A).var(0)
    print('chunk', a, round(time.time() - t0), 's', flush=True)
for g in G: print(g, 'median share', round(float(np.median(out[f'var_drop_{g}'] / np.maximum(out['var_full'], 1e-6))), 3))
np.savez_compressed(f'{OUTD}/node_attribution.npz', days=days.values.astype('datetime64[D]').astype(str), **out)
print('done', round(time.time() - t0), 's')
