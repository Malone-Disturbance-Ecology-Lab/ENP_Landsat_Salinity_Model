"""Grouped permutation importance of the final model, held out in space (10-block CV): for each block, fit on the other nine,
then shuffle one predictor group jointly across the held-out rows (3 repeats) and record the rise in RMSE and fall in R2,
overall and by station salinity zone (mean < 5, 5-25, > 25 psu)."""
import sys, numpy as np, pandas as pd, pickle, time; sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', 'final_model'))
from harness import HGB, PRM, F0
E = 'work/exp'; MS = ['msl_1m', 'msl_3m', 'msl_12m']; FE = F0 + MS
TA = pd.read_parquet(f'{E}/TA_index.parquet'); T = pd.read_parquet(f'{E}/T.parquet'); Fr = pd.read_parquet(f'{E}/fresh_station_features.parquet')
ref = T.drop_duplicates('date').set_index('date')
for c in MS: Fr[c] = Fr.date.map(ref[c])
Fr['block'] = TA[TA.is_new].block.values
A = pd.concat([T, Fr[T.columns.intersection(Fr.columns)]], ignore_index=True); A = A[A.salinity.notna()].reset_index(drop=True)
G = {'water level': ['stage_cell_cm30', 'stage_5km_cm30', 'stage_5km_cm90'], 'flows': ['SRS90', 'TSC11190', 'src_flow90_w'],
     'climate': ['P30', 'P90', 'P180', 'P365', 'R30', 'T90'], 'anchor stations': [c for c in FE if c.startswith('anc_')],
     'sea level': MS, 'season': ['doy_cos', 'doy_sin']}
G['site variables'] = [c for c in FE if not any(c in v for v in G.values())]
print({k: len(v) for k, v in G.items()}, sum(len(v) for v in G.values()), len(FE))
zone = A.groupby('station').salinity.transform('mean'); A['zone'] = pd.cut(zone, [-1, 5, 25, 99], labels=['< 5 psu', '5–25 psu', '> 25 psu'])
rng = np.random.default_rng(0); X = A[FE].astype('float32').values; y = A.salinity.values; col = {c: i for i, c in enumerate(FE)}
base = np.full(len(A), np.nan); perm = {g: np.full((3, len(A)), np.nan) for g in G}; t0 = time.time()
for b in range(10):
    tr = (A.block != b).values; te = ~tr; m = HGB(**PRM).fit(X[tr], y[tr]); Xt = X[te]; base[te] = np.clip(m.predict(Xt), 0, None)
    for g, cols in G.items():
        idx = [col[c] for c in cols]
        for r in range(3):
            Xp = Xt.copy(); Xp[:, idx] = Xt[rng.permutation(len(Xt))][:, idx]; perm[g][r, te] = np.clip(m.predict(Xp), 0, None)
    print('block', b, round(time.time() - t0), 's', flush=True)
def sc(p, s):
    e = p[s] - y[s]; return np.sqrt(np.mean(e ** 2)), 1 - np.sum(e ** 2) / np.sum((y[s] - y[s].mean()) ** 2)
rows = []
for zn, s in [('all', np.ones(len(A), bool))] + [(z, (A.zone == z).values) for z in A.zone.cat.categories]:
    r0, q0 = sc(base, s)
    for g in G:
        d = [sc(perm[g][r], s) for r in range(3)]
        rows.append(dict(zone=zn, group=g, n_pred=len(G[g]), rmse_base=r0, r2_base=q0, rmse_rise=np.mean([x[0] for x in d]) - r0, r2_drop=q0 - np.mean([x[1] for x in d])))
R = pd.DataFrame(rows); R.round(4).to_csv('permutation_importance_groups.csv', index=False)
print(R.pivot(index='group', columns='zone', values='rmse_rise').round(2).to_string())
np.savez_compressed('perm_preds.npz', base=base)
