"""Cross-calibrated coverage of the 90 % prediction intervals (Sect. 4.2).

The released half-widths (final_model/q_levels.json) are the 90th percentile of the
absolute block-CV residuals in four predicted-salinity bins, pooled over all ten
spatial blocks. Using them to score the same residuals is optimistic, so here each
block is scored with half-widths recomputed from the residuals of the other nine
blocks only (stations + 13 fresh sites in the calibration pool).

Writes interval_coverage_crosscal.csv (overall + by bin, model stations) and
interval_coverage_crosscal_detail.csv (by era, fresh sites, per-station).
"""
import os, json
import numpy as np, pandas as pd

WORK = os.environ.get('ENP_WORK', 'work')
C = pd.read_parquet(f'{WORK}/final/cv_predictions_final.parquet')
C = C[C.salinity.notna()].copy()
Q = json.load(open(f'{WORK}/final/q_levels.json'))
bins = np.array(Q['bins'])
C['bin'] = np.clip(np.digitize(C.pred_cv, bins[1:-1]), 0, 3)
C['ae'] = (C.salinity - C.pred_cv).abs()

# leave-one-block-out half-widths
C['cov'] = False
for k in sorted(C.block.unique()):
    pool = C[C.block != k]
    q = np.array([np.quantile(pool.ae[pool.bin == b], 0.9) for b in range(4)])
    m = C.block == k
    C.loc[m, 'cov'] = C.ae[m].values <= q[C.bin[m].values]

M = C[~C.fresh_site]
out = pd.Series({'overall': M['cov'].mean(),
                 **{f'bin{b}': M['cov'][M.bin == b].mean() for b in range(4)}})
out.to_csv('interval_coverage_crosscal.csv')

era = np.where(M.date.dt.year <= 2012, '2000-2012', '2013-2025')
st = M.groupby('station')['cov'].agg(['mean', 'size'])
st = st[st['size'] >= 30]
detail = {**{f'era_{e}': v for e, v in M.groupby(era)['cov'].mean().items()},
          'fresh_sites': C['cov'][C.fresh_site].mean(),
          'stations_lt70pct': int((st['mean'] < 0.7).sum()),
          'stations_scored': len(st)}
pd.Series(detail).to_csv('interval_coverage_crosscal_detail.csv')
print(out.round(4).to_string()); print(pd.Series(detail).to_string())
