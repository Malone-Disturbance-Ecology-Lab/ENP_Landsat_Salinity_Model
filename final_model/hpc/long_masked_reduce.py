"""Combine long_masked_final tiles (final one-step model, 2026-10-13): pixel results (with BH-FDR q), trend tables by region/ecosystem, event and series tables, station series."""
import sys, glob, os, numpy as np, pandas as pd
R = sys.argv[1]; P = f'{R}/long_masked_final/parts'; D = f'{R}/long_masked_final/results'; os.makedirs(D, exist_ok=True)
idx = np.load(f'{R}/events_masked/in/px30_index.npy'); NP = len(idx); eco = np.load(f'{R}/events_masked/in/eco30.npy'); reg = np.load(f'{R}/events_masked/in/region30.npy')
SETS = ['monthly', 'drymonths', 'wetmonths']; PK = ['p95', 'p50', 'events_per_yr', 'wet_days', 'hyper_days_per_yr'] + [f'{q}_{s}' for s in SETS for q in ('slope', 'z', 'p')]
PIX = {k: np.full(NP, np.nan, 'float32') for k in PK}; ACC = {}; months = None
files = sorted(glob.glob(f'{P}/tile_*.npz')); assert len(files) == 75, len(files)
for f in files:
    z = np.load(f)
    if int(z['n']) == 0: continue
    s = z['sel']; months = z['months']
    for k in PK: PIX[k][s] = z[k]
    for k in z.files:
        if k.startswith(('reg_', 'eco_')): ACC[k] = ACC.get(k, 0) + z[k]
def bh(p):
    ok = np.isfinite(p); q = np.full(p.shape, np.nan, 'float32'); ps = p[ok]; o = np.argsort(ps); n = len(ps)
    qq = ps[o] * n / np.arange(1, n + 1); qq = np.minimum.accumulate(qq[::-1])[::-1]; t = np.empty(n, 'float32'); t[o] = np.minimum(qq, 1); q[ok] = t; return q
for s in SETS: PIX[f'q_{s}'] = bh(PIX[f'p_{s}'])
np.savez_compressed(f'{D}/pixel_long_masked_2000_2025.npz', **PIX)
ECO = ['Freshwater marsh & prairie', 'Swamp forest & scrub', 'Salt marsh', 'Mangrove', 'Inland & estuarine water', 'Open bay', 'Other']
REG = ['Big Cypress', 'Shark River Slough', 'Taylor Slough', 'C-111', 'Gulf coastal waters', 'Florida Bay']; A = 0.000818
rows = []
for gname, cat, labels in (('ecosystem', eco, ECO[:6]), ('region', reg, REG)):
    for i, l in enumerate(labels):
        m = cat == i; d = dict(group=gname, name=l, pixels=int(m.sum()))
        for s in SETS:
            sl = PIX[f'slope_{s}'][m]; q = PIX[f'q_{s}'][m]; pp = PIX[f'p_{s}'][m]
            d[f'{s}_median_slope'] = float(np.nanmedian(sl)); d[f'{s}_p10_slope'] = float(np.nanpercentile(sl, 10)); d[f'{s}_p90_slope'] = float(np.nanpercentile(sl, 90))
            d[f'{s}_pct_sig_up'] = 100 * float(((q < 0.05) & (sl > 0)).mean()); d[f'{s}_pct_sig_down'] = 100 * float(((q < 0.05) & (sl < 0)).mean())
            d[f'{s}_km2_sig_up'] = round(float(((q < 0.05) & (sl > 0)).sum() * A), 1); d[f'{s}_km2_sig_down'] = round(float(((q < 0.05) & (sl < 0)).sum() * A), 1)
            d[f'{s}_pct_p05'] = 100 * float((pp < 0.05).mean())
        rows.append(d)
pd.DataFrame(rows).round(4).to_csv(f'{D}/monthly_anomaly_trends_by_group_2000_2025.csv', index=False)
Y, Ms = [], []
for g, labels in (('reg', REG), ('eco', ECO)):
    npx = ACC[f'{g}_npx']
    for i, l in enumerate(labels):
        for y in range(26):
            Y.append(dict(group='region' if g == 'reg' else 'ecosystem', name=l, year=2000 + y, events_per_pixel=ACC[f'{g}_ev_year'][i, y] / npx[i],
                          event_days_per_pixel=ACC[f'{g}_evd_year'][i, y] / npx[i], pct_area_with_event=100 * ACC[f'{g}_px_event_year'][i, y] / npx[i],
                          hyper_days_mean=ACC[f'{g}_hyp_year'][i, y] / npx[i]))
        for j, mo in enumerate(months):
            c = ACC[f'{g}_mcnt'][i, j]; Ms.append(dict(group='region' if g == 'reg' else 'ecosystem', name=l, month=mo, mean_salinity_wet=ACC[f'{g}_msum'][i, j] / c if c else np.nan,
                                                       pct_cells_with_surface_water=100 * c / npx[i]))
pd.DataFrame(Y).round(3).to_csv(f'{D}/high_salinity_events_by_group_year_2000_2025.csv', index=False)
pd.DataFrame(Ms).round(3).to_csv(f'{D}/monthly_series_by_group_2000_2025.csv', index=False)
st = [pd.read_csv(f, index_col=0, parse_dates=True) for f in sorted(glob.glob(f'{P}/stations_tile_*.csv.gz'))]
if st:
    S = pd.concat(st, axis=1); S.index.name = 'date'; S.stack().rename('model_salinity_30m_onestep_masked').reset_index().rename(columns={'level_1': 'station'}).to_csv(f'{D}/station_daily_onestep_masked_2000_2025.csv.gz', index=False)
T = pd.read_csv(f'{D}/monthly_anomaly_trends_by_group_2000_2025.csv'); print(T[['group', 'name', 'monthly_median_slope', 'monthly_pct_sig_up', 'monthly_pct_sig_down', 'monthly_pct_p05']].to_string())
