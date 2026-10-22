"""Daily structure-flow aggregates 1998-2013, built exactly like dbhydro_flows/build_flows_daily.py."""
import pandas as pd, glob
d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob('flows_raw/dbhydro_*_Q_1998_2013.csv')) if open(f).read().count('\n') > 1])
d['date'] = pd.to_datetime(d.date)
d['col'] = d.site + '_' + d.parameter.str.extract(r'\((\w+)\)')[0].str[0]
w = d.pivot_table(index='date', columns='col', values='value', aggfunc='mean').sort_index()
w = w.rename(columns=lambda c: c.split('_')[0] if c != 'S178_W' else 'S178_weir')
w['S178'] = w['S178'].fillna(0) + w.pop('S178_weir').fillna(0)
w = w.mask(w.abs() > 5000)
for c in ['S12B','S12C','S12D','S333','S333N','S334','S332B','S332C','S332D','S199','S200','S18C','S197']:
    if c not in w: w[c] = float('nan')
z = w.fillna(0)
w['SRS_inflow_cfs'] = z[['S12B','S12C','S12D','S333','S333N']].sum(1) - z['S334']
w['TS_C111_detention_inflow_cfs'] = z[['S332B','S332C','S332D','S199','S200']].sum(1)
w['C111_S18C_cfs'] = z['S18C']; w['C111_S197_cfs'] = z['S197']
w.index.name = 'date'; w.round(3).to_csv('DBHydro_structure_flows_daily_1998_2013.csv')
o = pd.read_csv('../dbhydro_flows/DBHydro_structure_flows_daily_2013_2026.csv', parse_dates=['date']).set_index('date')
j = pd.concat([w.SRS_inflow_cfs.rename('new'), o.SRS_inflow_cfs.rename('old')], axis=1).dropna(); print('SRS overlap days', len(j), 'max abs diff', round((j.new - j.old).abs().max(), 3))
j = pd.concat([w.TS_C111_detention_inflow_cfs.rename('new'), o.TS_C111_detention_inflow_cfs.rename('old')], axis=1).dropna(); print('TSC111 overlap days', len(j), 'max abs diff', round((j.new - j.old).abs().max(), 3))
y = w[['SRS_inflow_cfs','TS_C111_detention_inflow_cfs','S12A']].groupby(w.index.year).mean().round(0); print(y.to_string())
print('training 2013-2025 mean SRS', round(o.SRS_inflow_cfs.mean()), 'TSC111', round(o.TS_C111_detention_inflow_cfs.mean()))
