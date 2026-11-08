# Build a daily structure-flow table from DBHydro Insights downloads (raw_downloads/*.csv).
# Downloaded 2026-09-25 from https://insights.api.sfwmd.gov/v1/insights/sites/<SITE>/data?startDate=2013-01-01&endDate=2026-02-28&hydroParams=Q&wqParams=0
import pandas as pd, glob
fs=sorted(glob.glob('raw_downloads/dbhydro_*_Q_2013_2026.csv'))
d=pd.concat([pd.read_csv(f) for f in fs])
d['date']=pd.to_datetime(d.date)
d['col']=d.site+'_'+d.parameter.str.extract(r'\((\w+)\)')[0].str[0]
w=d.pivot_table(index='date',columns='col',values='value',aggfunc='mean').sort_index()
w=w.rename(columns=lambda c: c.split('_')[0] if c!='S178_W' else 'S178_weir')
w['S178']=w['S178'].add(w.pop('S178_weir'),fill_value=0)   # NaN only when both are missing
# QC: physically impossible spikes (S12D 53,658 cfs on 2021-01-15 and 40,896 on 2022-08-01)
w=w.mask(w.abs()>5000)
# Aggregates. A structure with no value that day (gap, or masked spike) contributes 0 to the sum, as before,
# but a day on which *every* component is missing is NaN rather than 0, so the 90-day rolling means
# downstream (min_periods=72) can see the gap. n_missing_* counts the absent components per day.
SRS=['S12B','S12C','S12D','S333','S333N']; TSC=['S332B','S332C','S332D','S199','S200']
w['SRS_inflow_cfs']=w[SRS].sum(axis=1,min_count=1)-w['S334'].fillna(0)   # Shark River Slough: Tamiami Trail structures minus S334 (east to L-31N); S12A not available
w['TS_C111_detention_inflow_cfs']=w[TSC].sum(axis=1,min_count=1)          # pumps into ENP east-boundary detention areas (Taylor Slough / C-111 spreader)
w['C111_S18C_cfs']=w['S18C']      # lower C-111 toward Florida Bay / Manatee Bay (NaN kept)
w['C111_S197_cfs']=w['S197']      # C-111 outlet to Manatee Bay / Barnes Sound (NaN kept)
w['n_missing_SRS']=w[SRS+['S334']].isna().sum(axis=1); w['n_missing_TSC111']=w[TSC].isna().sum(axis=1)
w.index.name='date'
w.round(3).to_csv('DBHydro_structure_flows_daily_2013_2026.csv')
print(w.shape, w.SRS_inflow_cfs.max())
