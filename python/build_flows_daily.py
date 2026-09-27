# Build a daily structure-flow table from DBHydro Insights downloads (raw_downloads/*.csv).
# Downloaded 2026-09-25 from https://insights.api.sfwmd.gov/v1/insights/sites/<SITE>/data?startDate=2013-01-01&endDate=2026-02-28&hydroParams=Q&wqParams=0
import pandas as pd, glob
fs=sorted(glob.glob('raw_downloads/dbhydro_*_Q_2013_2026.csv'))
d=pd.concat([pd.read_csv(f) for f in fs])
d['date']=pd.to_datetime(d.date)
d['col']=d.site+'_'+d.parameter.str.extract(r'\((\w+)\)')[0].str[0]
w=d.pivot_table(index='date',columns='col',values='value',aggfunc='mean').sort_index()
w=w.rename(columns=lambda c: c.split('_')[0] if c!='S178_W' else 'S178_weir')
w['S178']=w['S178'].fillna(0)+w.pop('S178_weir').fillna(0)
# QC: physically impossible spikes (S12D 53,658 cfs on 2021-01-15 and 40,896 on 2022-08-01)
w=w.mask(w.abs()>5000)
z=w.fillna(0)
# Aggregates (missing daily values treated as 0 in sums)
w['SRS_inflow_cfs']=z[['S12B','S12C','S12D','S333','S333N']].sum(1)-z['S334']   # Shark River Slough: Tamiami Trail structures minus S334 (east to L-31N); S12A not available
w['TS_C111_detention_inflow_cfs']=z[['S332B','S332C','S332D','S199','S200']].sum(1)  # pumps into ENP east-boundary detention areas (Taylor Slough / C-111 spreader)
w['C111_S18C_cfs']=z['S18C']      # lower C-111 toward Florida Bay / Manatee Bay
w['C111_S197_cfs']=z['S197']      # C-111 outlet to Manatee Bay / Barnes Sound
w.index.name='date'
w.round(3).to_csv('DBHydro_structure_flows_daily_2013_2026.csv')
print(w.shape, w.SRS_inflow_cfs.max())
