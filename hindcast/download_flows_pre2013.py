"""Daily structure flows 1998-01-01..2013-12-31 from DBHydro Insights (same endpoint as dbhydro_flows/raw_downloads).
Writes flows_raw/dbhydro_<SITE>_Q_1998_2013.csv with columns site,date,parameter,value (like the 2013-2026 files).
Also pulls S332, S174, S175 and S12A (not in the model's sums) to document how deliveries changed."""
import json, urllib.request, time, os, sys, pandas as pd
OUT = sys.argv[1] if len(sys.argv) > 1 else 'flows_raw'
SITES = ['S12B','S12C','S12D','S176','S177','S178','S18C','S197','S199','S200','S331','S332B','S332C','S332D','S333N','S333','S334',
         'S332','S174','S175','S12A']
for s in SITES:
    fn = f'{OUT}/dbhydro_{s}_Q_1998_2013.csv'
    if os.path.exists(fn): continue
    parts = []
    for a, b in [('1998-01-01', '2005-12-31'), ('2006-01-01', '2013-12-31')]:
        for k in range(4):
            try:
                j = json.load(urllib.request.urlopen(f'https://insights.api.sfwmd.gov/v1/insights/sites/{s}/data?startDate={a}&endDate={b}&hydroParams=Q&wqParams=0', timeout=300)); break
            except Exception as e: print(' retry', s, e); time.sleep(10)
        t = pd.DataFrame(j.get('tabledata', []))
        if len(t): parts.append(pd.DataFrame({'site': s, 'date': t.timestamp.str[:10], 'parameter': t.parameter, 'timeseriesId': t.timeseriesId, 'value': t.value, 'qc': t.qc}))
        time.sleep(1)
    d = pd.concat(parts) if parts else pd.DataFrame(columns=['site','date','parameter','timeseriesId','value','qc'])
    d.to_csv(fn, index=False); print(s, len(d), d.date.min() if len(d) else '', d.parameter.unique()[:4] if len(d) else '')
print('done')
