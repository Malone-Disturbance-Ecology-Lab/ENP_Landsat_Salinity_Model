"""Combine the pre-2013 DBHydro downloads into one daily table, processed like 01_harmonize_dbhydro.R + the model QC:
continuous series (all DBKEYs of a station) and grab samples (test 98) -> one value per station-day (mean of readings);
values < 0 or > 80 psu removed. Output: salinity_daily_pre2013.csv (station, date, salinity, grab, n_readings)."""
import pandas as pd, glob, io, os
rows = []
for f in sorted(glob.glob('raw_cont/*.csv')):
    s = open(f, errors='ignore').read()
    if not s.strip(): continue
    d = pd.read_csv(io.StringIO('\n'.join(l for l in s.splitlines() if not l.startswith('#'))), dtype={'STATION': str})
    rows.append(pd.DataFrame({'station': d.STATION.str.lstrip('0'), 'date': pd.to_datetime(d.TIMESTAMP).dt.normalize(),
                              'value': pd.to_numeric(d.VALUE, errors='coerce'), 'grab': 0}))
for f in sorted(glob.glob('raw_grab/*.csv')):
    s = open(f, errors='ignore').read()
    d = pd.read_csv(io.StringIO('\n'.join(l for l in s.splitlines() if not l.startswith('#'))), dtype=str)
    if d.empty: continue
    d = d[d.parameter.str.upper().eq('SALINITY') & d.sampleType.eq('SAMP')]
    rows.append(pd.DataFrame({'station': d.station, 'date': pd.to_datetime(d.collectDate).dt.normalize(),
                              'value': pd.to_numeric(d.value, errors='coerce'), 'grab': 1}))
A = pd.concat(rows).dropna(subset=['value'])
bad = (A.value < 0) | (A.value > 80); print('removed out-of-range values:', int(bad.sum()), sorted(A.value[bad].round(1).tolist())[:20])
A = A[~bad]
D = A.groupby(['station', 'date', 'grab']).value.agg(['mean', 'size']).reset_index().rename(columns={'mean': 'salinity', 'size': 'n_readings'})
D = D[D.date < '2013-04-01']
D.to_csv('salinity_daily_pre2013.csv', index=False)
S = D.groupby(['station', 'grab']).agg(start=('date', 'min'), end=('date', 'max'), days=('date', 'size')).reset_index()
S.to_csv('station_period_of_record_pre2013.csv', index=False)
print(len(D), 'station-days at', D.station.nunique(), 'stations'); print(S.groupby('grab').agg(stations=('station', 'size'), median_start=('start', 'median')))
print(D.assign(y=D.date.dt.year).groupby('y').station.nunique().to_string())
