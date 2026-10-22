"""Download DBHydro field specific conductance (water-quality test 9, all projects and stations) for 2000-2025.
The 13 freshwater training sites are selected from these data by fresh_sites.py.
DBHydro Insights chemistry report (same API as hindcast/download_dbhydro_pre2013.py), in four periods to keep requests small.
Usage: python download_conductivity.py [out_dir]     (default: $ENP_DATA_ROOT/dbhydro_conductivity)"""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import os, sys, time, json, urllib.request
OUT = sys.argv[1] if len(sys.argv) > 1 else f'{ROOT}/dbhydro_conductivity'
API = 'https://insightsdata.api.sfwmd.gov/v1/insights-data/chem/report/data'
PERIODS = [('20000101', '20061231'), ('20070101', '20131231'), ('20140101', '20191231'), ('20200101', '20251231')]
os.makedirs(OUT, exist_ok=True)
q = json.dumps({'query': {'testNumbers': ['9']}}).encode()
for a, b in PERIODS:
    fn = f'{OUT}/cond_test9_{a[:4]}_{b[:4]}.csv'
    if os.path.exists(fn): print('exists', fn); continue
    url = f'{API}?reportType=timeseries&format=csv&startDate={a}&endDate={b}&excludeFlaggedResults=N'
    for k in range(4):
        try:
            req = urllib.request.Request(url, data=q, headers={'Content-Type': 'application/json'})
            body = urllib.request.urlopen(req, timeout=900).read(); break
        except Exception as e:
            print('  retry', k + 1, e); time.sleep(30 * (k + 1))
    else:
        raise RuntimeError(url)
    open(fn, 'wb').write(body); print(fn, round(len(body) / 1e6, 1), 'MB'); time.sleep(5)
