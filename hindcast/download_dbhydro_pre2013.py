"""Download DBHydro salinity before the current record (1980-01-01 to 2013-03-31) for the stations used by the model.
Continuous series: the same DBHydro Insights timeseries IDs as 00_download_dbhydro.R (one request per ID, zipped CSV).
Grab samples: water-quality test 98 (SALINITY) for the 43 grab stations.
Run from anywhere:  python3 download_dbhydro_pre2013.py <out_dir>"""
import sys, os, io, zipfile, time, json, urllib.request
OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
START, END = '19800101', '20130331'
API = 'https://insightsdata.api.sfwmd.gov/v1/insights-data'
IDS = ('89867,89865,89866,AL985,AL970,90111,90112,38108,63572,63578,63581,63584,63588,39533,63592,63596,63604,63608,63617,63621,'
       '63625,63630,63634,63638,63642,63646,63650,63654,63663,63667,63671,63675,63679,63683,63687,63691,AL633,63695,63699,63703,'
       '63708,39208,39209,39212,39215,AN690,38066,38067,66106,66107,66112,66113,66125,66126,90464,38010,38012,38011,66108,66109,'
       '66115,66117,66128,66129,38019,38021,38020,38023,38094,66110,66111,66123,66124,66130,66131').split(',')
GRAB = ['FLAB%02d' % k for k in (5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,27,29,30,31,33,34,35,36,37,38,39,40,41,43,44,47,48)] + \
       ['TTI51', 'TTI51B', 'TTI52', 'TTI55', 'TTI57', 'TTI59']

def get(url, data=None, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'} if data else {})
            with urllib.request.urlopen(req, timeout=300) as r: return r.read()
        except Exception as e:
            print('  retry', k + 1, e); time.sleep(10 * (k + 1))
    raise RuntimeError(url)

os.makedirs(f'{OUT}/raw_cont', exist_ok=True); os.makedirs(f'{OUT}/raw_grab', exist_ok=True)
for i in IDS:
    fn = f'{OUT}/raw_cont/{i}.csv'
    if os.path.exists(fn): continue
    b = get(f'{API}/cont/data?timeseriesIds={i}&reportType=timeseries&format=csv&startDate={START}&endDate={END}')
    if b[:2] == b'PK':
        z = zipfile.ZipFile(io.BytesIO(b)); b = z.read(z.namelist()[0]) if z.namelist() else b''   # empty zip = no data in period
    open(fn, 'wb').write(b); print(i, len(b)); time.sleep(2)
for s in GRAB:
    fn = f'{OUT}/raw_grab/{s}.csv'
    if os.path.exists(fn): continue
    q = json.dumps({'query': {'locations': [{'name': s, 'type': 'STATION'}], 'testNumbers': ['98']}}).encode()
    b = get(f'{API}/chem/report/data?reportType=timeseries&format=csv&startDate={START}&endDate={END}&excludeFlaggedResults=N', q)
    open(fn, 'wb').write(b); print(s, len(b)); time.sleep(2)
print('done')
