"""Download NOAA monthly mean sea level for the three tide gauges used as a predictor of the final model.
Key West 8724580, Vaca Key 8723970, Virginia Key 8723214 (CO-OPS API, product monthly_mean, station datum, metres).
The model uses the MSL column of each station relative to its own 2000-2025 mean (harness.py), so the datum does not matter.
Usage: python download_sea_level.py [out_dir]     (default: $ENP_DATA_ROOT/sea_level)
Writes noaa_monthly_mean_<station>.csv with the CO-OPS columns (Year, Month, Highest, MHHW, MHW, MSL, ...)."""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import os, sys, time, urllib.request
OUT = sys.argv[1] if len(sys.argv) > 1 else f'{ROOT}/sea_level'
STATIONS = {'8724580': 'Key West', '8723970': 'Vaca Key', '8723214': 'Virginia Key'}
API = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter'
os.makedirs(OUT, exist_ok=True)
for sid, name in STATIONS.items():
    url = (f'{API}?product=monthly_mean&station={sid}&begin_date=19980101&end_date=20251231&datum=STND&units=metric'
           f'&time_zone=GMT&format=csv&application=ENP_Landsat_Salinity_Model')
    for k in range(4):
        try:
            b = urllib.request.urlopen(url, timeout=120).read(); break
        except Exception as e:
            print('  retry', k + 1, e); time.sleep(10 * (k + 1))
    else:
        raise RuntimeError(url)
    if not b.startswith(b'Year'):
        raise RuntimeError(f'{sid}: unexpected response {b[:200]!r}')
    open(f'{OUT}/noaa_monthly_mean_{sid}.csv', 'wb').write(b)
    print(sid, name, b.count(b'\n') - 1, 'months')
