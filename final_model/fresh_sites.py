"""Freshwater training sites from DBHydro field specific conductance (test 9), 2000-2025.
1. Keep surface-water samples (matrix SW; sample types SAMP and RS) without a quality flag, with conductance >= 20 uS/cm
   (lower values are unit errors) and sampling depth <= 1 (as reported; blank depth kept).
2. Convert to practical salinity at 25 degC with TEOS-10 (gsw.SP_from_C, conductance in mS/cm) and average per station-day.
3. Station summary: number of sample days, mean and median salinity, first and last date, whether the station lies in the
   30 m map domain (px30_index.npy), and the distance to the nearest model station (great circle, km).
4. Training sites: inside the domain, at least 24 sample days, and more than 0.3 km from a model station (13 sites).
Usage: python fresh_sites.py [raw_dir] [out_dir]
  raw_dir default $ENP_DATA_ROOT/dbhydro_conductivity (cond_*.csv from download_conductivity.py)
  out_dir default work/exp; writes cond_daily.parquet, cond_station_summary.csv, fresh13_coords.csv.
Needs $ENP_DATA_ROOT/site_variables/station_site_variables.csv and $ENP_DATA_ROOT/best_two_step/salinity30m/hpc_inputs/px30_index.npy."""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import os, sys, glob, csv, numpy as np, pandas as pd, gsw
RAW = sys.argv[1] if len(sys.argv) > 1 else f'{ROOT}/dbhydro_conductivity'
OUT = sys.argv[2] if len(sys.argv) > 2 else f'{WORK}/exp'
IDXF = os.environ.get('PX30_INDEX', f'{ROOT}/best_two_step/salinity30m/hpc_inputs/px30_index.npy')
os.makedirs(OUT, exist_ok=True)

def read_insights(fn):
    """DBHydro Insights CSV: skip the '#' disclaimer block. A few rows have badly escaped quotes in the free-text comment
    columns at the end, which split them into extra fields; the measurement columns come first, so keep the first n fields."""
    with open(fn, encoding='utf-8', errors='replace', newline='') as f:
        lines = f.read().splitlines(True)
    h = next(i for i, l in enumerate(lines) if l.startswith('station,'))
    rows = list(csv.reader(lines[h:])); cols = rows[0]
    return pd.DataFrame([x[:len(cols)] for x in rows[1:] if len(x) >= len(cols)], columns=cols)

r = pd.concat([read_insights(f) for f in sorted(glob.glob(f'{RAW}/cond_*.csv'))], ignore_index=True).drop_duplicates()
v = pd.to_numeric(r.value, errors='coerce'); dep = pd.to_numeric(r.depth, errors='coerce')
keep = (r.matrix.str.strip() == 'SW') & r.sampleType.str.strip().isin(['SAMP', 'RS']) & (r.flag.str.strip() == '') & (v >= 20) \
       & (dep.isna() | (dep <= 1))
x = r[keep].copy(); x['sal'] = gsw.SP_from_C(v[keep].values / 1000.0, 25.0, 0.0)
x['date'] = pd.to_datetime(x.collectDate).dt.normalize()
x['lat'] = x.latitude.astype(float); x['lon'] = -x.longitude.astype(float).abs()      # DBHydro reports west longitude as positive
d = x.groupby(['station', 'date']).agg(sal=('sal', 'mean'), lat=('lat', 'first'), lon=('lon', 'first')).reset_index()
d = d[(d.date >= '2000-01-01') & (d.date <= '2025-12-31')]
d.to_parquet(f'{OUT}/cond_daily.parquet')

S = d.groupby('station').agg(n=('sal', 'size'), mean=('sal', 'mean'), med=('sal', 'median'), lat=('lat', 'first'), lon=('lon', 'first'),
                              first=('date', 'min'), last=('date', 'max'))
H, W = 3839, 4171; D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873
rr = np.floor((TOP - S.lat.values) / D30).astype(int); cc = np.floor((S.lon.values - LEFT) / D30).astype(int)
ok = (rr >= 0) & (rr < H) & (cc >= 0) & (cc < W)
S['inside'] = ok & np.isin(np.where(ok, rr * W + cc, -1), np.load(IDXF))
M = pd.read_csv(f'{ROOT}/site_variables/station_site_variables.csv')
la0, lo0 = np.radians(S.lat.values)[:, None], np.radians(S.lon.values)[:, None]
la1, lo1 = np.radians(M.latitude.values)[None], np.radians(M.longitude.values)[None]
S['km_to_model_station'] = (2 * 6371.0088 * np.arcsin(np.sqrt(np.sin((la1 - la0) / 2) ** 2 + np.cos(la0) * np.cos(la1) * np.sin((lo1 - lo0) / 2) ** 2))).min(1)
S.to_csv(f'{OUT}/cond_station_summary.csv')
new = S[S.inside & (S.n >= 24) & (S.km_to_model_station > 0.3)]
new[['lat', 'lon', 'n', 'mean']].to_csv(f'{OUT}/fresh13_coords.csv')
print(len(d), 'station-days at', len(S), 'stations;', len(new), 'training sites:', ', '.join(new.index))
