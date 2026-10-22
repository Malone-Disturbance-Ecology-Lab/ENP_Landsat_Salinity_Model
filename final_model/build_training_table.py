"""Training table of the final model: station table 2000-2025 + 13 freshwater sites, with sea level and spatial blocks.
- T.parquet: station-day inputs and salinity at the non-anchor gauges, 2000-04 .. 2025-12 (harness.load(): sea level, blocks).
- TA_index.parquet: station, date, salinity, block and is_new for T plus the freshwater sites. Blocks are the 10 k-means
  blocks of the gauges; each freshwater site goes to the nearest block centroid.
Inputs: $ENP_DATA_ROOT/dbhydro_hindcast/features/station_features_2000_2025_v2.parquet, $ENP_DATA_ROOT/sea_level/,
$ENP_DATA_ROOT/site_variables/station_site_variables.csv, and $ENP_WORK/exp/{fresh_station_features.parquet,
cond_station_summary.csv} from fresh_sites.py and fresh_features.py. Writes to $ENP_WORK/exp/."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from sklearn.cluster import KMeans
from harness import load, U
WORK = os.environ.get('ENP_WORK', 'work'); E = f'{WORK}/exp'; os.makedirs(E, exist_ok=True)
T = load(); T.to_parquet(f'{E}/T.parquet'); print('T', T.shape, int(T.salinity.notna().sum()), 'observed station-days', flush=True)
F = pd.read_parquet(f'{E}/fresh_station_features.parquet')
S = pd.read_csv(f'{U}/site_variables/station_site_variables.csv').set_index('station'); st = sorted(set(T.station))
km = KMeans(10, n_init=10, random_state=0).fit(S.loc[st, ['longitude', 'latitude']].values)
assert (km.predict(S.loc[st, ['longitude', 'latitude']].values) == T.drop_duplicates('station').set_index('station').loc[st, 'block'].values).all()
NS = pd.read_csv(f'{E}/cond_station_summary.csv', index_col=0); F['block'] = km.predict(NS.loc[F.station, ['lon', 'lat']].values)
TA = pd.concat([T, F[T.columns.intersection(F.columns)]], ignore_index=True); TA['is_new'] = TA.station.isin(F.station.unique())
TA[['station', 'date', 'salinity', 'block', 'is_new']].to_parquet(f'{E}/TA_index.parquet')
print('TA', TA.shape, 'fresh-site blocks:', F.groupby('station').block.first().to_dict())
