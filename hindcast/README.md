# dbhydro_hindcast — pre-2013 gauge data and the no-spectra hindcast (2026-10-09)

> Scripts in this folder build the pre-2013 gauge, flow and climate inputs and the 2000–2025 station feature table used by the final model. The data files listed below are on the lab server in `dbhydro_hindcast/`; the hindcast evaluation scripts of the earlier model are in `deprecated/hindcast_released_model/`.

| path | content |
|---|---|
| `raw_cont/<DBKEY>.csv` | DBHydro continuous salinity, 1980 → 2013-03-31 (same timeseries IDs as `00_download_dbhydro.R`; empty = no data) |
| `raw_grab/<STATION>.csv` | DBHydro grab samples (test 98, SALINITY) for the 43 grab stations |
| `salinity_daily_pre2013.csv` | daily station salinity (mean of readings; QC 0–80 psu), the same processing as the training data |
| `station_period_of_record_pre2013.csv` | first/last date and number of days per station |
| `flows_raw/`, `DBHydro_structure_flows_daily_1998_2013.csv` | structure flows 1998–2013 and the model aggregates (as `dbhydro_flows/build_flows_daily.py`) |
| `features/gridmet_raw/` | gridMET THREDDS subsets (pr, srad, tmmn, tmmx; 1998–2013; same 22 × 32 crop) + `gridmet_pre2013.py` |
| `features/gridmet_filled_1998_2013.npz` | filled gridMET climate windows 1998–2013 |
| `features/daily_drivers_1998_2025.csv` | SRS90, TSC11190 and anchor-station 3/30-day means, 1998–2025 |
| `features/static_at_stations.csv` | the 27 static layers at each station's 30 m cell |
| `features/station_features_2000_2025_v2.parquet` | one-step model inputs + observed salinity at all 90 stations, 2000–2025 (use this one) |
| `features/station_features_2000_2025_v1_superseded.parquet` | first version: pre-2013 grab samples were missing (bug); superseded by v2 |
| `features/nodes250_static.npz` | static layers at the 119,019 250 m map nodes (not used by the final maps) |
| `eden_anyv.npy`, `eden_station_weights.npz` | EDEN domain mask and exact harmonic-fill weights at the station cells (`scripts/eden_weights.py`) |
| `results/` | station hindcast: daily predictions 2000–2025 (released model, block CV, block CV without flows), skill by era, figure |
| `scripts/` | download, harmonise, feature, hindcast and figure scripts; `m3.so` is the one-step model compiled with treelite/tl2cgen |
| `transfer/` | gzipped EDEN copies used to move data to the cloud workspace; safe to delete (the files are in `data/EDEN_GRIDS_Q`) |

EDEN v3 quarterly NetCDFs for 2000–2012 were added to `data/EDEN_GRIDS_Q/` as `YYYY_qN.nc`.
The maps made from this hindcast are in `best_two_step/hindcast_onestep/`.
