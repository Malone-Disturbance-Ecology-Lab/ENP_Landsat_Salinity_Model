# ENP_Landsat_Salinity_Model

Daily surface-water salinity for Everglades National Park (Florida, USA), April 2000 – December 2025, on the 30 m Landsat grid, with a daily surface-water mask and a 250 m version carrying calibrated 90 % prediction intervals. A Landsat-grid successor to https://github.com/Malone-Disturbance-Ecology-Lab/ENP_Salinity_Model.

**Manuscript:** *Daily 30 m surface-water salinity for Everglades National Park, 2000–2025, from gauge records and hydrologic drivers* (in preparation for *Earth System Science Data*). Working draft: [ESSD manuscript (Google Doc)](https://docs.google.com/document/d/1XrcMQsavZ3_EjXpiEht6dzYCoFI50tBKd_fPUGjt-tA/edit). Access is by request until the paper is submitted.

**Data:** Environmental Data Initiative, DOI to be added on release (CC BY 4.0).

This repository holds only the workflow that produces the manuscript and its data products. Earlier models and analyses are in [`deprecated/`](#deprecated).

---

## The model

One gradient-boosting model (scikit-learn `HistGradientBoostingRegressor`) predicts daily salinity directly.

- **Training data:** 84 DBHydro gauges (41 continuous sondes, 43 grab-sample stations) plus 13 freshwater sites with DBHydro field specific conductance converted to practical salinity (TEOS-10). In total, 355 156 station-days, April 2000 – December 2025, with equal weights.
- **56 predictors.** Static predictors (27) are on the 30 m grid: signed distance to the coast, vegetation composition within 250 m and 1 km, hydrologic connectivity, and least-cost distances and weights for the freshwater source basins. Daily drivers:
  - gridMET precipitation (30/90/180/365 days), 30-day shortwave radiation and 90-day temperature;
  - EDEN v3 water level at the cell (30 days) and within 5 km (30 and 90 days);
  - water-control structure flows (SRS, Taylor Slough/C-111, basin-weighted source flow);
  - 3- and 30-day salinity at six anchor stations (ENPWP, ENPSR, ENPWW, ENPGB, TAYLORS3, ENPWB);
  - NOAA monthly mean sea level (Key West, Vaca Key, Virginia Key; 1-, 3- and 12-month means);
  - day-of-year sine and cosine.
- **Hyperparameters:** learning rate 0.05, at most 400 iterations, 31 leaves, at least 200 samples per leaf, L2 1.0, `random_state` 0; scikit-learn 1.8.
- **Validation:** 10 spatial blocks (k-means on station coordinates), each predicted by a model trained without it. Held-out R² 0.85, RMSE 5.3 psu, bias −0.1 psu over 352 714 station-days; R² 0.82–0.86 in every period.
- **Intervals:** split-conformal 90 % half-widths by predicted level (0–5, 5–15, 15–30, > 30 psu): 4.3, 10.3, 10.0, 5.9 psu. Coverage is 89.5 % when calibrated on nine blocks and tested on the tenth.
- **Surface-water mask:** daily wet/dry codes from EDEN water levels and the EDEN elevation models with uncertainty bounds (Simons et al., 2024). Salinity is set to missing on dry (0) and never-wet (9) days. The mask agrees with EDEN gauge dry flags on 94.5 % of committed days.

## Data products

All products are on the MaloneLab server under `ENP_Landsat_Salinity_Model/best_two_step/` and will be released on EDI. They share one grid: EPSG:4326, 3839 × 4171 cells of 0.00027092271° (upper left −81.51875°, 25.88965°), with a 250 m subgrid of every 8th cell. Salinity is stored as int16 psu × 10 (divide by 10), with −1 as the declared nodata value.

| product | server folder | files |
|---|---|---|
| **Daily salinity, masked (released record)** | `daily30m_masked_final/` | `YYYY/ENP_salinity30m_masked_daily_YYYY_MM.tif`, 309 files, one band per day, 12.3 GB |
| Daily surface-water mask | `wetdry30m/` | `YYYY/ENP_wetdry30m_daily_YYYY_MM.tif`, uint8 codes 0 dry, 1 wet, 2 uncertain, 3 permanent water, 4 tidal coastal wetland, 5 not assessed, 9 never wet |
| Wet-day monthly salinity | `monthly30m_wetdays/` | band 1 mean over days with surface water; band 2 number of those days |
| Climatology | `climatology_final/` | wet-day mean, dry and wet season, May, October, interannual SD, share of wet days, 95th percentile; annual means |
| Daily salinity at 250 m with 90 % bounds | `onestep_final_products/daily250/`, `monthly250/` | daily maps; monthly means of salinity and its bounds |
| Trends and high-salinity events | `long_masked_final/results/` | per-pixel seasonal Kendall / Theil–Sen trends with BH-FDR, event days, region and ecosystem tables |
| Unmasked daily salinity (every cell, every day) | `daily30m_onestep_final/` | model output before masking |
| Fitted model and CV predictions | `onestep_final/` | `final_model_onestep.pkl`, `m3_final.so`, `cv_predictions_final.parquet`, `q_levels.json`, `msl_by_date.csv` |
| Paper figures and tables | `paper_final/` | Figs 1–10 and the statistics in the manuscript |

Reading the files: GDAL-based tools (QGIS, R `terra`, Python `rasterio`/`rioxarray`) read −1 as missing. For example, `terra::rast(f) / 10`. In the masked files, −1 marks both cells outside the domain and days without surface water; use the mask files to tell them apart or to drop uncertain days (code 2). See [`surface_water_mask/apply_mask/README.md`](surface_water_mask/apply_mask/README.md).

---

## Repository layout

| path | step | contents |
|---|---|---|
| `00_*.R`, `01_harmonize_*.R`, `03_extract_and_harmonize_everything.R`, `03b_clean_and_split_data.R` | 1 | DBHydro salinity and HLS Landsat download and harmonisation (R), the station × date table and its cleaning (2013–2025) |
| `hindcast/` | 2 | DBHydro salinity, structure flows and gridMET before 2013; EDEN weights; the station feature table for 2000–2025 |
| `python/` | 3 | 30 m site variables, drivers, connectivity and flow basins, station tables, 30 m static layers, region and ecosystem classes |
| `final_model/` | 4 | freshwater-site features, the model-fitting harness, final fit with spatial-block CV and conformal intervals, gauge events |
| `final_model/hpc/`, `hpc/` | 5 | Bouchet (Yale) SLURM jobs: compile the model, daily 30 m maps, 250 m and monthly products, climatology, masked trends and events |
| `surface_water_mask/` | 6 | daily wet/dry mask, gauge validation, and the script that writes the masked daily files |
| `paper/` | 7 | figure, statistics and model-driver scripts for the manuscript |
| `docs/` | – | README of the derived products; project notes on the manuscript |
| `appeears_landsat_lookup/` | – | HLS Fmask lookup tables |
| `deprecated/` | – | superseded models and analyses |

### Running the workflow

Requirements: R (data download and harmonisation); Python ≥ 3.10 with numpy, pandas, scipy, scikit-learn 1.8, rasterio, pyproj, pyshp, matplotlib, pyarrow, and treelite/tl2cgen to compile the model.

**Paths.** Scripts read the server project folder as `/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model` (the plotting scripts accept `ENP_DATA_ROOT`). Intermediate files are written to a local working folder, written in the scripts as `work/`:

| `work/` subfolder | holds |
|---|---|
| `work/map/` | working folder for `fresh_features.py` |
| `work/exp/` | the station tables built by `final_model/harness.py` and `fresh_features.py` |
| `work/final/` | the fitted model and CV predictions (copied to `best_two_step/onestep_final/`) |
| `work/drivers/`, `work/mask/`, `work/trend/` | intermediate files of the driver analysis, the mask and the masked trends |

The Python scripts in `python/` run from one working directory with three links: `SERVER/` → the server folder, `sitevars/` → `SERVER/site_variables`, and `downloads/` → a local folder holding the gridMET NetCDFs (`agg_met_{pr,srad,tmmn,tmmx}_1979_CurrentYear_CONUS.nc`) and the EDEN intermediates. EDEN v3 quarterly NetCDFs (`YYYY_qN.nc`, 2000–2025) are read from `SERVER/data/EDEN_GRIDS_Q/`. They are available from `https://sflthredds.er.usgs.gov/thredds/fileServer/eden/surfaces/YYYY_qN.nc`.

**1. Gauges and Landsat (R)**
- `00_download_dbhydro.R`, `01_harmonize_dbhydro.R`: DBHydro salinity (continuous and grab) for 2013–2026 and station locations.
- `00_download_landsat.R`, `01_harmonize_landsat.R`: HLS L30/S30 stacks. The final model uses only the Landsat water frequency (from Fmask) and the 30 m grid; vegetation spectra are not predictors.
- `03_extract_and_harmonize_everything.R`, `03b_clean_and_split_data.R`: the station × date table `DBSAL.csv` and `cleaned_data/DBSAL_daily_clean.csv` (0–80 psu, one row per station per day).

**2. Before 2013** (`hindcast/`; outputs in `SERVER/dbhydro_hindcast/`)
- `download_dbhydro_pre2013.py`, `harmonize_pre2013.py`: DBHydro salinity back to 1980, processed the same way.
- `download_flows_pre2013.py`, `build_flows_pre2013.py`: structure flows from 1998.
- `gridmet_pre2013.py`: filled gridMET windows from 1998.
- `eden_weights.py`: EDEN domain mask and harmonic-fill weights at the stations.
- `station_features_hindcast.py`: model inputs and observed salinity at all stations, 2000–2025 (`station_features_2000_2025_v2.parquet`). `check_features.py` compares them with the 2013–2025 tables.

**3. Site variables and drivers** (`python/`)
1. `01_water_mask_and_signed_distcoast.py`, `02_rasterize_vegmap.py`, `03_extract_station_site_variables.py`: Landsat water frequency, signed distance to the open coast, the rasterised NPS vegetation map, and station site variables.
2. `build_flows_daily.py`, `extract_eden_grid_at_stations.py`, `eden_monthly_grids.py`, `eden_mean_grids.py`, `gridmet_prep.py`, `eden_daily.py`: daily structure flows, EDEN stage at stations and on the grid, and gridMET windows (2013–2025).
3. `connectivity.py`, `conn_cv.py`, `basins.py`, `smooth_offshore_v3.py`: along-channel distance, channel density, flow basins traced up the EDEN water surface, least-cost distances, soft basin weights, and offshore smoothing.
4. `build_station_feature_table.py`, `prep_v3_tables.py` (with `fillers.py`, `spatial_cv.py`): station-day feature tables. `veg_spectra_rows_L30.py` and `pack_intermediates.py` are needed only because `prep30.py` reads their outputs.
5. `prep30.py`: the 30 m static predictors and coordinates (`static30.npy`, `px30_index.npy`, …). `classes30.py`: region and ecosystem class of each 30 m cell. `hpc/make_hpc_inputs.py`, `hpc/make_grid30.py`: the remaining HPC inputs.

**4. Model** (`final_model/`)
- Sea level: monthly mean sea level for NOAA stations 8724580, 8723970 and 8723214, downloaded from https://tidesandcurrents.noaa.gov to `SERVER/sea_level/noaa_monthly_mean_<id>.csv`.
- Freshwater sites: DBHydro field specific conductance for project EVER (test 9, surface samples ≤ 1 m), converted with TEOS-10; site summary in `SERVER/model_experiments_2026-10-13/dbhydro_EVER_conductivity_station_summary.csv`. `fresh_features.py` builds the model inputs at the 13 sites with the same code as the maps (`onestep_maps.py`).
- `harness.py`: loads the 2000–2025 station table, adds sea level and the spatial blocks, and runs the block CV.
- `fit_final.py`: block CV of the final configuration, conformal half-widths (`q_levels.json`), the final fit on all station-days, and the compiled model.
- `gauge_events.py`, `gauge_trend_long_fn.py`: trends and high-salinity events at the gauges for the checks in Sect. 4.4.
- `onestep_maps.py` (53 predictors, used by `fresh_features.py`) and `onestep_maps_final.py` (final model with sea level, used by the driver attribution) compute the model inputs at any 30 m cell and day. They read the HPC inputs from `final_model/in/` (a copy of `best_two_step/salinity30m/hpc_inputs/` plus `final_model_onestep.pkl`, `m3_final.so` and `msl_by_date.csv`). `fillers.py` is a copy of `python/fillers.py`.

**5. Maps on Bouchet** (`final_model/hpc/`; run from `~/enp_salinity/hpc`; `hpc/setup_bouchet.slurm` builds the `sal30` environment and the daily EDEN fields once)
- `submit_final.sh`: `compile_final.slurm`, then `run_daily30_onestep.slurm` (309 monthly tasks of `daily30_onestep_hpc.py`) → `daily30m_onestep_final/`.
- `submit_derive_final.sh`: `derive_final.py`, `reduce_monthly250.py` → monthly 30 m means, daily 250 m maps, and monthly 250 m means with 90 % bounds.
- `submit_clim_final.sh`: `wetdays_final.py`, `climatology_final.py` → wet-day monthly means and the climatology.
- `submit_long_final.sh`: `long_masked_hpc.py` (75 tiles), `long_masked_reduce.py` → masked trends, high-salinity events, region and ecosystem tables, station series.

**6. Surface-water mask** (`surface_water_mask/`)
- `scripts/mask_static.py`: static classes and 30 m DEM bounds. `scripts/wetdry30_hpc.py` with `run_wetdry30.slurm`: daily codes, 309 monthly tasks. `scripts/validate_gauges.py`: comparison with EDEN gauge dry flags.
- `apply_mask/mask_daily30.py` (`run_batch.sh` runs it in parallel): writes the released masked daily files.

**7. Manuscript** (`paper/scripts/`)
- `validation_final.py`: skill by station, region and period, and interval coverage.
- Figures: `fig01_study_area.py` (uses helpers in `python/make_figures.py`), `fig_workflow.py`, `fig_validation.py`, `fig_stations.py`, `fig_surface_water_mask.py`, `fig_spatial.py`, `fig_temporal.py`, `fig_trends.py`, `fig_events.py`, `fig_drivers.py`.
- Drivers (Sect. 5.6): `drivers_perm_importance.py`, `drivers_perm_combined.py` (grouped permutation importance with blocks held out), `drivers_group_response.py`, `drivers_partial_dependence.py`, `drivers_attrib_nodes.py` (variance attribution at 500 m nodes).

---

## Data sources

DBHydro (South Florida Water Management District); EDEN v3 water surfaces and EDEN DEMs with elevation uncertainty, v2 (USGS; Simons et al., 2024, doi:10.5066/P1XTHZ7D); gridMET (Abatzoglou, 2013); NOAA Tides and Currents monthly mean sea level; Harmonized Landsat Sentinel-2 L30 (NASA); NPS vegetation map of Everglades National Park (2021).

## Related documents

| document | link |
|---|---|
| ESSD manuscript (working version) | [Google Doc](https://docs.google.com/document/d/1XrcMQsavZ3_EjXpiEht6dzYCoFI50tBKd_fPUGjt-tA/edit) |
| Earlier ESSD draft (Claude doc) and project notes | Claude project "ENP_Landsat_Salinity_Model" |

---

## Deprecated

Kept in `deprecated/` for reference; not used for the manuscript or the released products.
- `two_step_hls_model/`: the September 2026 two-step model (site mean from HLS vegetation spectra plus a daily departure, blended with a one-step model) for 2013–2025, its 250 m and 30 m maps, analyses and HPC job. It was replaced by the single model above, which covers 2000–2025 without spectra.
- `released_onestep_2026-10-10/`: the earlier one-step model's daily 30 m HPC job (2013–2025 training), replaced by `final_model/hpc/`.
- `model_experiments_2026-10-13/`: the experiments that set the final configuration (training years, freshwater sites, sea level, weights).
- `hindcast_released_model/`: hindcast evaluation of the earlier one-step model before 2013.
- `rf_model/`, `exploratory/`: random-forest and other early models.
- `superseded_data_prep/`: data-prep scripts replaced by newer sources (unsigned coast distance, elevatr elevation, gridMET resampled to 30 m).
- `band_and_target_screening/`: Landsat band and salinity-target screening.
