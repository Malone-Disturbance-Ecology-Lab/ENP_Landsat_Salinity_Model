# ENP_Landsat_Salinity_Model
A Landsat version of https://github.com/Malone-Disturbance-Ecology-Lab/ENP_Salinity_Model

Daily and monthly park-wide surface-water salinity maps for Everglades National Park (ENP), 2013–2025, on the native 30 m Landsat (HLS L30) grid, with a 250 m daily version. The maps are built from DBHydro station salinity, HLS Landsat vegetation spectra, EDEN water levels, water-management structure flows and gridMET climate.

IMPORTANT: Ctrl+F search for "CHANGE AS NEEDED" in each R script to look for sections where you can change and customize your own user settings. The Python scripts use the working-directory layout described in [Running the Python pipeline](#running-the-python-pipeline).

---

## Current model (September 2026): two-step salinity model

Salinity at a site is split into a slowly changing **site mean** and a **daily departure** from that mean. Each part is predicted from the data best suited to it.

**Step 1: time-varying site mean.**
- Target: B = centered 12-month mean salinity at each station and month.
- Main predictors: the vegetation spectral signature around the site, from HLS L30, 2013–2025.
  - Only clear pixels are used (Fmask cloud, cirrus, shadow, adjacent-to-cloud and snow bits removed).
  - Vegetated pixels are those with Landsat water frequency < 0.2.
  - Indices: NDVI, EVI2, NDMI and SWIR1. Each is summarised within 1 km and 3 km of the site as a centered 12-month composite, a 3-month composite, the long-term median, and the 12-month departure from that median.
- Other predictors:
  - static site features: signed distance to coast, vegetation-map fractions, along-channel distance to the coast, channel density, flow-basin weights and distance to each freshwater source;
  - slow hydro-climate: 365- and 180-day precipitation, the 90-day EDEN stage averaged over 5 km, and 180-day SRS and Taylor Slough/C-111 flows.
- Model: a gradient-boosted tree model (scikit-learn `HistGradientBoostingRegressor`) with very shallow trees (4 leaves, at least 400 samples per leaf). The predicted mean is smoothed with a 12-month rolling window.

**Step 2: daily departure from the mean (salinity − B).** Only daily-resolution drivers are used, so the predictions do not depend on how often Landsat sees a site:
- gridMET precipitation windows (30/90/180/365 days), solar radiation and temperature;
- EDEN daily water-surface stage at the nearest cell and averaged over 5 km (daily, 30-day and 90-day);
- structure flows: SRS, Taylor Slough/C-111, S-18C and S-197 (30/90/180 days);
- six anchor stations (ENPWP, ENPSR, ENPWW, ENPGB, TAYLORS3, ENPWB): 3-day and 30-day salinity;
- for each driver, its departure from its own centered 365-day mean;
- day of year and the static site features.

**Final prediction** = 0.7 × (step-1 mean + step-2 departure) + 0.3 × the one-step v3 model. The v3 model is a gradient-boosted model on the same site, connectivity, climate, EDEN, flow and anchor features. Predictions are clipped at 0.

**90% intervals**: split-conformal intervals calibrated by predicted-salinity level (0–5, 5–15, 15–30 and >30 psu), using out-of-block residuals from 2013–2024.

### Validation

All scores are for stations held out in spatial-block cross-validation. The 83–84 non-anchor stations are grouped into 10 spatial blocks with KMeans, and each block is predicted by a model trained without it. Models are trained on 2013–2024 and also scored on 2025, which is held out in both space and time.

| model | CV R² | CV RMSE (psu) | 2025 R² | 2025 RMSE (psu) | median within-station r |
|---|---|---|---|---|---|
| v3 one-step GBM | 0.786 | 6.23 | 0.701 | 6.74 | 0.855 |
| two-step only | 0.797 | 6.08 | 0.774 | 5.87 | 0.838 |
| **0.7 × two-step + 0.3 × v3 (final)** | **0.808** | **5.91** | **0.779** | **5.79** | **0.863** |
| two-step with the *observed* mean (upper bound) | 0.920 | 3.82 | – | – | 0.883 |

- **Vegetation spectra improve the site mean.** Step-1 R² for the 12-month mean is 0.84 in CV and 0.83 in 2025 with spectra, versus 0.79 in both without them.
- **Most remaining error is in the site mean at unmonitored sites.** Given the observed mean, the daily step alone reaches R² 0.92.

**Daily model vs. gauge** (held-out predictions at 83 stations):
- median within-station r is 0.86 (0.84 in 2025);
- median KGE is 0.61 (0.66 in 2025); 69% of stations have KGE > 0.5 and 94% have KGE > 0;
- median RMSE is 5.1 psu and median absolute bias 3.0 psu.

**90% intervals (2025)**: 90% coverage, with a median width of 20.6 psu.

**Known limitations**
- **Fresh interior marsh gauges are over-predicted.** Sites that read about 0.3 psu are predicted at 4–8 psu (TAYSLOUWET, EDEN3, 22908295, ENPCN). Log or square-root targets reduced this but lowered overall skill (2025 R² 0.69–0.78).
- **Wet-season lows are too high at some mangrove creeks** (ENPHR, ENPLO, ENPGI).
- **The northwest Gulf coast (Ten Thousand Islands) is predicted too fresh.** Few training stations are nearby, so treat that area as extrapolation.
- **The spectra were extracted with some compromises:** Landsat L30 only (no Sentinel-2 S30), rows sampled about 480 m apart, and no green or SWIR2 bands, so NDWI and NBR2 are unused. These choices kept reads from the server tractable. Filling these gaps is the most direct way to improve step 1.

### Outputs (server: `ENP_Landsat_Salinity_Model/best_two_step/`)
- `best_two_step_cv.csv`: CV and 2025 scores for every candidate combination.
- `final_interval_calibration.csv`: coverage and width for each interval-calibration scheme.
- `timeseries/`: per-station daily model vs. gauge.
  - `station_daily_model_vs_gauge.csv.gz`: held-out predictions.
  - `station_timeseries_metrics.csv`: r, RMSE, bias, NSE and KGE by station.
  - Figures: one per basin, plus a summary.
- `daily_maps/`:
  - `daily/YYYY/ENP_salinity_daily_YYYY_MM.tif`: one band per day; salinity in psu × 10 (int16, nodata −1).
  - `monthly/ENP_salinity_monthly_{mean,lower90,upper90}_YYYY.tif`: 12 bands per year.
  - `station_daily_model_series_2013_2025.csv.gz`: continuous daily model series at every station.
  - Preview figures.
- `scripts/`: copies of the scripts in `python/`.

Grid (250 m products): every 8th pixel of the 30 m HLS grid, EPSG:4326, covering the ENP Landsat footprint including Florida Bay.

### 30 m products (server: `best_two_step/salinity30m/`)
All inputs are resampled to the native HLS L30 30 m grid, so output pixels match Landsat pixels.
- **Grid:** 3839 × 4171 cells, 0.00027092271°, 7,622,366 cells in the domain.
- **Inputs:** site variables are native 30 m. Spectral composites (240 m row grid), EDEN (400 m) and gridMET (4 km) are bilinearly interpolated to each cell.

| folder or file | contents |
|---|---|
| `monthly/ENP_salinity30m_monthly_YYYY_MM.tif` | monthly mean salinity, 156 files, psu × 10 (int16, nodata −1) |
| `hpc/`, `hpc_inputs/` | Grace SLURM array job that writes the full-model **daily** 30 m archive (`YYYY/ENP_salinity30m_daily_YYYY_MM.tif`, about 11–15 GB); see `hpc/README_HPC.md` |
| `change30/` | 30 m region/ecosystem summaries, trends (annual means and monthly anomalies), high-salinity events, station checks |
| `figures/` | Figures 1–7 of the methods/results document (made by `make_figures.py`) |
| `ENP_region_ecosystem_classes_30m.tif` | region (band 1) and ecosystem (band 2) for each 30 m cell |
| `station_daily_model_series_2013_2025_250m_30m.csv.gz` | daily model at 90 stations (250 m and 30 m) with gauge values |
| `scripts/` | copies of the 30 m scripts |

**How the monthly maps are built.** The full model was run at every 30 m cell on the 15th of each month. Daily 30 m salinity is the interpolated 250 m daily map plus the 30 m detail (full 30 m model − 250 m map), interpolated in time. The monthly maps average it.

**Checks.**
- The 30 m and 250 m maps agree at the 250 m nodes (mean absolute difference ≤ 0.08 psu).
- A full-model daily test month on the HPC code matched the monthly map to 0.04 psu.
- Gauge skill is the same as at 250 m (in-sample daily R² 0.94).
- Within 250 m cells, 30 m salinity varies little: median SD 0.2 psu in mangroves, up to 1–2 psu along shorelines, creeks and vegetation edges.

---

## Where things are

### Data and outputs (MaloneLab server)
Root: `/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model` (`smb://corellia.environment.yale.edu/MaloneLab/Research/ENP_Landsat_Salinity_Model`)

| what | location |
|---|---|
| Station × date table (salinity + Landsat + site) | `DBSAL.csv` |
| Cleaned daily data and train/test split | `cleaned_data/` |
| HLS Landsat/Sentinel stacks (per band, all dates) | `harmonized_appeears_landsat_data/L30/`, `S30/` |
| DBHydro station locations (shapefile) | `ENP_DBHydro_sf/`, `DBHydro_lonlat.csv` |
| Structure flows (daily) | `dbhydro_flows/DBHydro_structure_flows_daily_2013_2026.csv` |
| EDEN daily water surfaces (52 quarterly NetCDFs) | `data/EDEN_GRIDS_Q/` |
| EDEN stage at stations | `eden_grid/` |
| Site variables (signed distance, water frequency, mainland mask, vegetation map, station table) | `site_variables/` |
| Vegetation spectral composites (HLS L30 rows) | `two_step_model/veg_rows_L30_part1-3.npz` |
| **Final model**: CV results, interval calibration, scripts | `best_two_step/` |
| **Gauge vs model time series** (held-out), per-station metrics, figures | `best_two_step/timeseries/` |
| **Daily salinity maps** 2013–2025 (monthly GeoTIFFs, one band per day, psu × 10) | `best_two_step/daily_maps/daily/YYYY/` |
| Monthly mean and 90% bound maps | `best_two_step/daily_maps/monthly/` |
| Daily modeled series at 90 stations | `best_two_step/daily_maps/station_daily_model_series_2013_2025.csv.gz` |
| **30 m salinity maps** (monthly GeoTIFFs), 30 m class raster, 30 m station series | `best_two_step/salinity30m/` |
| **Daily 30 m archive** (written by the HPC job) and its inputs | `best_two_step/salinity30m/hpc/`, `hpc_inputs/` (output: `best_two_step/daily30m/`) |
| Region/ecosystem summaries, trends, high-salinity events (30 m, current) | `best_two_step/salinity30m/change30/` |
| Figures 1–7 of the methods/results document | `best_two_step/salinity30m/figures/` |
| Same analyses at 250 m (previous version) | `best_two_step/daily_maps/change_analysis/` |
| Superseded model versions (v1–v3, connectivity and basin CV) | `spatial_cv/`, `connectivity_uncertainty/`, `basins_monthly/`, `basins_monthly_v3/`, `two_step_model/` |
| Legacy RF model outputs | `model_validity_results_timeframe_diff/`, `deprecated/` |

gridMET daily NetCDFs (`agg_met_{pr,srad,tmmn,tmmx}_1979_CurrentYear_CONUS.nc`) were downloaded to a local machine and are not yet on the server; copy them to `data/gridMET/` before re-running `gridmet_prep.py`.

### Documents (claude.ai; private to the owner until shared)
| document | link |
|---|---|
| Methods and results write-up (model, validation, 30 m mapping, salinity across regions/ecosystems, change 2013–2025, high-salinity events; Figures 1–7) | https://claude.ai/code/artifact/98493237-64e4-4392-a3a6-7b893b116bc3 |
| EDI data-package metadata (for publishing the daily 30 m layers) | https://claude.ai/code/artifact/62565e15-a16b-4721-ab45-4ed2dd245386 |
| Working notes for each analysis step (dated `claude/*.md` docs) | Claude project "ENP_Landsat_Salinity_Model" |

---

## Repository layout

| path | contents |
|---|---|
| `00_*.R`, `01_harmonize_*.R`, `03_*.R`, `03b_clean_and_split_data.R` | Data acquisition and harmonization (R): DBHydro salinity, HLS Landsat/Sentinel, the station × date table `DBSAL.csv`, then cleaning and the train/test split |
| `python/` | Site variables, drivers, connectivity/basins, the two-step model, validation and daily maps (Python) |
| `appeears_landsat_lookup/` | HLS Fmask lookup tables |
| `deprecated/` | Superseded models, analyses and data-prep scripts, kept for reference (see [Deprecated](#deprecated)) |

### R scripts (data acquisition and cleaning)
- **00_download_dbhydro.R**: finds and downloads DBHydro salinity (continuous and grab) from 04/01/2013 to 2/13/2026. Outputs: `sfwmd-data-2026XXXXXXXX.zip`.
- **00_download_landsat.R**: AppEEARS request for HLS L30/S30 surface reflectance and Fmask over ENP. Outputs: `HLS{L,S}30.020_BXX_doyXXXXXXX_aid0001_17N.tif`.
- **01_harmonize_dbhydro.R**: DBHydro salinity to a clean data frame. Outputs: `ENP_DBHydro_salinity_df.csv`, `DBHydro_lonlat.csv`, `ENP_DBHydro_sf.shp`.
- **01_harmonize_landsat.R**: one multi-date stack per band. Outputs: `L30_ENP_B01…B07.tif`, `S30_ENP_B01…B04, B8A, B11, B12.tif`, `L30_ENP_Fmask.tif`, `S30_ENP_Fmask.tif`.
- **03_extract_and_harmonize_everything.R**: station × date table of salinity, Landsat bands and indices, meteorology and site variables. Output: `DBSAL.csv`.
  - Some of its inputs come from scripts now in `deprecated/superseded_data_prep/`. Those outputs are already on the server; the current model replaces them with the Python site variables and native gridMET.
- **03b_clean_and_split_data.R**:
  - removes salinity values > 80 psu or < 0;
  - masks Fmask cloud, cirrus, shadow, adjacent and snow pixels, and out-of-range reflectance;
  - averages to one row per station per day;
  - splits the data into train (≤ 2024) and test (2025).
  - Outputs in `cleaned_data/`: `DBSAL_daily_clean.csv`, `DBSAL_train_2013_2024.csv`, `DBSAL_test_2025.csv`, `cleaning_summary.csv`, `removed_salinity_values.csv`.

### Running the Python pipeline
Requirements: Python ≥ 3.10 with numpy, pandas, scipy, scikit-learn, scikit-image, rasterio, pyproj, pyshp and matplotlib.

All scripts in `python/` run from one working directory (e.g. `python/` itself). They read and write intermediate files (`*.pkl`, `*.npz`, `*.csv`) there. Create these three entries in that directory:

| entry | points to |
|---|---|
| `SERVER/` | the server project folder, e.g. `ln -s /Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model SERVER` |
| `sitevars/` | `SERVER/site_variables` |
| `downloads/` | a fast local folder holding the gridMET NetCDFs (`agg_met_{pr,srad,tmmn,tmmx}_1979_CurrentYear_CONUS.nc`) and the EDEN intermediates |

The EDEN quarterly NetCDFs (`YYYY_qN.nc`, 2013–2025) are read from `SERVER/data/EDEN_GRIDS_Q/`. They are available at `https://sflthredds.er.usgs.gov/thredds/fileServer/eden/surfaces/YYYY_qN.nc`.

Run order:

**1. Site variables** (30 m HLS grid; outputs to `SERVER/site_variables/`)
1. `01_water_mask_and_signed_distcoast.py`: Landsat water frequency, mainland mask, and signed distance to the open coast (negative offshore).
2. `02_rasterize_vegmap.py`: rasterizes the NPS EVER vegetation map.
3. `03_extract_station_site_variables.py`: station distance and vegetation-class fractions. Writes `station_site_variables.csv` and `stations_coast.csv`.

**2. Drivers**
4. `build_flows_daily.py`: daily DBHydro structure flows (S-12A–D + S-333, S-332D/S-174/S-175, S-18C, S-197).
5. `extract_eden_grid_at_stations.py`: daily EDEN stage at each station's nearest valid cell and averaged over 5 km.
6. `eden_monthly_grids.py`: EDEN 30- and 90-day mean grids at mid-month. `eden_mean_grids.py`: 2025 dry- and wet-season grids used to trace basins.
7. `gridmet_prep.py`: gridMET P30/P90/P180/P365, R30 and T90 on the native 4 km grid.
8. `veg_spectra_rows_L30.py <worker> <n_workers>`: monthly vegetation-index medians from `L30_ENP_*.tif`, averaged in 240 m blocks along every 16th row. Run the workers in parallel; the script resumes if interrupted.
9. `pack_intermediates.py`: packs steps 5, 6 and 8 into single files.
10. `build_station_feature_table.py`: station-day feature table (`all_features.pkl`).

**3. Connectivity and flow basins**
11. `connectivity.py`: along-channel distance to the coast, channel density, and distance to each freshwater source.
12. `conn_cv.py`: drainage regions and station connectivity table (also runs the connectivity CV).
13. `basins.py`: flow basins traced up the EDEN water surface from inflow seed zones, then assigned along water by least cost, with soft basin weights.
14. `smooth_offshore_v3.py`: smooths the offshore basin, cost and source fields (4.3 km Gaussian, ramped in over 3 km). This removes wedge artifacts in Florida Bay.

**4. Model**
15. `prep_v3_tables.py`: builds the v3 station tables. Station climate comes from harmonically filled native gridMET, and offshore EDEN stage from harmonic extrapolation. Uses `fillers.py`, `spatial_cv.py` and `conn_cv.py`.
16. `two_stage.py`: vegetation-spectra features at stations and the step-1/step-2 tables (`two_stage_tables.pkl`), plus the first round of two-step CV. `two_stage_v2.py` is the second round (optional).
17. `s1_tune3.py`: step-1 tuning; shallow trees generalize best. Uses `ts_setup.py`.
18. `best_two_step.py`: spatial-block CV of the candidate combinations. Saves held-out daily predictions.
19. `ts_compare.py`: daily model vs. gauge time series, metrics and figures.
20. `fit_final.py`: final models on all 2013–2025 data, plus interval calibration (`final_models.pkl`).

**5. Daily maps**
21. `eden_daily.py`: daily EDEN fields (cell, 30-day, 5 km 30/90-day, centered 365-day means). Extended offshore by harmonic extrapolation and cropped to the map domain.
22. `daily_maps.py`: daily 250 m predictions for 2013–2025, about 2.5 s per day on 2 cores. The script resumes from `daily_maps/progress.txt`.
23. `export_daily.py`: monthly GeoTIFFs of the daily maps, monthly mean and 90% bounds, station series and preview figures.

**6. 30 m maps and analyses**

24. `prep30.py`: builds the 30 m static predictors, coordinates and the spectral sampling grid (`static30.npy`, `spec_grid.npy`, …).
25. `hpc/compile_models.py`: compiles the fitted models with treelite/tl2cgen. Predictions are identical and about 12× faster.
26. `run30.py`: runs the full model at every 30 m cell on the 15th of each month and stores the 30 m detail layer (`map30/C_int16.npy`). It takes about 3 h on 2 cores and resumes by row tile.
27. `export30.py`: writes the monthly 30 m GeoTIFFs.
28. `classes30.py`, `analyze30.py`: 30 m region and ecosystem classes, summaries, and Theil–Sen / Mann–Kendall trends with BH-FDR.
29. `events30.py`: high-salinity events from daily 30 m salinity, and station series. Takes about 1.5 h.
30. `station_check30.py`: gauge checks (daily skill, annual trends, event days).
31. `trend30_seasonal.py`, `station_trend_monthly30.py`: trends in monthly anomalies at every 30 m pixel and at the gauges. Uses the seasonal Kendall test with the Hirsch–Slack correction and the seasonal Theil–Sen slope, with BH-FDR across pixels. Takes about 10 min.
32. **`make_figures.py`**: makes every figure in the methods/results document (Figures 1–7) from the saved outputs. Run `python make_figures.py` for all of them, or `python make_figures.py 4 5` for a subset. Wherever colour encodes salinity, red = more saline.
33. `hpc/daily30_hpc.py` + `hpc/run_daily30.slurm`: the full-model daily 30 m archive on Grace, one array task per month (see `hpc/README_HPC.md`).

---

## Deprecated
Kept in `deprecated/` for reference; they are not part of the current model.
- `deprecated/rf_model/` (`04`–`07`): random-forest temporal model validity, R² cutoffs, visualizations and sensitivity analysis.
  - That model was evaluated with 2025 inside its training windows, so its skill was partly in-sample.
  - It did not beat a station × month climatology.
- `deprecated/exploratory/` (`99_*`): early RF, XGBoost, neural-net and VSURF exploration, and spatial and temporal validity tests.
- `deprecated/superseded_data_prep/`: scripts replaced by newer data sources.
  - `01_calculate_dist_coast`: unsigned distance to any shoreline. Replaced by the signed distance to the open coast.
  - `01_download_elevation_slope`: elevatr elevation, unreliable over water and mangrove (5–8 m at coastal stations). Not used.
  - `01_download_meteorology` and `02_reformat_meteorology`: gridMET resampled to 30 m with terra. Replaced by native 4 km gridMET, which carries the same information without about 130 GB of rasters.
- `deprecated/band_and_target_screening/`: band-screening and salinity-target analyses (`03b_build_salinity_targets.R`, `08_screen_salinity_bands.R`, results notes and tests). They showed that Landsat bands add little *temporal* skill at a site; that finding led to using vegetation spectra for the site mean instead.

Intermediate Python model versions (v1–v3 one-step GBMs, connectivity and basin CV, the conformal-by-zone maps) are not in this repository. Their results remain on the server in their original folders (`spatial_cv/`, `connectivity_uncertainty/`, `basins_monthly*/`, `two_step_model/`).
