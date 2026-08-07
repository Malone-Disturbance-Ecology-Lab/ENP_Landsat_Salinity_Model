# ENP_Landsat_Salinity_Model
A Landsat version of https://github.com/Malone-Disturbance-Ecology-Lab/ENP_Salinity_Model

IMPORTANT: Ctrl+F search for "CHANGE AS NEEDED" in each script to look for sections where you can change and customize your own user settings. 

## Script Explanations
- **00_download_dbhydro.R**: This script finds and downloads salinity data (continuous and grab measurements) from DBHydro from 04/01/2013 to 2/13/2026. Outputs:
   - sfwmd-data-2026XXXXXXXX.zip files for continuous and grab data
  
- **00_download_landsat.R**: This script creates a download request for the Harmonized Landsat and Sentinel-2 Land Surface Reflectance rasters (HLSL30.020) using the NASA AppEEARS API. The Landsat rasters will cover the Everglades National Park (ENP). Outputs:
   - HLSL30.020_BXX_doyXXXXXXX_aid0001_17N.tif files for each band
   - HLSS30.020_BXX_doyXXXXXXX_aid0001_17N.tif files for each band
   - HLSL30.020_Fmask_doyXXXXXXX_aid0001_17N.tif files for each band
   - HLSS30.020_Fmask_doyXXXXXXX_aid0001_17N.tif files for each band
  
- **01_calculate_dist_coast.R**: This script calculates distance to the coast. NOTE: run the accompanying shell script **01_calculate_dist_coast.sh** on Grace cluster for fast computation. Outputs:
   - ENP_DistCoast.tif 
  
- **01_download_elevation_slope.R**: This script downloads elevation and slope data for Everglades National Park using `elevatr::get_elev_raster()` and `terra::terrain()`. Outputs:
   - ENP_Elevation.tif
   - ENP_Slope.tif 
  
- **01_download_meteorology.R**: This script downloads meteorology data for Everglades National Park using `climateR::getGridMET()`. NOTE: run on Grace cluster for fast computation. Outputs:
   - ENP_Precipitation.tif
   - ENP_SolarRadiation.tif
   - ENP_AverageTemperature.tif 
  
- **01_harmonize_dbhydro.R**: This script harmonizes DBHydro salinity data (continuous and grab measurements) into a clean data frame. Outputs:
   - ENP_DBHydro_salinity_df.csv
   - DBHydro_lonlat.csv
   - ENP_DBHydro_sf.shp 
  
- **01_harmonize_landsat.R**: This script harmonizes Landsat rasters into a clean tif file for each band. Outputs:
   - L30_ENP_B01.tif, L30_ENP_B02.tif, L30_ENP_B03.tif, L30_ENP_B04.tif, L30_ENP_B05.tif, L30_ENP_B06.tif, L30_ENP_B07.tif
   - S30_ENP_B01.tif, S30_ENP_B02.tif, S30_ENP_B03.tif, S30_ENP_B04.tif, S30_ENP_B8A.tif, S30_ENP_B11.tif, S30_ENP_B12.tif
   - L30_ENP_Fmask.tif, S30_ENP_Fmask.tif
  
- **02_reformat_meteorology.R**: This script reformats meteorology data for Everglades National Park by changing the resolution to match Landsat. NOTE: run the accompanying shell script **02_reformat_meteorology.sh** on Grace cluster for fast computation. Outputs:
   - ENP_AverageTemperature_XXXX_XXXX.tif files
   - ENP_Precipitation_XXXX_XXXX.tif files
   - ENP_SolarRadiation_XXXX_XXXX.tif files
  
- **03_extract_and_harmonize_everything.R**: This script extracts and harmonizes these data sources into a clean data frame: Landsat, precipitation, solar radiation, temperature, elevation, slope, distance to coast, DBHydro salinity. Outputs:
   - DBSAL.csv  

- **04_model_validity_temporal_rf_timeframe.R**: This script checks for the validity of random forest models by calculating R-squared for observed vs predicted for every station. This tests how well does the model capture temporal patterns in the last year of the training set across different timeframe lengths. Outputs:
   - obs_vs_pred_XXXX_XXXX_full_results.csv files
   - obs_vs_pred_XXXX_XXXX_summary.csv files

- **05a_model_validity_rsquared_cutoff.R**: This script checks for the ideal R-squared cutoff value to be considered a "good" value in order to find out which stations the model consistently performs well for. Outputs:
   - results_0.33_cutoff.csv 

- **05b_model_validity_visualizations_rsq_cutoff.R**: This script plots the selected model variables for each station to see if there are any patterns for "good" or "bad" stations. Outputs:
   - {VARIABLE}_plot.png for all model variables for both R-squared 0.33 and 0.5 cutoffs

- **05c_model_validity_visualizations_check_pred.R**: This script checks if my predictions are over or underestimating and plots the difference between predicted & actual salinity value. Outputs:
   - pred_year_vs_rsq_positive.png
   - pred_year_vs_rsq.png
   - pred_year_vs_sal.png

- **05d_model_validity_visualizations_explore.R**: This script plots:
  - A map of good stations
  - Salinity bins (diff_sal_bins_XXXX.png files)
  - Average predictions vs prediction year, colored by good or bad stations (pred_vs_pred_year.png)
  - Good vs. bad, years and stations
  - Salinity vs. average daily prediction for all prediction years (sal_vs_avgpred_XXXX.png files)

- **06_sensitivity_analysis.R**: This script is for performing sensitivity analysis on our model. It investigates how the model predicts when all other variables are held constant (aside from our variable of interest). Outputs:
   - {VARIABLE}_{LEVEL}_obs_vs_pred_XXXX_XXXX_full_results.csv for all model variables and percentile levels
   - {VARIABLE}_{LEVEL}_obs_vs_pred_XXXX_XXXX_summary.csv
