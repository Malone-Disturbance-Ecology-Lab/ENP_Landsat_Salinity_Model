# ENP_Landsat_Salinity_Model
A Landsat version of https://github.com/Malone-Disturbance-Ecology-Lab/ENP_Salinity_Model

## Script Explanations
- **00_download_dbhydro.R**: This script finds and downloads salinity data (continuous and grab measurements) from DBHydro from 04/01/2013 to 2/13/2026.
  
- **00_download_landsat.R**: This script creates a download request for the Harmonized Landsat and Sentinel-2 Land Surface Reflectance rasters (HLSL30.020) using the NASA AppEEARS API. The Landsat rasters will cover the Everglades National Park (ENP).
  
- **01_calculate_dist_coast.R**: This script calculates distance to the coast. NOTE: run the accompanying shell script **01_calculate_dist_coast.sh** on Grace cluster for fast computation.
  
- **01_download_elevation_slope.R**: This script downloads elevation and slope data for Everglades National Park using `elevatr::get_elev_raster()` and `terra::terrain()`.
  
- **01_download_meteorology.R**: This script downloads meteorology data for Everglades National Park using `climateR::getGridMET()`. NOTE: run on Grace cluster for fast computation.
  
- **01_harmonize_dbhydro.R**: This script harmonizes DBHydro salinity data (continuous and grab measurements) into a clean data frame.
  
- **01_harmonize_landsat.R**: This script harmonizes Landsat rasters into a clean tif file for each band.
  
- **02_reformat_meteorology.R**: This script reformats meteorology data for Everglades National Park by changing the resolution to match Landsat. NOTE: run the accompanying shell script **02_reformat_meteorology.sh** on Grace cluster for fast computation.
  
- **03_extract_and_harmonize_everything.R**: This script extracts and harmonizes these data sources into a clean data frame: Landsat, precipitation, solar radiation, temperature, elevation, slope, distance to coast, DBHydro salinity. 

- **04_model_validity_temporal_rf_timeframe.R**: This script checks for the validity of random forest models by calculating R-squared for observed vs predicted for every station. This tests how well does the model capture temporal patterns in the last year of the training set across different timeframe lengths.

- **05a_model_validity_rsquared_cutoff.R**: This script checks for the ideal R-squared cutoff value to be considered a "good" value in order to find out which stations the model consistently performs well for.

- **05b_model_validity_visualizations_rsq_cutoff.R**: This script plots the selected model variables for each station to see if there are any patterns for "good" or "bad" stations.

- **05c_model_validity_visualizations_check_pred.R**: This script checks if my predictions are over or underestimating and plots the difference between predicted & actual salinity value.

- **06_sensitivity_analysis.R**: This script is for performing sensitivity analysis on our model. It investigates how the model predicts when all other variables are held constant (aside from our variable of interest).