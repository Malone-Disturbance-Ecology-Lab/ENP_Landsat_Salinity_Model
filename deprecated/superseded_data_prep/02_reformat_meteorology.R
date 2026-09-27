## --------------------------------------------- ##
#              Reformat Meteorology
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script reformats meteorology data for Everglades National Park
## by changing the resolution to match Landsat.
## NOTE: run on Bouchet cluster for fast computation.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(tidyverse)
library(terra)
library(sf)
library(climateR)

# CHANGE AS NEEDED --------------------------------

# Point to the Landsat Salinity folder on cluster
my_folder <- '/home/ac3656/GitHub/ENP_Landsat_Salinity_Model'

# Point to the project folder on cluster
my_project <- '/nfs/roberts/project/pi_sm3466/ac3656'

# -------------------------------------------------

# Create new folders to store rasters
dir.create(path = file.path(my_project, "ENP_Precipitation_Landsat_res"), showWarnings = F)
dir.create(path = file.path(my_project, "ENP_SolarRadiation_Landsat_res"), showWarnings = F)
dir.create(path = file.path(my_project, "ENP_AverageTemperature_Landsat_res"), showWarnings = F)

# Read in one raster to use as a template 
template_raster <- terra::rast(file.path(my_folder, "appeears_landsat_data", "L30", "B01", "HLSL30.020_B01_doy2013111_aid0001_17N.tif"))

# Read in rasters
#precip_r <- terra::rast(file.path(my_folder, "ENP_Precipitation.tif"))
#srad_r <- terra::rast(file.path(my_folder, "ENP_SolarRadiation.tif"))
tavg_r <- terra::rast(file.path(my_folder, "ENP_AverageTemperature.tif"))

## --------------------------------------------- ##
#             Change Resolution -----
## --------------------------------------------- ##

# Change resolution for precip, solar radiation, average temp as needed
data_type <- tavg_r
target_folder <- file.path(my_project, "ENP_AverageTemperature_Landsat_res")
target_name <- "ENP_AverageTemperature_"

# For 1 through 48...
for (i in 1:48){
  
  # If i is 48...
  if(i == 48){
    # Subset the last 2 rasters
    sub <- data_type[[4701:nlyr(data_type)]]
    # Change the resolution to match Landsat
    changed_res <- terra::resample(sub, template_raster)
    # Export
    terra::writeRaster(changed_res, file.path(target_folder, paste0(target_name, 4701, "_", nlyr(data_type), ".tif")), overwrite = T)
    # Stop for loop
    break()
  }
  
  # Set the increase
  increase <- (i-1) * 99
  
  # Starting number is i + increase
  num_1 <- i + increase
  # Ending number is i * 100
  num_2 <- i * 100
  
  message(num_1, ",", num_2)
  
  # Subset the data
  sub <- data_type[[num_1:num_2]]
  # Change the resolution to match Landsat
  changed_res <- terra::resample(sub, template_raster)
  
  # If num_1 is not a 4-digit number, add zeroes to the front to pad it
  if (num_1 < 10){
    num_1_name <- paste0("000", num_1)
  } else if (num_1 < 100){
    num_1_name <- paste0("00", num_1)
  } else if (num_1 < 1000){
    num_1_name <- paste0("0", num_1)
  } else {
    num_1_name <- num_1
  }
  
  # If num_2 is not a 4-digit number, add zeroes to the front to pad it
  if (num_2 < 10){
    num_2_name <- paste0("000", num_2)
  } else if (num_2 < 100){
    num_2_name <- paste0("00", num_2)
  } else if (num_2 < 1000){
    num_2_name <- paste0("0", num_2)
  } else {
    num_2_name <- num_2
  }
  
  # Export
  terra::writeRaster(changed_res, file.path(target_folder, paste0(target_name, num_1_name, "_", num_2_name, ".tif")), overwrite = T)
  
}

message("Export to cluster storage complete. If you want this file on the MaloneLab Server, please manually upload it there.")