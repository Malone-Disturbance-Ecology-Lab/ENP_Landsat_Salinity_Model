## --------------------------------------------- ##
#              Download Meteorology
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script downloads meteorology data for Everglades National Park
## using climateR::getGridMET().

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(tidyverse)
library(terra)
library(sf)
library(climateR)

# Point to the folder with the ENP shapefile
shapefile_folder <- file.path("/", "corellia.environment.yale.edu", "MaloneLab", "Research", "ENP", "shapefiles")

# Read it in
enp <- sf::read_sf(file.path(shapefile_folder, "Everglades_NP_4326.shp"))

# Read in one raster to use as a template 
template_raster <- terra::rast("HLSL30.020_B01_doy2013111_aid0001_17N.tif")

## --------------------------------------------- ##
#                Get Meteorology -----
## --------------------------------------------- ##

# Grab meteorology data
climate_rast <-  climateR::getGridMET(enp, c("pr", "tmmn", "tmmx", "srad"),
                            startDate = "2013-04-01",
                            endDate = "2026-02-13")

# Pick out the precip, shortwave radiation, and temp rasters
precip_r <- climate_rast$precipitation_amount
srad_r <- climate_rast$daily_mean_shortwave_radiation_at_surface
tmin_r <- climate_rast$daily_minimum_temperature
tmax_r <- climate_rast$daily_maximum_temperature
tavg_r <- mean(tmin_r, tmax_r)

# Resample to the template raster's resolution
precip <- terra::resample(precip_r, template_raster)
srad <- terra::resample(srad_r, template_raster)
tavg <- terra::resample(tavg_r, template_raster)

# Export precip, shortwave radiation, and temp rasters
terra::writeRaster(precip, "ENP_Precipitation.tif", overwrite = T)
terra::writeRaster(srad, "ENP_SolarRadiation.tif", overwrite = T)
terra::writeRaster(tavg, "ENP_AverageTemperature.tif", overwrite = T)
