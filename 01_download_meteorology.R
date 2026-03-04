## --------------------------------------------- ##
#              Download Meteorology
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script downloads meteorology data for Everglades National Park using climateR::getGridMET().
## NOTE: run on Grace cluster for fast computation.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(tidyverse)
library(terra)
library(sf)
library(climateR)

# Point to the Landsat Salinity folder
my_folder <- '/home/ac3656/ENP_Landsat_Salinity_Model'

# Read in ENP shapefile
# Can also be found on file.path("/", "Volumes", "malonelab", "Research", "ENP", "shapefiles", "Everglades_NP_4326.shp")
enp <- sf::read_sf(file.path(my_folder, "Everglades_NP_4326", "Everglades_NP_4326.shp"))

# Read in one raster to use as a template 
template_raster <- terra::rast(file.path(my_folder, "appeears_landsat_data", "B01", "HLSL30.020_B01_doy2013111_aid0001_17N.tif"))

# Make the CRS of the ENP boundary shapefile the same as the template raster just in case
enp_v2 <- enp %>%
  sf::st_transform(sf::st_crs(template_raster))

## --------------------------------------------- ##
#                Get Meteorology -----
## --------------------------------------------- ##

# Grab meteorology data
climate_rast <-  climateR::getGridMET(enp_v2, c("pr", "tmmn", "tmmx", "srad"),
                                      startDate = "2013-04-01",
                                      endDate = "2026-02-13")

# Pick out the precip, shortwave radiation, and temp rasters
precip_r <- climate_rast$precipitation_amount
srad_r <- climate_rast$daily_mean_shortwave_radiation_at_surface
tmin_r <- climate_rast$daily_minimum_temperature
tmax_r <- climate_rast$daily_maximum_temperature
tavg_r <- mean(tmin_r, tmax_r)

precip_r_v2 <- precip_r %>%
  terra::project(terra::crs(template_raster))

srad_r_v2 <- srad_r %>%
  terra::project(terra::crs(template_raster))

tavg_r_v2 <- tavg_r %>%
  terra::project(terra::crs(template_raster))

# Export rasters
terra::writeRaster(precip_r_v2, file.path(my_folder, "ENP_Precipitation.tif"), overwrite = T)
terra::writeRaster(srad_r_v2, file.path(my_folder,"ENP_SolarRadiation.tif"), overwrite = T)
terra::writeRaster(tavg_r_v2, file.path(my_folder,"ENP_AverageTemperature.tif"), overwrite = T)

