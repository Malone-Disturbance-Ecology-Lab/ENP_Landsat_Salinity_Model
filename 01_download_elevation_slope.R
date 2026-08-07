## --------------------------------------------- ##
#       Download Elevation and Slope Data
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script downloads elevation and slope data for Everglades National Park
## using elevatr::get_elev_raster() and terra::terrain().

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(tidyverse)
library(terra)
library(sf)
library(elevatr)

# CHANGE AS NEEDED --------------------------------

# Export to server? 0 for no, 1 for yes
export_server <- 0

# -------------------------------------------------

if (export_server == 1){
  # Point to the Landsat Salinity Model folder
  landsat_salinity_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP_Landsat_Salinity_Model") 
}

# Point to the folder with the ENP shapefile
shapefile_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP", "shapefiles")

# Read it in
enp <- sf::read_sf(file.path(shapefile_folder, "Everglades_NP_4326.shp"))

# Read in one raster to use as a template 
template_raster <- terra::rast(file.path("appeears_landsat_data", "L30", "B01", "HLSL30.020_B01_doy2013111_aid0001_17N.tif"))

# Make the CRS of the ENP boundary shapefile the same as the template raster just in case
enp_v2 <- enp %>%
  sf::st_transform(sf::st_crs(template_raster))

## --------------------------------------------- ##
#                Get Elevation -----
## --------------------------------------------- ##

# Grab elevation data for ENP
# Set z = 12 in order to get 30 meter resolution later
ele <- elevatr::get_elev_raster(enp_v2, z = 12) %>%
  # Convert to terra object
  terra::rast()

# Fix name
names(ele) <- "elevation"

# Project and resample to the template raster's CRS and resolution, respectively
ele_project <- terra::project(ele, template_raster)

# Mask elevation to just the ENP boundary
elevation <- terra::mask(ele_project, enp_v2)

# Export elevation raster
terra::writeRaster(elevation, file = "ENP_Elevation.tif", overwrite = T)

if (export_server == 1){
  terra::writeRaster(elevation, file = file.path(landsat_salinity_folder, "ENP_Elevation.tif"), overwrite = T)
}

## --------------------------------------------- ##
#                  Get Slope -----
## --------------------------------------------- ##

# Grab slope data for ENP
slope <- terra::terrain(elevation, "slope")

# Fix name
names(slope) <- "slope"

# Export slope raster
terra::writeRaster(slope, file = "ENP_Slope.tif", overwrite = T)

if (export_server == 1){
  terra::writeRaster(slope, file = file.path(landsat_salinity_folder, "ENP_Slope.tif"), overwrite = T)
}