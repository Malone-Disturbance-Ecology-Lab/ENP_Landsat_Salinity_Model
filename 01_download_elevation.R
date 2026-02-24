## --------------------------------------------- ##
#            Download Elevation Data
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script downloads elevation data for Everglades National Park using elevatr.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(tidyverse)
library(terra)
library(sf)
library(elevatr)

# Point to the folder with the ENP shapefile
shapefile_folder <- file.path("/", "corellia.environment.yale.edu", "MaloneLab", "Research", "ENP", "shapefiles")

# Read it in
enp <- sf::read_sf(file.path(shapefile_folder, "Everglades_NP_4326.shp"))

# Read in one raster to use as a template 
template_raster <- terra::rast("HLSL30.020_B01_doy2013111_aid0001_17N.tif")

## --------------------------------------------- ##
#                Downloading -----
## --------------------------------------------- ##

# Grab elevation data for ENP
ele <- elevatr::get_elev_raster(enp, z = 12) %>%
  # Convert to terra object
  terra::rast()

# Fix name
names(ele) <- "elevation"

# Resample to the template raster's resolution
ele_resample <- terra::resample(ele, template_raster)

# Mask elevation to just the ENP boundary
elevation <- terra::mask(ele_resample, enp)

# Export elevation raster
writeRaster(elevation, file = "ENP_Elevation.tif")