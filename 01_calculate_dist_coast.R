## --------------------------------------------- ##
#            Find Distance to Coast
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script calculates distance to the coast.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(tidyverse)
library(terra)
library(sf)

# Point to Margo's Salinity Model folder
salinity_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP_Salinity_Model") 

# Point to the folder with the ENP shapefile
shapefile_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP", "shapefiles")

# Read in FL coastline shapefile
FL_coastline <- sf::st_read(file.path(salinity_folder, "florida_shoreline", "Florida_Shoreline_(1_to_12%2C000_Scale).shp"))

# Read it ENP shapefile
enp <- sf::st_read(file.path(shapefile_folder, "Everglades_NP_4326.shp"))

# Read in one raster to use as a template 
template_raster <- terra::rast(file.path("appeears_landsat_data", "B01", "HLSL30.020_B01_doy2013111_aid0001_17N.tif"))

## --------------------------------------------- ##
#               Calculating -----
## --------------------------------------------- ##

# Change CRS to match ENP shapefile
ENP_coastline <- sf::st_transform(FL_coastline, crs = sf::st_crs(enp)) %>% 
  # Grab the coastline that intersects with ENP
  sf::st_intersection(enp)

# ggplot() +
#   geom_sf(data = ENP_coastline)

ENP_coastline_line <- sf::st_cast(ENP_coastline, "LINESTRING") %>% 
  #as("Spatial") %>% 
  terra::vect()

# terra::plot(ENP_coastline_line)

coastline <- terra::rasterize(ENP_coastline_line, template_raster, field = 1, touches = T) # make all fields with line string = 1
# coastline %>% terra::plot()

# Calculate Distance to the coast:
distance.Coast <- terra::distance(coastline, unit="m", method="geo") # calculate distance of all other fields from 1 values
plot(distance.Coast)
names(distance.Coast) <- "distCoast"