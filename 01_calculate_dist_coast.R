## --------------------------------------------- ##
#            Find Distance to Coast
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script calculates distance to the coast.
## NOTE: run on Bouchet cluster for fast computation.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(tidyverse)
library(terra)
library(sf)

# Point to the Landsat Salinity folder on cluster
my_folder <- '/home/ac3656/GitHub/ENP_Landsat_Salinity_Model'

# Read in FL coastline shapefile
# Can also be found on file.path("/", "Volumes", "malonelab", "Research", "ENP_Salinity_Model", "florida_shoreline", "Florida_Shoreline_(1_to_12%2C000_Scale).shp") 
FL_coastline <- sf::st_read(file.path(my_folder, "florida_shoreline", "Florida_Shoreline_(1_to_12%2C000_Scale).shp"))

# Read in ENP shapefile
# Can also be found on file.path("/", "Volumes", "malonelab", "Research", "ENP", "shapefiles", "Everglades_NP_4326.shp")
enp <- sf::st_read(file.path(my_folder, "Everglades_NP_4326", "Everglades_NP_4326.shp"))

# Read in one raster to use as a template 
template_raster <- terra::rast(file.path(my_folder, "appeears_landsat_data", "L30", "B01", "HLSL30.020_B01_doy2013111_aid0001_17N.tif"))

# Make the CRS of the ENP boundary shapefile the same as the template raster just in case
enp_v2 <- enp %>%
  sf::st_transform(sf::st_crs(template_raster))

## --------------------------------------------- ##
#               Calculating -----
## --------------------------------------------- ##

# Change CRS to match ENP shapefile
ENP_coastline <- sf::st_transform(FL_coastline, crs = sf::st_crs(enp_v2)) %>% 
  # Grab the coastline that intersects with ENP
  sf::st_intersection(enp_v2)

# ggplot() +
#   geom_sf(data = ENP_coastline)

# Convert to a linestring geometry
ENP_coastline_line <- sf::st_cast(ENP_coastline, "LINESTRING") %>% 
  # Convert to terra object
  terra::vect() %>%
  # Make sure it has the same CRS as our template raster
  terra::project(terra::crs(template_raster))

# terra::plot(ENP_coastline_line)

# Make all fields with linestring = 1 
coastline <- terra::rasterize(ENP_coastline_line, template_raster, field = 1, touches = T) 

# terra::plot(coastline)

# Calculate distance to the coast
distance.Coast <- terra::distance(coastline, unit="m") # calculate distance of all other fields from 1 values

# terra::plot(distance.Coast)

# Fix name
names(distance.Coast) <- "distCoast"

# Mask to ENP boundary
distance.Coast.ENP <- terra::mask(distance.Coast, enp_v2)

# Export
terra::writeRaster(distance.Coast.ENP, file.path(my_folder, "ENP_DistCoast.tif"), overwrite = T)