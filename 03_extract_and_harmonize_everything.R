## --------------------------------------------- ##
#     Extract to Station Points and Harmonize
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script extracts and harmonizes these data sources into a clean data frame:
## Landsat, precipitation, solar radiation, temperature,
## elevation, slope, distance to coast, DBHydro salinity 

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(tidyverse)
library(terra)
library(sf)
library(tsibble)

# Point to Margo's Salinity Model folder
salinity_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP_Salinity_Model") 

# Read in Landsat files
B01 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B01.tif"))
B02 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B02.tif"))
B03 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B03.tif"))
B04 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B04.tif"))

# Read in precip, solar radiation, temp files


# Read in elevation, slope, distance to coast files
ele_r <- terra::rast(file.path("ENP_Elevation.tif"))
slope_r <- terra::rast(file.path("ENP_Slope.tif"))
dist_r <- terra::rast(file.path("ENP_DistCoast.tif"))

# Read in DBHydro files
DBHydro_sf <- sf::st_read(file.path("ENP_DBHydro_sf", "ENP_DBHydro_sf.shp"))
DBHydro_df <- readr::read_csv(file.path("ENP_DBHydro_salinity_df.csv"))
# Make the DBHydro station shapefile have the same CRS as Landsat just in case
DBHydro_sf_v2 <- DBHydro_sf %>%
  sf::st_transform(sf::st_crs(B01))

# Make sure all CRS are the same as the Landsat
terra::crs(B01) == terra::crs(DBHydro_sf_v2)
terra::crs(B01) == terra::crs(ele_r)
terra::crs(B01) == terra::crs(slope_r)
terra::crs(B01) == terra::crs(dist_r)

## --------------------------------------------------------- ##
#                 Extraction: Landsat -----
## --------------------------------------------------------- ##

extract_from_raster <- function(raster, point_shapefile, new_column_name){
  # Extract raster data for the DBHydro station points
  r_points <- terra::extract(raster, point_shapefile) %>%
    # Create a station column
    dplyr::mutate(station = point_shapefile$station)
  
  # Add 1 to the number of layers in the raster
  ending_num <- length(time(raster)) + 1
  
  # Rename the columns to the measurement dates
  names(r_points)[2:ending_num] <- time(raster) %>% as.character()
  
  # Pivot the data frame longer
  r_points_long <- r_points %>%
    tidyr::pivot_longer(
      cols = matches("[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}"),
      names_to = "date",
      values_to = new_column_name
    )
  
  return(r_points_long)
}

B01_points <- extract_from_raster(raster = B01, 
                                  point_shapefile = terra::vect(DBHydro_sf_v2), 
                                  new_column_name = "B01")

B02_points <- extract_from_raster(raster = B02, 
                                  point_shapefile = terra::vect(DBHydro_sf_v2), 
                                  new_column_name = "B02")

B03_points <- extract_from_raster(raster = B03, 
                                  point_shapefile = terra::vect(DBHydro_sf_v2), 
                                  new_column_name = "B03")

B04_points <- extract_from_raster(raster = B04, 
                                  point_shapefile = terra::vect(DBHydro_sf_v2), 
                                  new_column_name = "B04")

landsat_points_list <- list(B01_points, B02_points, B03_points, B04_points)

landsat_points <- landsat_points_list %>%
  # Join all extracted Landsat points by ID, station, date columns
  purrr::reduce(dplyr::full_join, by = c("ID", "station", "date")) %>% 
  # Create a YearWeek column
  dplyr::mutate(YearWeek = tsibble::yearweek(date)) %>% 
  # Drop the old date column
  dplyr::select(-date)

## --------------------------------------------------------- ##
#       Extraction: Precip, Solar Radiation, Temp -----
## --------------------------------------------------------- ##





## --------------------------------------------------------- ##
#    Extraction: Elevation, Slope, Distance to Coast -----
## --------------------------------------------------------- ##

# Stack elevation, slope, distance to coast rasters
ele_slope_dist_stack <- c(ele_r, slope_r, dist_r) 

# Extract raster data for the DBHydro station points
ele_slope_dist_points <- terra::extract(ele_slope_dist_stack, terra::vect(DBHydro_sf_v2)) %>%
  # Create a station column
  dplyr::mutate(station = DBHydro_sf_v2$station)

## --------------------------------------------------------- ##
#                       Harmonizing -----
## --------------------------------------------------------- ##

DBHydro_sal_df <- DBHydro_df %>% dplyr::select("station", "collectDate", "value", "grab") %>%
  dplyr::rename(salinity = value)

DBHydro_sal_df$day <- as.Date(DBHydro_sal_df$collectDate)



