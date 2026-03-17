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

# Point to the Landsat Salinity Model folder
landsat_salinity_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP_Landsat_Salinity_Model") 

# Read in Landsat files
B01 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B01.tif"))
B02 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B02.tif"))
B03 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B03.tif"))
B04 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B04.tif"))
B05 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B05.tif"))
B06 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B06.tif"))
B07 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B07.tif"))
B09 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B09.tif"))
B10 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B10.tif"))
B11 <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B11.tif"))
Fmask <- terra::rast(file.path("harmonized_appeears_landsat_data", "ENP_Landsat_Fmask.tif"))

# Point to precip, solar radiation, temp files
precip_files <- dir(file.path(landsat_salinity_folder, "ENP_Precipitation_Landsat_res"), pattern = ".tif", full.names = T)
solar_rad_files <- dir(file.path(landsat_salinity_folder, "ENP_SolarRadiation_Landsat_res"), pattern = ".tif", full.names = T)
temp_files <- dir(file.path(landsat_salinity_folder, "ENP_AverageTemperature_Landsat_res"), pattern = ".tif", full.names = T)

# Read in one precip, solar rad, and temp raster to check
test_precip_r <- terra::rast(precip_files[1])
test_solar_rad_r <- terra::rast(solar_rad_files[1])
test_temp_r <- terra::rast(temp_files[1])

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
terra::crs(B01) == terra::crs(test_precip_r)
terra::crs(B01) == terra::crs(test_solar_rad_r)
terra::crs(B01) == terra::crs(test_temp_r)
terra::crs(B01) == terra::crs(ele_r)
terra::crs(B01) == terra::crs(slope_r)
terra::crs(B01) == terra::crs(dist_r)
terra::crs(B01) == terra::crs(DBHydro_sf_v2)

## --------------------------------------------------------- ##
#                Create Extraction Function -----
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

## --------------------------------------------------------- ##
#                 Extraction: Landsat -----
## --------------------------------------------------------- ##

# List our current bands together
landsat_bands_list <- list(B01, B02, B03, B04, B05,
                           B06, B07, B09, B10, B11,
                           Fmask)
# List their names
landsat_names_list <- list("B01", "B02", "B03", "B04", "B05",
                           "B06", "B07", "B09", "B10", "B11",
                           "Fmask")

# Create an empty list to store our extracted points
landsat_points_list <- list()

# For every band...
for (i in seq_along(landsat_bands_list)){
  # Extract the Landsat data for the points
  band_points <- extract_from_raster(raster = landsat_bands_list[[i]],
                                     point_shapefile = terra::vect(DBHydro_sf_v2),
                                     new_column_name = landsat_names_list[[i]])
  
  # Save to list
  landsat_points_list[[i]] <- band_points
}

landsat_points <- landsat_points_list %>%
  # Join all extracted Landsat points by ID, station, date columns
  purrr::reduce(dplyr::full_join, by = c("ID", "station", "date")) %>% 
  # Drop the rows where all the bands have NA values
  # These empty rows come from Landsat rasters that only cover a small portion of ENP on that day
  dplyr::filter(!dplyr::if_all(starts_with("B"), is.na)) %>% 
  # Apply the Landsat scale factors to convert to original values
  dplyr::mutate(dplyr::across(.cols = B01:B09, .fns = ~.x * 0.0001)) %>%
  dplyr::mutate(dplyr::across(.cols = B10:B11, .fns = ~.x * 0.01)) %>%
  # Calculate indices
  dplyr::mutate(NDVI = (B05 - B04) / (B05 + B04),
                SI = (B03*B04)^0.5,
                NLI = (B05^2 - B04)/(B05^2 + B04),
                SRSI = ((NDVI - 1)^2 + SI^2)^0.5,
                S7 = (B06 - B07)/(B06 + B07),
                CRSI = ((B05*B04-B03*B02)/((B05*B04+B03*B02)))^0.5,
                NDSI = (B05 - B06) / (B05 + B06)) %>%
  # Rename date column to landsat_date
  dplyr::rename(landsat_date = date)

# Create a formatted date column
landsat_points$formatted_date <- as.Date(landsat_points$landsat_date)

## --------------------------------------------------------- ##
#       Extraction: Precip, Solar Radiation, Temp -----
## --------------------------------------------------------- ##

# Precipitation -----------------------------------------------

# Create an empty list to store extracted points
precip_points_list <- list()
# For every precip raster...
for (i in 1:length(precip_files)){
  # Read in the raster
  precip_r <- terra::rast(precip_files[i])
  # Extract the precip data for the points
  extracted_precip_points <- extract_from_raster(raster = precip_r, 
                                                 point_shapefile = terra::vect(DBHydro_sf_v2), 
                                                 new_column_name = "precip")
  # Save to list
  precip_points_list[[i]] <- extracted_precip_points
}

# Combine all extracted points
precip_points <- precip_points_list %>%
  purrr::map_dfr(.f = select, everything())

# Solar Radiation ---------------------------------------------

# Create an empty list to store extracted points
solar_rad_points_list <- list()
# For every solar radiation raster...
for (i in 1:length(solar_rad_files)){
  # Read in the raster
  solar_rad_r <- terra::rast(solar_rad_files[i])
  # Extract the solar radiation data for the points
  extracted_solar_rad_points <- extract_from_raster(raster = solar_rad_r, 
                                                    point_shapefile = terra::vect(DBHydro_sf_v2), 
                                                    new_column_name = "srad")
  # Save to list
  solar_rad_points_list[[i]] <- extracted_solar_rad_points
}

# Combine all extracted points
solar_rad_points <- solar_rad_points_list %>%
  purrr::map_dfr(.f = select, everything())

# Average Temperature -----------------------------------------

# Create an empty list to store extracted points
temp_points_list <- list()
# For every temperature raster...
for (i in 1:length(temp_files)){
  # Read in the raster
  temp_r <- terra::rast(temp_files[i])
  # Extract the temperature data for the points
  extracted_temp_points <- extract_from_raster(raster = temp_r, 
                                               point_shapefile = terra::vect(DBHydro_sf_v2), 
                                               new_column_name = "tavg")
  # Save to list
  temp_points_list[[i]] <- extracted_temp_points
}

# Combine all extracted points
temp_points <- temp_points_list %>%
  purrr::map_dfr(.f = select, everything())

# Combining all meteorology data ------------------------------

# List our extracted meteorology data
met_points_list <- list(precip_points, solar_rad_points, temp_points)

met_points <- met_points_list %>%
  # Join all extracted meteorology points by ID, station, date columns
  purrr::reduce(dplyr::full_join, by = c("ID", "station", "date")) %>% 
  # Order rows by date
  dplyr::arrange(date) 

# Create a formatted date column
met_points$formatted_date <- as.Date(met_points$date)

met_points <- met_points %>%
  # Drop old date column
  dplyr::select(-date)
  
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

DBHydro_sal_df <- DBHydro_df %>% 
  # Grab only relevant columns
  dplyr::select("station", "collectDate", "value", "grab") %>%
  # Rename value column to "salinity"
  dplyr::rename(salinity = value)

# Create a formatted date column
DBHydro_sal_df$formatted_date <- as.Date(DBHydro_sal_df$collectDate)

DBHydro_sal_df <- DBHydro_sal_df %>%
  # Drop old date column
  dplyr::select(-collectDate)

sal_met_ele_slope_dist <- DBHydro_sal_df %>%
  # Left join salinity with extracted met points
  dplyr::left_join(met_points, by = c("station", "formatted_date")) %>%
  # Left join salinity+met points with elevation+slope+distance
  dplyr::left_join(ele_slope_dist_points, by = c("ID", "station")) 

# Finally full join salinity+met+elevation+slope+distance with extracted Landsat points
DBSAL <- dplyr::full_join(sal_met_ele_slope_dist, landsat_points, by = c("ID", "station", "formatted_date")) %>%
  # Drop redundant ID column
  dplyr::select(-ID) %>%
  # Reorder columns
  dplyr::relocate(formatted_date, .after = station) %>%
  dplyr::relocate(landsat_date, .after = formatted_date) %>%
  dplyr::relocate(grab, .after = landsat_date) %>%
  # Create a flag column for rows with salinity measurements
  dplyr::mutate(has_salinity = dplyr::case_when(
    !is.na(salinity) ~ 1,
    T ~ 0
  ), .after = formatted_date) %>%
  # Create a flag column for rows with Landsat measurements
  dplyr::mutate(has_landsat = dplyr::case_when(
    !is.na(landsat_date) ~ 1,
    T ~ 0
  ), .after = has_salinity) %>%
  # Drop landsat_date column
  dplyr::select(-landsat_date) %>%
  # Order by salinity_date and station
  dplyr::arrange(formatted_date, station)

# Export all salinity + extracted data as CSV
readr::write_csv(DBSAL, "DBSAL.csv")

