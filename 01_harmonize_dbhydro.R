## --------------------------------------------- ##
#           Harmonize Salinity Data
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script harmonizes DBHydro salinity data (continuous and grab measurements) 
## into a clean data frame.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(tidyverse)
library(sf)

# Create new folder to store shapefile
dir.create(path = file.path("ENP_DBHydro_sf"), showWarnings = F)

# Point to Margo's Salinity Model folder
salinity_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP_Salinity_Model") 

# Point to the folder with the ENP shapefile
shapefile_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP", "shapefiles") 

# Read ENP shapefile
ENP <- sf::read_sf(file.path(shapefile_folder, "Everglades_NP_4326.shp"))

## --------------------------------------------- ##
#         Read in cont and grab data -----
## --------------------------------------------- ##

# List grab files
grab_files <- dir(file.path("dbhydro_salinity", "grab"), pattern = ".csv", full.names = T)

# Combine grab files
grab_combined <- grab_files %>%
  purrr::map(.f = ~readr::read_csv(., skip = 27)) %>%
  purrr::map_dfr(.f = select, everything())

# List cont files
cont_files <- dir(file.path("dbhydro_salinity", "cont"), pattern = ".csv", full.names = T)
# Remove cont files that need different skip numbers
cont_files <- setdiff(cont_files, c(file.path("dbhydro_salinity", "cont", "sfwmd-data-202602171356.csv"), 
                                    file.path("dbhydro_salinity", "cont", "sfwmd-data-202602171419.csv"),
                                    file.path("dbhydro_salinity", "cont", "sfwmd-data-202602181244.csv")))

# Combine cont files that need to skip 44 lines                      
cont_1 <- cont_files %>%
  purrr::map(.f = ~readr::read_csv(., col_types = cols(TIMESERIESID = col_character(), SOURCE_TIMESERIES = col_character()),
                                   skip = 44)) %>%
  purrr::map_dfr(.f = select, everything())

# Read in rest of cont files with their respective skip lines
cont_2 <- readr::read_csv(file.path("dbhydro_salinity", "cont", "sfwmd-data-202602171356.csv"), 
                          col_types = cols(TIMESERIESID = col_character(), SOURCE_TIMESERIES = col_character()),
                          skip = 48)
cont_3 <- readr::read_csv(file.path("dbhydro_salinity", "cont", "sfwmd-data-202602171419.csv"), 
                          col_types = cols(TIMESERIESID = col_character(), SOURCE_TIMESERIES = col_character()),
                          skip = 53)
cont_4 <- readr::read_csv(file.path("dbhydro_salinity", "cont", "sfwmd-data-202602181244.csv"), 
                          col_types = cols(TIMESERIESID = col_character(), SOURCE_TIMESERIES = col_character()),
                          skip = 57)

# Combine all cont data frames
cont_combined <- list(cont_1, cont_2, cont_3, cont_4) %>%
  purrr::map_dfr(.f = select, everything())

# Read in latlon data frame
cont_latlon <- readr::read_csv(file.path(salinity_folder, "sfwmd-data-continuous", "latlon.csv")) %>% 
  dplyr::rename(station = STATION,
                latitude = LAT,
                longitude = LONG)

## --------------------------------------------- ##
#               Harmonizing -----
## --------------------------------------------- ##

grab_v1 <- grab_combined %>%
  # Create a grab column
  dplyr::mutate(grab = 1) %>%
  # Fix longitude column
  dplyr::mutate(longitude = dplyr::case_when(
    longitude > 0 ~ longitude * -1,
    T ~ longitude
  )) %>% 
  # Select required columns
  dplyr::select("station", "latitude", "longitude", "collectDate", "value", "grab")

cont_v1 <- cont_combined %>%
  # Select and rename station, collectDate, value columns
  dplyr::select(STATION, TIMESTAMP, VALUE) %>%
  dplyr::rename(station = STATION,
                collectDate = TIMESTAMP,
                value = VALUE) %>%
  # Create a grab column and
  # Remove the leading 0 from station column
  dplyr::mutate(grab = 0,
                station = sub("^0+", "", station))

cont_v2 <- cont_v1 %>%
  # Join cont data frame with latlon
  dplyr::left_join(cont_latlon, by = "station") %>% 
  # Select required columns
  dplyr::select("station", "latitude", "longitude", "collectDate", "value", "grab")

# Finally combine cont and grab data
DBHydro_df <- dplyr::bind_rows(cont_v2, grab_v1)

## --------------------------------------------- ##
#               Exporting -----
## --------------------------------------------- ##

# Export combined cont and grab data as CSV
readr::write_csv(DBHydro_df, "ENP_DBHydro_salinity_df.csv")

DBHydro_points_df <- DBHydro_df %>%
  # Select required columns
  dplyr::select(station, longitude, latitude, grab) %>% 
  # Get distinct station points
  dplyr::distinct()

# Export distinct station points as CSV
readr::write_csv(DBHydro_points_df, "DBHydro_lonlat.csv")

# Convert distinct station points to shapefile
DBHydro_sf <- sf::st_as_sf(DBHydro_points_df, coords = c("longitude", "latitude"), crs = 4326)
DBHydro_sf <- DBHydro_sf[sf::st_within(DBHydro_sf, ENP, sparse = FALSE), ]

# Export distinct station points as shapefile
sf::st_write(DBHydro_sf, file.path("ENP_DBHydro_sf", "ENP_DBHydro_sf.shp"), append = FALSE)
