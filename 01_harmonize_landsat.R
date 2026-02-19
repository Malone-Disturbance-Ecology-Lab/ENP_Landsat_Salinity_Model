## --------------------------------------------- ##
#           Harmonize Landsat Data
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script harmonizes Landsat rasters into a clean tif file for each band.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(tidyverse)
library(terra)

# Create new folder to store rasters
dir.create(path = file.path("harmonized_appeears_landsat_data"), showWarnings = F)

# Get list of relevant tif files
B01_files_v0 <- dir(file.path("appeears_landsat_data", "B01"), pattern = "B01", full.names = T)

## --------------------------------------------- ##
#   Fixing rasters with different extents -----
## --------------------------------------------- ##

# Create empty list for rasters that have a different extent
diff_extent <- list()
# For every tif file...
for (i in seq_along(B01_files_v0)){
  # Grab its metadata
  size <- terra::describe(B01_files_v0[i])
  # If the extent is different...
  if (stringr::str_detect(size[3], "4171, 3839") == FALSE){
    # Add to list
    diff_extent <- append(diff_extent, B01_files_v0[i])
  }
}

# Get list of relevant tif files without the ones that have a different extent
B01_files_v1 <- setdiff(B01_files_v0, unlist(diff_extent))
# Read in one good raster to use as a template
template <- terra::rast(B01_files_v1[1])

# Create empty list to store fixed rasters
fixed_rasters_list <- list()
# For every raster with a different extent...
for (i in seq_along(unlist(diff_extent))){
  # Read in raster
  wrong_extent <- terra::rast(unlist(diff_extent)[i])
  # Fix extent 
  # resample() will automatically divide by the scaling factor, so 
  # multiply by the scaling factor to undo this
  fixed_raster <- terra::resample(wrong_extent, template) * (0.0001)
  # Add to list
  fixed_rasters_list[[i]] <- fixed_raster
}

# Combine fixed rasters into one SpatRaster
fixed_rasters <- terra::rast(fixed_rasters_list) 
# Fix varnames
varnames(fixed_rasters) <- names(fixed_rasters)

# Combine rest of the rasters into one SpatRaster
B01 <- terra::rast(B01_files_v1)

# Combine every raster together
B01_combined <- c(B01, fixed_rasters)
# Sort raster layers by alphabetical names
B01_combined_sorted <- terra::subset(B01_combined, order(names(B01_combined)))
# Fix varnames
varnames(B01_combined_sorted) <- names(B01_combined_sorted)
# Add min/max statistics to metadata for completeness
setMinMax(B01_combined_sorted)

## --------------------------------------------- ##
#      Adding dates to raster metadata -----
## --------------------------------------------- ##

# Extract the year and doy info
year_doy <- stringr::str_extract(names(B01_combined_sorted), "[:digit:]{7}")

# Create empty list to store formatted dates
formatted_dates <- list()
# For every extracted year and doy string...
for (i in seq_along(year_doy)){
  # Grab the year
  year <- stringr::str_sub(year_doy[i], 1, 4)
  # Grab the doy
  doy <- as.numeric(stringr::str_sub(year_doy[i], 5, 7))
  
  # Set the origin date as the start of the year
  origin_date <- as.Date(paste0(year, "-01-01"))
  # Convert doy to formatted date
  formatted_date <- as.Date(doy - 1, origin = origin_date)
  
  # Add to list
  formatted_dates[[i]] <- formatted_date
}

# Add formatted dates to metadata for completeness
time(B01_combined_sorted) <- do.call("c", formatted_dates)

# Export harmonized raster
terra::writeRaster(B01_combined_sorted, file.path("harmonized_appeears_landsat_data", "ENP_Landsat_B01.tif"))
