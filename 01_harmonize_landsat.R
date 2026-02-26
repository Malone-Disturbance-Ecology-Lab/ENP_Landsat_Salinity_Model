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

## --------------------------------------------- ##
#                   Function:
#   Fixing rasters with different extents -----
## --------------------------------------------- ##

fix_extent <- function(band_name){
  
  # Get list of relevant tif files
  band_files_v0 <- dir(file.path("appeears_landsat_data", band_name), pattern = band_name, full.names = T)
  
  # Create empty list for rasters that have a different extent
  diff_extent <- list()
  # For every tif file...
  for (i in seq_along(band_files_v0)){
    # Grab its metadata
    size <- terra::describe(band_files_v0[i])
    # If the extent is different...
    if (stringr::str_detect(size[3], "4171, 3839") == FALSE){
      # Add to list
      diff_extent <- append(diff_extent, band_files_v0[i])
    }
  }
  
  message(paste("Number of rasters with differing extents:", length(diff_extent)))
  
  # Get list of relevant tif files without the ones that have a different extent
  band_files_v1 <- setdiff(band_files_v0, unlist(diff_extent))
  # Read in one good raster to use as a template
  template <- terra::rast(band_files_v1[1])
  
  message("Fixing raster extents")
  
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
  
  # Combine the rasters that didn't need fixing into one SpatRaster
  band <- terra::rast(band_files_v1)
  
  message("Combining all rasters")
  
  # Combine all rasters together
  band_combined <- c(band, fixed_rasters)
  # Sort raster layers by alphabetical names
  band_combined_sorted <- terra::subset(band_combined, order(names(band_combined)))
  # Fix varnames
  varnames(band_combined_sorted) <- names(band_combined_sorted)
  
  message("Adding min/max statistics to metadata")
  
  # Add min/max statistics to metadata for completeness
  setMinMax(band_combined_sorted)
}

## --------------------------------------------- ##
#                   Function:
#      Adding dates to raster metadata -----
## --------------------------------------------- ##

add_dates <- function(band_name, harmonized_band){
  # Extract the year and doy info
  year_doy <- stringr::str_extract(names(harmonized_band), "[:digit:]{7}")
  
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
  time(harmonized_band) <- do.call("c", formatted_dates)
  
  message("Exporting raster, please be patient")
  
  # Export harmonized raster
  terra::writeRaster(harmonized_band, 
                     file.path("harmonized_appeears_landsat_data", paste0("ENP_Landsat_", band_name, ".tif")),
                     overwrite = T)
}

## --------------------------------------------- ##
#               Harmonizing -----
## --------------------------------------------- ##

# Needed layers: "B01", "B02", "B03", "B04", "B05", "B06", "B07", "B09", "B10", "B11"
# Harmonize as needed

my_band <- "B06"

band_fix_extent <- fix_extent(band_name = my_band)
add_dates(band_name = my_band, harmonized_band = band_fix_extent)
