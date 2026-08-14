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

# CHANGE AS NEEDED (also see towards bottom of script) ----

# Export to server as well as locally? 0 for no, 1 for yes
export_server <- 0

# ---------------------------------------------------------

if (export_server == 1){
  # Point to the Landsat Salinity Model folder
  landsat_salinity_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP_Landsat_Salinity_Model") 
  
  # Create new folder to store rasters
  dir.create(path = file.path(landsat_salinity_folder, "harmonized_appeears_landsat_data"), showWarnings = F)
  
  # Create new folder to store rasters
  dir.create(path = file.path(landsat_salinity_folder, "harmonized_appeears_landsat_data", "L30"), showWarnings = F)
  
  # Create new folder to store rasters
  dir.create(path = file.path(landsat_salinity_folder, "harmonized_appeears_landsat_data", "S30"), showWarnings = F)
  
}

# Create new folder to store rasters
dir.create(path = file.path("harmonized_appeears_landsat_data"), showWarnings = F)

# Create new folder to store rasters
dir.create(path = file.path("harmonized_appeears_landsat_data", "L30"), showWarnings = F)

# Create new folder to store rasters
dir.create(path = file.path("harmonized_appeears_landsat_data", "S30"), showWarnings = F)

## --------------------------------------------- ##
#                   Function:
#   Fixing rasters with different extents -----
## --------------------------------------------- ##

fix_extent <- function(band_name, type){
  
  # Get list of relevant tif files
  band_files_v0 <- dir(file.path("appeears_landsat_data", type, band_name), pattern = band_name, full.names = T)
  
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
  
  # Note the scale factor for different bands
  # Bands 10 and 11 have a scale factor of 0.01
  # Fmask has a scale factor of 1
  # All other bands have a scale factor of 0.0001
  
  if (type == "L30"){
    if (band_name == "B10" | band_name == "B11"){
      scale_factor <- 0.01 
    } else if (band_name == "Fmask"){
      scale_factor <- 1
    } else {
      scale_factor <- 0.0001
    }
  }
  
  # All S30 bands (except for Fmask) have a scale factor of 0.0001
  
  if (type == "S30"){
    if (band_name == "Fmask"){
      scale_factor <- 1
    } else {
      scale_factor <- 0.0001
    }
  }
  
  message(paste("Scale factor for this band:", scale_factor))
  
  # Create empty list to store fixed rasters
  fixed_rasters_list <- list()
  # For every raster with a different extent...
  for (i in seq_along(unlist(diff_extent))){
    # Read in raster
    wrong_extent <- terra::rast(unlist(diff_extent)[i])
    # Fix extent 
    # resample() will automatically divide by the scaling factor, so 
    # multiply by the scaling factor to undo this
    if (band_name == "Fmask"){
      fixed_raster <- terra::resample(wrong_extent, template, method = "near") * scale_factor
    } else {
      fixed_raster <- terra::resample(wrong_extent, template) * scale_factor
    }
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

add_dates <- function(band_name, harmonized_band, type){
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
                     file.path("harmonized_appeears_landsat_data", type, paste0(type, "_ENP_", band_name, ".tif")),
                     overwrite = T)
  
  if (export_server == 1){
    terra::writeRaster(harmonized_band, 
                       file.path(landsat_salinity_folder, "harmonized_appeears_landsat_data", type, paste0(type, "_ENP_", band_name, ".tif")),
                       overwrite = T)
  }
}

## --------------------------------------------- ##
#             Harmonizing Each Type -----
## --------------------------------------------- ##

# CHANGE AS NEEDED --------------------------------

# Needed layers L30: "B01", "B02", "B03", "B04", "B05", "B06", "B07", "Fmask"
# Needed layers S30: "B01", "B02", "B03", "B04", "B8A", "B11", "B12", "Fmask"
# Harmonize as needed

my_type <- "S30"
my_band <- "Fmask"

# -------------------------------------------------

band_fix_extent <- fix_extent(band_name = my_band, type = my_type)
add_dates(band_name = my_band, harmonized_band = band_fix_extent, type = my_type)
