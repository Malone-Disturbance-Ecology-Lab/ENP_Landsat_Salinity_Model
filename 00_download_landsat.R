## --------------------------------------------- ##
#             Download Landsat Data
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script creates a download request for the Harmonized Landsat and Sentinel-2 Land Surface Reflectance rasters (HLSL30.020)
## using the NASA AppEEARS API. 
## The Landsat rasters will cover the Everglades National Park (ENP). 

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(sf)
library(appeears)

# Specify your NASA Earth Data username 
my_user <- "anchen14"

# Enter your NASA Earth Data credentials 

# NOTE:
# There are 2 ways to enter your Earth Data credentials

# Method 1
# Run the following in the console with your own Earth Data username and password:
# rs_set_key(
#   user = "earth_data_user",
#   password = "XXXXXXXXXXXXXXXXXXXXXX"
# )

# Method 2 (if Method 1 fails)
# If this is your first time setting up your credentials, run:
# options(keyring_backend = "file")
# Then on the console, run the following with your own Earth Data username and password:
# rs_set_key(
#   user = "earth_data_user",
#   password = "XXXXXXXXXXXXXXXXXXXXXX"
# )
# And then set an additional local keyring password you can remember when prompted

# If you follow all the steps correctly, you can unlock your credentials 
# at the start of every session with just this line:
options(keyring_backend = "file")

# For more help, see https://github.com/bluegreen-labs/appeears?tab=readme-ov-file#setup

## --------------------------------------------- ##
#        Use ENP boundary shapefile -----
## --------------------------------------------- ##

# Point to the folder with the ENP shapefile
shapefile_folder <- file.path("/", "corellia.environment.yale.edu", "MaloneLab", "Research", "ENP", "shapefiles")

# Read it in
enp <- sf::read_sf(file.path(shapefile_folder, "Everglades_NP_4326.shp"))

## --------------------------------------------- ##
#           Create AppEEARS task -----
## --------------------------------------------- ##

# Needed layers: "B01", "B02", "B03", "B04", "B05", "B06", "B07", "B09", "B10", "B11", "Fmask"
# Create AppEEARS requests as needed

my_task <- "surf_reflect_ENP_B11"
my_layers <- "B11"

# Create a dataframe for your AppEEARS task
df <- data.frame(
  task = my_task, # name of task
  subtask = "subtask", # name of subtask 
  start = "2013-04-01", # start date for data
  end = "2026-02-13", # end date for data
  product = "HLSL30.020", # data product ID
  layer = my_layers # name of specific band(s)
)

# Build the area-based task
task <- appeears::rs_build_task(
  df = df,
  roi = enp,
  format = "geotiff"
)

# Request the task to be executed
# If prompted, enter your local keyring password
appeears::rs_request(
  request = task,
  user = my_user,
  transfer = FALSE,
  verbose = TRUE
)

# Check your email for a link to download the requested AppEEARS files manually
