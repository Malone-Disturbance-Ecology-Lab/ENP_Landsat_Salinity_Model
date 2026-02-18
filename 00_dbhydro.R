## --------------------------------------------- ##
#             Download Salinity Data
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script finds and downloads salinity data (continuous and grab measurements) 
## from DBHydro from 04/01/2013 to 2/13/2026.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary packages
library(tidyverse)

# Point to Margo's Salinity Model folder
folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP_Salinity_Model") 

## --------------------------------------------- ##
#          Download continuous data -----
## --------------------------------------------- ##

# Read in continuous salinity data
cont_ls <- list(
  DBHydro_cont_1 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-2025120200122.csv"), skip = 34)[ , c("STATION", "TIMESTAMP", "VALUE")],
  DBHydro_cont_2 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-2025120200121.csv"), skip = 36)[ , c("STATION", "TIMESTAMP", "VALUE")],
  DBHydro_cont_3 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-2025120200103.csv"), skip = 54)[ , c("STATION", "TIMESTAMP", "VALUE")],
  DBHydro_cont_4 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-2025120200102.csv"), skip = 51)[ , c("STATION", "TIMESTAMP", "VALUE")],
  DBHydro_cont_5 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-2025120200101.csv"), skip = 52)[ , c("STATION", "TIMESTAMP", "VALUE")],
  DBHydro_cont_6 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-2025120200062.csv"), skip = 51)[ , c("STATION", "TIMESTAMP", "VALUE")],
  DBHydro_cont_7 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-2025120200061.csv"), skip = 53)[ , c("STATION", "TIMESTAMP", "VALUE")],
  DBHydro_cont_8 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-202512020012.csv"), skip = 46)[ , c("STATION", "TIMESTAMP", "VALUE")],
  DBHydro_cont_9 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-202512020011.csv"), skip = 34)[ , c("STATION", "TIMESTAMP", "VALUE")],
  DBHydro_cont_10 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-202512020010.csv"), skip = 51)[ , c("STATION", "TIMESTAMP", "VALUE")],
  DBHydro_cont_11 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-202512020007.csv"), skip = 59)[ , c("STATION", "TIMESTAMP", "VALUE")],
  DBHydro_cont_12 <- read.csv(file.path(folder,"sfwmd-data-continuous/sfwmd-data-202512020006.csv"), skip = 62)[ , c("STATION", "TIMESTAMP", "VALUE")])

# Combine into one data frame
DBHydro_continuous <- cont_ls %>% 
  lapply(\(x) x %>% mutate(STATION = as.character(STATION))) %>% bind_rows()

# Find required stations
sort(unique(DBHydro_continuous$STATION))

# "02290930"   "22908295"   "EASTSIDE"   "EDEN3"      "ENPBA"      "ENPBK"     
# "ENPBN"     "ENPBR"      "ENPBS"      "ENPBSC"     "ENPCA"      "ENPCN"      
# "ENPCW"      "ENPDK"     "ENPGB"      "ENPGI"      "ENPHC"      "ENPHR"      
# "ENPJK"      "ENPLB"      "ENPLM"     "ENPLN"      "ENPLO"      "ENPLR"      
# "ENPLS"      "ENPMK"      "ENPNR"      "ENPPK"     "ENPSR"      "ENPTB"     
# "ENPTC"      "ENPTE"      "ENPTR"      "ENPWB"      "ENPWE"     "ENPWP"      
# "ENPWW"      "HIGHWAY_CR" "JBTS"       "JOEBAY2E"   "MCCORMICK"  "MUD_CRKM"  
# "SEVENPALM"  "TAYLORS3"   "TAYLORUPS"  "TAYSLOUWET" "TROUT CR_B"

# Click on download links below
# https://insightsdata.sfwmd.gov/#/shared-reports?timeseriesIds=89867,89865,89866,AL985,AL970,90111,90112,38108,63572,63578&reportType=timeseries&format=csv&startDate=20130401&endDate=20260213&module=continuous
# https://insightsdata.sfwmd.gov/#/shared-reports?timeseriesIds=63581,63584,63588,39533,63592,63596&reportType=timeseries&format=csv&startDate=20130401&endDate=20260213&module=continuous
# https://insightsdata.sfwmd.gov/#/shared-reports?timeseriesIds=63604,63608,63617,63621,63625,63630&reportType=timeseries&format=csv&startDate=20130401&endDate=20260213&module=continuous
# https://insightsdata.sfwmd.gov/#/shared-reports?timeseriesIds=63634,63638,63642,63646,63650,63654&reportType=timeseries&format=json&startDate=20130401&endDate=20260213&reportFormat=json&module=continuous
# https://insightsdata.sfwmd.gov/#/shared-reports?timeseriesIds=63663,63667,63671,63675,63679,63683&reportType=timeseries&format=csv&startDate=20130401&endDate=20260213&module=continuous
# https://insightsdata.sfwmd.gov/#/shared-reports?timeseriesIds=63687,63691,AL633,63695,63699,63703&reportType=timeseries&format=csv&startDate=20130401&endDate=20260213&module=continuous
# https://insightsdata.sfwmd.gov/#/shared-reports?timeseriesIds=63708,39208,39209,39212,39215,AN690,38066,38067,66106,66107,66112,66113,66125,66126,90464&reportType=timeseries&format=csv&startDate=20130401&endDate=20260213&module=continuous
# https://insightsdata.sfwmd.gov/#/shared-reports?timeseriesIds=38010,38012,38011,66108,66109,66115,66117,66128,66129,38019,38021,38020,38023,38094,66110,66111,66123,66124,66130,66131&reportType=timeseries&format=csv&startDate=20130401&endDate=20260213&module=continuous

## --------------------------------------------- ##
#            Download grab data -----
## --------------------------------------------- ##

# Required stations from 
# https://github.com/Malone-Disturbance-Ecology-Lab/ENP_Salinity_Model/blob/main/00_DBHydro.R

# "FLAB05"  "FLAB06"  "FLAB07"  "FLAB08"  "FLAB09"  "FLAB10"  
# "FLAB11"  "FLAB12"  "FLAB13"  "FLAB14"  "FLAB15"  "FLAB16"  
# "FLAB17"  "FLAB18"  "FLAB19"  "FLAB20"  "FLAB21"  "FLAB23" 
# "FLAB24"  "FLAB25"  "FLAB27"  "FLAB29"  "FLAB30"  "FLAB31"  
# "FLAB33"  "FLAB34"  "FLAB35"  "FLAB36"  "FLAB37"  "FLAB38"  
# "FLAB39"  "FLAB40"  "FLAB41"  "FLAB43"  "FLAB44"  "FLAB47"  
# "FLAB48"  "TTI51"   "TTI51B"  "TTI52"   "TTI55"   "TTI57"   "TTI59"

# Sadly no direct download links
# Go to https://insightsdata.sfwmd.gov/#/waterquality to download the grab data manually