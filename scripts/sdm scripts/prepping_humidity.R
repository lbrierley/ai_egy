## Adapted from script by Hayes et al., date 21/9/23

rm(list = ls())

library(terra)
library(raster)
library(rnaturalearth)

# Set area of interest
egy_map <- ne_countries(country = "Egypt", returnclass = "sf")

# Read in 2022 daily data from Copernicus

files_list <- list.files("data/env_vars/humidity/humidity_raw", pattern="\\.nc$")
files_list_full <- paste("data/env_vars/humidity/humidity_raw/",files_list, sep = "")

# Crop to map of interest and resave
terra::rast(files_list_full) %>% crop(egy_map) %>% writeRaster("data/env_vars/humidity/humidity_daily.tif")
