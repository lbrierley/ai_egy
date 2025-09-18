## Adapted from script by Hayes et al., date 21/3/24
rm(list = ls())

library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(rnaturalearth)

# Load data handling functions
source("scripts/functions.R")

# Set grid resolution for prediction
grid_res <- 0.1

# Set seasonal dates to run predictions at (must be given as month-day)
pred_dates <- c("01-01", "04-01", "07-01", "10-01")

## Create standard prediction grid
egy_map <- ne_countries(country = "Egypt", returnclass = "sf")
sf::st_write(egy_map, "egy_map.shp", delete_layer = TRUE)

egy_lims <- st_bbox(egy_map)

egy_ext <- terra::ext(24,37,21,32) # based on bounding box of Egypt above

egy_rast <- terra::rast(extent=egy_ext, 
                        res = grid_res,
                        crs = "EPSG:4326")

# Non-temporal covariates
# Resample to prediction grid by bilinear interpolation

pred_elev <- geodata::elevation_global(res = 0.5, path = "data/env_vars/elevation/geodata_global_res0.5") %>% 
  resample(egy_rast, method = "bilinear")

pred_cost <- terra::rast("data/env_vars/dist_to_coast/dist_to_coast.tif") %>% 
  project(egy_rast) %>%
  resample(egy_rast, method = "bilinear")

pred_wetl <- terra::rast("data/env_vars/lakes_and_wetlands/dist_wetland.tif") %>% 
  project(egy_rast) %>%
  resample(egy_rast, method = "bilinear")

pred_chik <- terra::rast("data/env_vars/livestock/chickens/6_Ch_2015_Aw.tif") %>% 
  resample(egy_rast, method = "bilinear")

pred_duck <- terra::rast("data/env_vars/livestock/ducks/6_Dk_2015_Aw.tif") %>% 
  resample(egy_rast, method = "bilinear")

# Resample land use by mode assignment (since categorical)

pred_land <- terra::rast("data/env_vars/landcover/landcover_type1_full_raster.tif") %>% 
  resample(egy_rast, method = "mode")

# Convert land use to binary rasters
# Identify used categories
categories <- pred_land %>% values %>% unique() %>% as.vector() %>% sort()
pred_land <- lapply(categories, function(i) as.numeric(pred_land == i)) %>%
  rast()

# Temporal covariates, matching raster layer of specified calendar date
# Resample to prediction grid by bilinear interpolation

# NDVI (16-day periods)
env_ndvi <- terra::rast("data/env_vars/ndvi/ndvi_lat_long.tif") 
ndvi_layer_index <- sapply(pred_dates, match_raster_layer, raster_ymd = env_ndvi)
ndvi_layers <- list()

for (i in 1:length(pred_dates)){
  ndvi_layers[[i]] <- env_ndvi[[ndvi_layer_index[i]]] %>%
    resample(egy_rast, method = "bilinear")
}

# Select date of interest to average for monthly covars and then do bilinear interpolation

# Mean relative humidity in the month prior
env_humd <- terra::rast("data/env_vars/humidity/humidity_daily.tif")
humd_layer_index <- sapply(pred_dates, match_raster_layer, raster_ymd = env_humd)
humd_layers <- list()

for (i in 1:length(pred_dates)){
  humd_layers[[i]] <- env_humd[[calendar_wraparound_month(humd_layer_index[i])]] %>%
    mean %>%
    resample(egy_rast, method = "bilinear")
}

# Minimum zero-degree isotherm in the month prior
env_ziso <- terra::rast("data/env_vars/isotherm/isocline_raw.grib") %>%
  project(egy_rast)
ziso_layer_index <- sapply(pred_dates, match_raster_layer, raster_ymd = env_ziso)
ziso_layers <- list()

for (i in 1:length(pred_dates)){
  ziso_layers[[i]] <- env_ziso[[calendar_wraparound_month(ziso_layer_index[i])]] %>%
    min %>%
    resample(egy_rast, method = "bilinear")
}

# Mean midday temperature in the month prior
env_temp <- terra::rast("data/env_vars/climate/temp_midday_daily.tif")
temp_layer_index <- sapply(pred_dates, match_raster_layer, raster_ymd = env_temp)
temp_layers <- list()

for (i in 1:length(pred_dates)){
  temp_layers[[i]] <- env_temp[[calendar_wraparound_month(temp_layer_index[i])]] %>%
    mean %>%
    resample(egy_rast, method = "bilinear")
}

# Mean diurnal temperature range in the month prior
env_diur <- terra::rast("data/env_vars/climate/temp_diurnal_daily.tif") %>% subst(from = -Inf, to = NA)
diur_layer_index <- sapply(pred_dates, match_raster_layer, raster_ymd = env_diur)
diur_layers <- list()

for (i in 1:length(pred_dates)){
  diur_layers[[i]] <- env_diur[[calendar_wraparound_month(diur_layer_index[i])]] %>%
    mean %>%
    resample(egy_rast, method = "bilinear")
}

# Total precipitation in the month prior
env_prec <- terra::rast("data/env_vars/climate/total_precip_daily.tif")
prec_layer_index <- sapply(pred_dates, match_raster_layer, raster_ymd = env_prec)
prec_layers <- list()

for (i in 1:length(pred_dates)){
  prec_layers[[i]] <- env_prec[[calendar_wraparound_month(prec_layer_index[i])]] %>%
    sum %>%
    resample(egy_rast, method = "bilinear")
}


## Save final dataset as seasonal rasters and csvs

for (i in 1:length(pred_dates)){
  
  temp_raster <- c(pred_elev,
                   pred_cost,
                   pred_wetl,
                   pred_chik,
                   pred_duck,
                   pred_land,
                   ndvi_layers[[i]],
                   humd_layers[[i]],
                   ziso_layers[[i]],
                   temp_layers[[i]],
                   diur_layers[[i]],
                   prec_layers[[i]]
                   
  )
  
  names(temp_raster) <- c("wc2.1_30s_elev",
                          "dist_to_coast", 
                          "dist_to_wetland",
                          "X6_Ch_2015_Aw",
                          "X6_Dk_2015_Aw",
                          paste0("dem", categories),
                          "ndvi",
                          "humd_m",
                          "ziso_m",
                          "temp_m",
                          "diur_m",
                          "prec_m"
  )
  
  temp_raster %>% terra::writeRaster(paste0("data/pred/pred_grid_",i,".tif"), overwrite = TRUE)
  
  temp_df <- temp_raster %>% as.data.frame(xy = TRUE)

  temp_df %>% write.csv(paste0("data/pred/pred_grid_",i,".csv"))
  
}

# Test plot
ggplot() +
  geom_raster(data = temp_raster, aes(x = x, y = y, fill = diur_m)) +
  geom_sf(data = egy_map[1], color = "black", alpha = 0)  +
  theme_minimal() +
  scale_fill_viridis_c(na.value = "transparent") +
  #  coord_sf(xlim = c(30,34), ylim = c(29,32)) +
  labs(x = "Longitude",
       y = "Latitude")


# Nicer plot

names(temp_raster) <- c("elevation",
                        "distance_coast", 
                        "distance_wetland",
                        "chicken_dens",
                        "duck_dens",
                        "land_category",
                        "ndvi",
                        "month_humidity",
                        "min_isotherm",
                        "month_temperature",
                        "month_diurn_range",
                        "month_total_rain"
)

temp_raster %>% mask(egy_map) %>% plot
