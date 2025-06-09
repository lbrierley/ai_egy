## Adapted from script by Hayes et al., date 21/3/24

# Load data handling functions
source("scripts/functions.R")

## Create map
egy_map <- ne_countries(country = "Egypt", returnclass = "sf")
sf::st_write(egy_map, "egy_map.shp", delete_layer = TRUE)

egy_lims <- st_bbox(egy_map)
egy_rast <- terra::rast(extent=terra::ext(24,37,21,32), res = 0.25)

point_data <- st_as_sf(egy_data, 
                       coords = c("Longitude", "Latitude"),
                       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")



ggplot() +
  geom_sf(data = egy_map[1], fill = "grey", color = "black", alpha = 0.4)  +
  geom_sf(data = point_data %>% arrange(Result, Source), aes(color = Result, shape = Source), size = 4, alpha = 0.7) +
  theme_minimal() +
  coord_sf(xlim = c(30,34), ylim = c(29,32)) +
  scale_shape_manual(values=c(19, 18, 15, 17)) +
  scale_color_manual(values=c('#56B4E9','#F8766D')) +
  labs(x = "Longitude",
       y = "Latitude",
       shape = "Source",
       color = "Result")



# How many unique spatial data points?

p <- egy_data %>%
  filter(Result == "positive") %>%
  select(Latitude, Longitude, Result) %>% 
  distinct()
n <- egy_data %>%
  filter(Result == "negative") %>%
  select(Latitude, Longitude, Result) %>% 
  distinct() %>%
  anti_join(p, by = c("Latitude", "Longitude")) # Remove any positives from same time/place

nrow(p)
nrow(n)

# # How many unique spatiotemporal data points?
# 
# p <- egy_data %>%
#   filter(Result == "positive") %>%
#   select(Latitude, Longitude, Date, Result) %>%
#   distinct()
# n <- egy_data %>%
#   filter(Result == "negative") %>%
#   select(Latitude, Longitude, Date, Result) %>%
#   distinct() %>%
#   anti_join(p, by = c("Latitude", "Longitude", "Date")) # Remove any positives from same time/place
# 
# nrow(p)
# nrow(n)

# Calendar day-wise

p <- egy_data %>%
  filter(Result == "positive") %>%
  mutate(Date = format(as.Date(Date), "%m-%d")) %>%
  select(Latitude, Longitude, Date, Result) %>% 
  distinct()
n <- egy_data %>%
  filter(Result == "negative") %>%
  mutate(Date = format(as.Date(Date), "%m-%d")) %>%
  select(Latitude, Longitude, Date, Result) %>% 
  distinct() %>%
  anti_join(p, by = c("Latitude", "Longitude", "Date")) # Remove any positives from same time/place

nrow(p)
nrow(n)

# # Month-wise
# 
# p <- egy_data %>%
#   filter(Result == "positive") %>%
#   mutate(Date = format(as.Date(Date), "%m")) %>%
#   select(Latitude, Longitude, Date, Result) %>%
#   distinct()
# n <- egy_data %>%
#   filter(Result == "negative") %>%
#   mutate(Date = format(as.Date(Date), "%m")) %>%
#   select(Latitude, Longitude, Date, Result) %>%
#   distinct() %>%
#   anti_join(p, by = c("Latitude", "Longitude", "Date")) # Remove any positives from same time/place
# 
# nrow(p)
# nrow(n)


# # Just make a degree based map rather than metre based map to start
# template <- rast(vect(egy_map),res=0.001)
# poly_raster <- rasterize(vect(egy_map), template)
# plot(poly_raster)
# 
# p <- point_data %>%
#   filter(Result == "positive") %>%
#   rasterize(y = poly_raster, fun = "max") %>%  # Consider positives as 1 no matter how many records in the cell
#   as.data.frame(xy = TRUE)  # Convert back to data frame
# n <- point_data %>%
#   filter(Result == "negative") %>%
#   rasterize(y = poly_raster, fun = "max") %>%  # Consider negatives as 1 no matter how many records in the cell
#   as.data.frame(xy = TRUE) %>%  # Convert back to data frame
#   anti_join(p)
# 
# nrow(p)
# nrow(n)
# # 22 vs 15 points at res = 0.001 deg
# # 22 vs 15 points at res = 0.0005 deg
# # Crashes at res = 0.00001 deg

#######################
# Assemble covariates #
#######################

train_points <- st_as_sf(bind_rows(p, n), 
                         coords = c("Longitude", "Latitude"),
                         crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

env_elev <- geodata::elevation_global(res = 0.5, path = "data/env_vars/elevation/geodata_global_res0.5")
env_land <- terra::rast("data/env_vars/landcover/landcover_type1_full_raster.tif")
env_cost <- read.csv("data/env_vars/dist_to_coast/dist_coast_sp.csv")
env_wetl <- read.csv("data/env_vars/lakes_and_wetlands/dist_wetland_sp.csv")
env_chik <- terra::rast("data/env_vars/livestock/chickens/6_Ch_2015_Aw.tif")
env_duck <- terra::rast("data/env_vars/livestock/ducks/6_Dk_2015_Aw.tif")


# Extract values for non-temporal covariates using simple point method:

train_points <- terra::extract(env_elev, train_points, method = "simple", bind = TRUE)
train_points <- terra::extract(env_land, train_points, method = "simple", bind = TRUE)
train_points$dist_to_coast <- env_cost$dist_to_coast
train_points$dist_to_wetland <- env_wetl$dist_to_wetland

######### CHECK SEARCH RADIUS HERE AND BIND ONLY SINGLE VALUE
train_points <- terra::extract(env_chik, train_points, method = "simple", bind = TRUE, search_radius = 5000) # gives some NAs - may be a way of asking for nearest cell
train_points <- terra::extract(env_duck, train_points, method = "simple", bind = TRUE, search_radius = 5000) # gives some NAs - may be a way of asking for nearest cell


# Extract values for temporal covariates by matching raster layer of appropriate calendar date

# NDVI (16-day periods)
env_ndvi <- terra::rast("data/env_vars/ndvi/ndvi_lat_long.tif")
ndvi_layer_index <- sapply(train_points$Date, match_raster_layer, raster_ymd = env_ndvi)

train_points$ndvi <- mapply(
  function(i, layer){
    extract(env_ndvi[[layer]], train_points[i, ], search_radius = 2000)[1, 2]   # use search radius to ensure coastal cells are assigned value of nearest-neighbour
  }, 
  i = seq_len(nrow(train_points)), layer = ndvi_layer_index
)


# Mean relative humidity (daily)
env_humd <- terra::rast("data/env_vars/humidity/humidity_daily.tif")
humd_layer_index <- sapply(train_points$Date, match_raster_layer, raster_ymd = env_humd)
train_points$humd <- mapply(
  function(i, layer){
    extract(env_humd[[layer]], train_points[i, ], search_radius = 2000)[1, 2]   # use search radius to ensure coastal cells are assigned value of nearest-neighbour
  }, 
  i = seq_len(nrow(train_points)), layer = humd_layer_index
)

# Mean relative humidity in the month prior
train_points$humd_m <- mapply(
  function(i, layer){
    env_humd[[calendar_wraparound_month(layer)]] %>% mean %>% extract(., train_points[i, ], search_radius = 2000) %>% .[1, 2]   # use search radius to ensure coastal cells are assigned value of nearest-neighbour
  }, 
  i = seq_len(nrow(train_points)), layer = humd_layer_index
)


# Zero-degree isotherm (daily)
env_ziso <- terra::rast("data/env_vars/isotherm/isocline_raw.grib")
ziso_layer_index <- sapply(train_points$Date, match_raster_layer, raster_ymd = env_ziso)
train_points$ziso <- mapply(
  function(i, layer){
    extract(env_ziso[[layer]], train_points[i, ], search_radius = 2000)[1, 2]   # use search radius to ensure coastal cells are assigned value of nearest-neighbour
  }, 
  i = seq_len(nrow(train_points)), layer = ziso_layer_index
)

# Minimum zero-degree isotherm in the month prior
train_points$ziso_m <- mapply(
  function(i, layer){
    env_ziso[[calendar_wraparound_month(layer)]] %>% min %>% extract(., train_points[i, ], search_radius = 2000) %>% .[1, 2]   # use search radius to ensure coastal cells are assigned value of nearest-neighbour
  }, 
  i = seq_len(nrow(train_points)), layer = ziso_layer_index
)


# Mean middday temperature (daily)
env_temp <- terra::rast("data/env_vars/climate/temp_midday_daily.tif")
temp_layer_index <- sapply(train_points$Date, match_raster_layer, raster_ymd = env_temp)
train_points$temp <- mapply(
  function(i, layer){
    extract(env_temp[[layer]], train_points[i, ], search_radius = 2000, na.rm=TRUE)[1, 2]   # use search radius to ensure coastal cells are assigned value of nearest-neighbour
  }, 
  i = seq_len(nrow(train_points)), layer = temp_layer_index
)

# Mean midday temperature in the month prior
train_points$temp_m <- mapply(
  function(i, layer){
    env_temp[[calendar_wraparound_month(layer)]] %>% mean %>% extract(., train_points[i, ], search_radius = 2000, na.rm=TRUE) %>% .[1, 2]   # use search radius to ensure coastal cells are assigned value of nearest-neighbour
  }, 
  i = seq_len(nrow(train_points)), layer = temp_layer_index
)


# Diurnal tdrnerature range (daily)
env_diur <- terra::rast("data/env_vars/climate/temp_diurnal_daily.tif")
diur_layer_index <- sapply(train_points$Date, match_raster_layer, raster_ymd = env_diur)
train_points$diur <- mapply(
  function(i, layer){
    extract(env_diur[[layer]], train_points[i, ], search_radius = 2000, na.rm=TRUE)[1, 2]   # use search radius to ensure coastal cells are assigned value of nearest-neighbour
  }, 
  i = seq_len(nrow(train_points)), layer = diur_layer_index
)

# Mean diurnal temperature range in the month prior
train_points$diur_m <- mapply(
  function(i, layer){
    env_diur[[calendar_wraparound_month(layer)]] %>% mean %>% extract(., train_points[i, ], search_radius = 2000, na.rm=TRUE) %>% .[1, 2]   # use search radius to ensure coastal cells are assigned value of nearest-neighbour
  }, 
  i = seq_len(nrow(train_points)), layer = diur_layer_index
)


# Precipitation (daily)
env_prec <- terra::rast("data/env_vars/climate/total_precip_daily.tif")
prec_layer_index <- sapply(train_points$Date, match_raster_layer, raster_ymd = env_prec)
train_points$prec <- mapply(
  function(i, layer){
    extract(env_prec[[layer]], train_points[i, ], search_radius = 2000, na.rm=TRUE)[1, 2]   # use search radius to ensure coastal cells are assigned value of nearest-neighbour
  }, 
  i = seq_len(nrow(train_points)), layer = prec_layer_index
)

# Total precipiation in the month prior
train_points$prec_m <- mapply(
  function(i, layer){
    env_prec[[calendar_wraparound_month(layer)]] %>% sum %>% extract(., train_points[i, ], search_radius = 2000, na.rm=TRUE) %>% .[1, 2]   # use search radius to ensure coastal cells are assigned value of nearest-neighbour
  }, 
  i = seq_len(nrow(train_points)), layer = prec_layer_index
)



# CHECK NAS
ggplot() +
  geom_spatraster(data = env_chik, maxcell = 5000) +
  geom_sf(data = train_points[217,], size = 4) +
  theme_minimal() +
  coord_sf(xlim = c(31.5, 32), ylim = c(31,31.5)) +
  labs(x = "Longitude",
       y = "Latitude",
       shape = "Source", color = "Result")

ggplot() +
  geom_sf(data = egy_map[1], fill = "grey", color = "black", alpha = 0.4)  +
  geom_sf(data = point_data %>% arrange(Result, Source), aes(color = Result, shape = Source), size = 4, alpha = 0.7) +
  theme_minimal() +
  coord_sf(xlim = c(30,34), ylim = c(29,32)) +
  scale_shape_manual(values=c(19, 18, 15, 17)) +
  scale_color_manual(values=c('#56B4E9','#F8766D')) +
  labs(x = "Longitude",
       y = "Latitude",
       shape = "Source",
       color = "Result")

ggplot() +
  #  geom_sf(data = egy_map[1], fill = "grey", color = "black", alpha = 0.4)  +
  geom_spatraster(data = env_elev, maxcell = 500) +
  geom_sf(data = train_points %>% st_as_sf, aes(shape = Result), size = 4, color = "grey70") +
  theme_minimal() +
  guides(fill=guide_legend(title="Elevation")) +
  coord_sf(xlim = c(25,36.5), ylim = c(22,32))

ggplot() +
  geom_sf(data = egy_map[1], fill = "grey", color = "black", alpha = 0.4)  +
  geom_sf(data = train_points %>% st_as_sf, aes(shape = Result, color = ziso), size = 4) +
  theme_minimal() +
  coord_sf(xlim = c(25,36.5), ylim = c(22,32))

# Examine NAs - what is overlap with raster where covar values cannot be extracted
ggplot() +
  geom_spatraster(data = env_duck) +
  geom_sf(data = egy_map[1], color = "black", alpha = 0.05)  +
  geom_sf(data = train_points %>% st_as_sf, aes(shape = Result, color = "red"), size = 4) +
  theme_minimal() +
  coord_sf(xlim = c(25,36.5), ylim = c(22,32))

# View raster - crashes when setting coord_sf?
ggplot() +
  geom_spatraster(data = env_land, maxcell = 5e+03) +
  geom_sf(data = egy_map[1], color = "black", alpha = 0.05)  +
  guides(fill=guide_legend(title="Land cover category")) +
  theme_minimal()  +
  coord_sf(xlim = c(25,36.5), ylim = c(22,32))
