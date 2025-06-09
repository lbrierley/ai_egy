## Adapted from script by Hayes et al., date 25/1/24

library(tidyverse)
library(terra)
library(sf)

# # Define scale of raster to calculate if calculating on raster
# egy_rast <- terra::rast(extent=terra::ext(24,37,21,32), res = 0.25)

## need a bigger map to avoid the boundary being classified as the sea
afr_map <- rnaturalearth::ne_download(scale = 10, type = "countries") %>% subset(., REGION_WB %in% c("Middle East & North Africa",  "Sub-Saharan Africa"))
plot(afr_map[1])

# However, we just want the outline so have to remove the internal lines

sf_use_s2(FALSE)
afrunion <- st_union(afr_map$geom)
afrunion <- nngeo::st_remove_holes(afrunion)
sf_use_s2(TRUE)
afrline <- st_cast(afrunion, "MULTILINESTRING")

# Calculate distance to coast using the line
## use intersects instead of intersection

##### Sample point version 
dist <- st_distance(afrline, train_points)

dist_coast_sp <- data.frame(dist = as.vector(dist)/1000,
                              st_coordinates(train_points)) %>%
  rename(dist_to_coast = dist)

write.csv(dist_coast_sp, "data/env_vars/dist_to_coast/dist_coast_sp.csv")

# ##### Raster version 
# egy_points <- raster::xyFromCell(egy_rast, 1:ncell(egy_rast)) %>% as.data.frame %>% st_as_sf(coords = c("x", "y"), crs = "WGS84")
# 
# afrgrid <- st_intersects(egy_points, 
#                          afrunion, 
#                          sparse = F)
#
# egy_points_land <- raster::xyFromCell(egy_rast[afrgrid,], 1:ncell(egy_rast)) %>% as.data.frame %>% st_as_sf(coords = c("x", "y"), crs = "WGS84")
# 
# # Calculate distance to coast
# dist <- st_distance(afrline, egy_points_land)
# 
# 
# dist_coast_rast <- data.frame(dist = as.vector(dist)/1000,
#                               st_coordinates(egy_points_land)) %>%
#   select(X, Y, dist) %>%
#   rast(type="xyz")
# 
# terra::writeRaster(dist_coast_rast, "data/env_vars/dist_to_coast/dist_to_coast.tif", overwrite=TRUE)

