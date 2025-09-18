## Adapted from script by Hayes et al., date 9/6/23

glwd_rast <- rast("data/env_vars/lakes_and_wetlands/glwd_3/w001001x.adf")
glwd_rast # resolution is supposed to be 30 arc seconds = 0.008 degrees = approx 1km

# Crop to reasonable boundaries to capture cross-continent wetlands
glwd_crop_prj <- terra::crop(x = glwd_rast, y = terra::ext(10,50,10,50))
plot(glwd_crop_prj) 
head(glwd_crop_prj)

# info on the key is available in the data documentation pdf in the data file
# Not sure we want to differentiate between the different types? Re-code so that they are all 
# either water or not-water?? 

table(values(glwd_crop_prj))
# This has some values of 9 in it. We don't want to include these as they are intermittent wetland
# We just want to include permanent wetlands.

glwd_crop_prj_combo <- terra::subst(glwd_crop_prj, 
                                    from = c(1:8), 
                                    to = c(rep(99,8)),  # recode category 1:8 as same category (99)
                                    others = NA)        # removes category 9
plot(glwd_crop_prj_combo)


# Calculate distance to nearest wetland

##### Sample point version 
dist_wetl <- crds(train_points) %>% as.data.frame
dist_wetl$dist_to_wetland <- as.vector(distance(train_points, as.polygons(glwd_crop_prj_combo))) # convert to vector and calc point-vector distance

write.csv(dist_wetl, "data/env_vars/lakes_and_wetlands/dist_wetland_sp.csv")

##### Raster version - may take a long time
# To points for distance calculation to centroid
egy_points <- raster::xyFromCell(egy_rast, 1:ncell(egy_rast)) %>% as.data.frame %>% st_as_sf(coords = c("x", "y"), crs = "WGS84")

sf_use_s2(FALSE)
egy_filt <- st_intersects(egy_points, egy_map, sparse = F) # returns TRUE if on land polygon, FALSE if in ocean
sf_use_s2(TRUE)

dist_wetl <- crds(vect(egy_points[egy_filt,])) %>% as.data.frame
dist_wetl$dist <- distance(vect(egy_points[egy_filt,]), as.polygons(glwd_crop_prj_combo)) # convert to vector and calc point-vector distance

# Rasterize again
dist_wetland_rast <- data.frame(dist_wetl) %>%
  select(x, y, dist) %>%
  rast(type="xyz")

terra::writeRaster(dist_wetland_rast, "data/env_vars/lakes_and_wetlands/dist_wetland.tif", overwrite=TRUE)

