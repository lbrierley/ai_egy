## Adapted from script by Hayes et al., date 8/8/23

library(luna)

# Retrieve NDVI data
egy_map <- ne_countries(country = "Egypt", returnclass = "sf")

tile_list <- luna::getNASA(
  product = "MOD13Q1", 
  start_date = "2022-01-01", 
  end_date = "2022-12-31", 
  aoi = st_bbox(egy_map),
  version = "061",
  download = FALSE
)

tile_data <- luna::getNASA(
  product = "MOD13Q1", 
  start_date = "2022-01-01", 
  end_date = "2022-12-31", 
  aoi = st_bbox(egy_map),
  version = "061",
  download = TRUE, 
  path = "data/env_vars/ndvi/raw",
  username="REPLACE WITH YOUR OWN USER EMAIL",
  password="REPLACE WITH YOUR OWN USER PASSWORD"
)

# Check tiles
terra::rast(list.files("data/env_vars/ndvi/raw", full.names = TRUE)[1])

# Keep only first layer
for (i in 1:length(list.files("data/env_vars/ndvi/raw", pattern=".hdf", full.names = TRUE))) {
  rr <- terra::rast(list.files("data/env_vars/ndvi/raw", pattern=".hdf", full.names = TRUE)[[i]], lyrs = "\"250m 16 days NDVI\"")
  terra::writeRaster(rr, paste0("data/env_vars/ndvi/raw/", sprintf("%02d", i), ".tif"), overwrite = T)
}


# Extract year and calendar day of NDVI records
ndvi_names <- list.files(path = "data/env_vars/ndvi/raw/", pattern = "*.hdf$") %>% 
  gsub("\\.hdf", "", .) %>% 
  gsub(".*\\.","", .) %>%
  substr(., start = 1, stop = 7) %>%
  unique

for (i in 1:length(ndvi_names)){
  ## combine tiles into rasters
  vrt(
    x = list.files(path = "data/env_vars/ndvi/raw/", pattern = "*.tif$", full.names = TRUE)[(4*i-3):(4*i)], 
    filename = "data/env_vars/ndvi/ndvi.vrt", overwrite = T
  )
  
  nrast <- rast("data/env_vars/ndvi/ndvi.vrt")
  terra::writeRaster(nrast * 0.00000001, paste0("data/env_vars/ndvi/ndvi_", ndvi_names[i],".tif"), overwrite=TRUE) # Change scale
}

## Combine into single raster, reproject to lat-long, and assign times
ndvi_all <- terra::rast(list.files(path = "data/env_vars/ndvi/", pattern = "*.tif$", full.names = TRUE)) %>% 
  terra::project(y = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
time(ndvi_all) <- ndvi_names %>% as.Date(format = "%Y%j")

ndvi_all %>% terra::writeRaster("data/env_vars/ndvi/ndvi_lat_long.tif")
