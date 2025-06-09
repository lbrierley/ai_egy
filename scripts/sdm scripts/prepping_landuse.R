## Adapted from script by Hayes et al., date 9/8/23

#install.packages("remotes")
#remotes::install_github("rspatial/luna")

#install.packages('luna', repos='https://rspatial.r-universe.dev')


rm(list = ls())

library(luna)

# Retrieve land use data
egy_map <- ne_download(scale = 10, type = "countries") %>% subset(., NAME_EN == "Egypt")

tile_list <- luna::getNASA(
  product = "MCD12Q1", 
  start_date = "2022-01-01", 
  end_date = "2022-12-31", 
  aoi = st_bbox(egy_map),
  version = "061",
  download = FALSE
)

tile_data <- luna::getNASA(
  product = "MCD12Q1", 
  start_date = "2022-01-01", 
  end_date = "2022-12-31", 
  aoi = st_bbox(egy_map),
  version = "061",
  download = TRUE, 
  path = "data/env_vars/landcover",
  username="sarahhayes",
  password="NASATigtogs43!"
)

# Check tiles
terra::rast(list.files("data/env_vars/landcover", full.names = TRUE)[1])

# Keep only first layer
for (i in 1:length(list.files("data/env_vars/landcover", pattern=".hdf", full.names = TRUE))) {
  rr <- terra::rast(list.files("data/env_vars/landcover", pattern=".hdf", full.names = TRUE)[[i]], lyrs = "LC_Type1")
  terra::writeRaster(rr, paste0("data/env_vars/landcover/processing/",i,".tif"), overwrite = T)
}

## combine them into a single raster
vrt(
  x = list.files(path = "data/env_vars/landcover/processing/",
                 pattern = "*.tif$", full.names = TRUE), 
  filename = "dem.vrt", overwrite = T
)

dem <- rast("dem.vrt")
dem_fact <- terra::as.factor(dem)
terra::writeRaster(dem_fact, "data/env_vars/landcover/landcover_type1_full_raster.tif", overwrite=TRUE)

