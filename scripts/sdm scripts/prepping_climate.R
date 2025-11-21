library(ecmwfr) 

wf_set_key(user = "REPLACE WITH YOUR OWN USER EMAIL", 
           key = "REPLACE WITH YOUR OWN API TOKEN")

# Download daily midday temperature
for (j in 1:12) {   # can only request 1 month at once
  file <- wf_request(
    user = "REPLACE WITH YOUR OWN USER EMAIL",
    request  =  list(
      "dataset_short_name" = "reanalysis-era5-land",
      "product_type" = "reanalysis",
      "variable" = "2m_temperature",
      "year" = 2018,
      "month" = c(sprintf("0%d",1:9),10,11,12)[j],
      "day" = c(sprintf("0%d",1:9), 10:31),
      "time" = "12:00",
      "area" = "32/24/22/37", ## bbox for Egypt
      "format" = "netcdf",
      "target" = paste0("temp_midday",
                        c(sprintf("0%d",1:9),10,11,12)[j],
                        ".nc")
    ),
    transfer = TRUE,
    path = "data\\env_vars\\climate\\climate_raw\\temp_midday"
  )
}

# Process daily midday temperature
zipped <- list.files(
  "data\\env_vars\\climate\\climate_raw\\temp_midday", 
  pattern = "\\.zip", 
  full.names = TRUE)

for (i in 1:length(zipped)) {
  name <- tools::file_path_sans_ext(basename(zipped[i]))
  unzip(zipped[i], exdir =  "data\\env_vars\\climate\\climate_raw\\temp_midday")
  file.rename("data\\env_vars\\climate\\climate_raw\\temp_midday\\data_0.nc", 
              paste0("data\\env_vars\\climate\\climate_raw\\temp_midday\\", name, ".nc"))
}

nc <- list.files(
  "data\\env_vars\\climate\\climate_raw\\temp_midday", 
  pattern = "\\.nc", 
  full.names = TRUE) 
rast_temp <- terra::rast(nc)

# Convert Kelvin to Celsius
rast_temp <- rast_temp - 273.15

# Add calendar dates
time(rast_temp, tstep = "days") <- as.Date("2018-01-01") + 0:364

rast_temp %>% terra::writeRaster("data/env_vars/climate/temp_midday_daily.tif")



# Download hourly temperature
for (j in 1:12) {   # can only request 1 month at once
  file <- wf_request(
    user = "REPLACE WITH YOUR OWN USER EMAIL",
    request  =  list(
      "dataset_short_name" = "reanalysis-era5-land",
      "product_type" = "reanalysis",
      "variable" = "2m_temperature",
      "year" = 2018,
      "month" = c(sprintf("0%d",1:9),10,11,12)[j],
      "day" = c(sprintf("0%d",1:9), 10:31),
      "time" = c(paste0("0",0:9,":00"),paste0(10:23,":00")),
      "area" = "32/24/22/37", ## bbox for Egypt
      "format" = "netcdf",
      "target" = paste0("temp_all",
                        c(sprintf("0%d",1:9),10,11,12)[j],
                        ".nc")
    ),
    transfer = TRUE,
    path = "data\\env_vars\\climate\\climate_raw\\temp_all"
  )
}

# Process hourly temperature and calc diurnal range
zipped <- list.files(
  "data\\env_vars\\climate\\climate_raw\\temp_all", 
  pattern = "\\.zip", 
  full.names = TRUE)

for (i in 1:length(zipped)) {
  name <- tools::file_path_sans_ext(basename(zipped[i]))
  unzip(zipped[i], exdir =  "data\\env_vars\\climate\\climate_raw\\temp_all")
  file.rename("data\\env_vars\\climate\\climate_raw\\temp_all\\data_0.nc", 
              paste0("data\\env_vars\\climate\\climate_raw\\temp_all\\", name, ".nc"))
}

nc <- list.files(
  "data\\env_vars\\climate\\climate_raw\\temp_all", 
  pattern = "\\.nc", 
  full.names = TRUE) 
rast_temp <- terra::rast(nc)

# Calculate diurnal range (range over max - min of day)
days <- split(1:8760, ceiling((1:8760) / 24)) # assign each layer to its calendar day
rast_daily <- lapply(days, function(x) {
  app(rast_temp[[x]], fun = function(y) max(y, na.rm = TRUE) - min(y, na.rm = TRUE))
})

rast_daily <- rast(rast_daily)

# Add calendar dates
time(rast_daily, tstep = "days") <- as.Date("2018-01-01") + 0:364

rast_daily %>% terra::writeRaster("data/env_vars/climate/temp_diurnal_daily.tif")



# Download hourly precipitation
for (j in 1:12) {   # can only request 1 month at once
  file <- wf_request(
    user = "REPLACE WITH YOUR OWN USER EMAIL",
    request  =  list(
      "dataset_short_name" = "reanalysis-era5-land",
      "product_type" = "reanalysis",
      "variable" = "total_precipitation",
      "year" = 2018,
      "month" = c(sprintf("0%d",1:9),10,11,12)[j],
      "day" = c(sprintf("0%d",1:9), 10:31),
      "time" = c(paste0("0",0:9,":00"),paste0(10:23,":00")),
      "area" = "32/24/22/37", ## bbox for Egypt
      "format" = "netcdf",
      "target" = paste0("total_precip_",
                        c(sprintf("0%d",1:9),10,11,12)[j],
                        ".nc")
    ),
    transfer = TRUE,
    path = "data\\env_vars\\climate\\climate_raw\\total_precip"
  )
}

# Process hourly precipitation and calc daily total
zipped <- list.files(
  "data\\env_vars\\climate\\climate_raw\\total_precip", 
  pattern = "\\.zip", 
  full.names = TRUE)

for (i in 1:length(zipped)) {
  name <- tools::file_path_sans_ext(basename(zipped[i]))
  unzip(zipped[i], exdir =  "data\\env_vars\\climate\\climate_raw\\total_precip")
  file.rename("data\\env_vars\\climate\\climate_raw\\total_precip\\data_0.nc", 
              paste0("data\\env_vars\\climate\\climate_raw\\total_precip\\", name, ".nc"))
}

nc <- list.files(
  "data\\env_vars\\climate\\climate_raw\\total_precip", 
  pattern = "\\.nc", 
  full.names = TRUE) 
rast_temp <- terra::rast(nc)

# Calculate daily total precipitation (sum)
days <- split(1:8760, ceiling((1:8760) / 24)) # assign each layer to its calendar day
rast_daily <- lapply(days, function(x) {
  sum(rast_temp[[x]], na.rm = TRUE)
})

rast_daily <- rast(rast_daily)

# Add calendar dates
time(rast_daily, tstep = "days") <- as.Date("2018-01-01") + 0:364

rast_daily %>% terra::writeRaster("data/env_vars/climate/total_precip_daily.tif")
