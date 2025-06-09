# Match month-day of training data to calendar period of rasters (assuming raster date labels denote LAST day of each data period)
# Return raster to pluck value from
match_raster_layer <- function(md, raster_ymd) {
  idx <- which((time(raster_ymd) %>% format("%m-%d")) >= md)
  if (length(idx) == 0) return(length(time(env_ndvi)))  # if no match (because at the end of calendar year) return last possible raster, assume this spills into next year
  return(idx[1])
}

# For given calendar date return daily indexes of all days back to 30 days prior, wrapping around where needed
calendar_wraparound_month <- function(n) {
  indices <- ((n - 30):n) %% 365
  indices[indices == 0] <- 365
  return(seq(1:365)[indices])
}