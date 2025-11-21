rm(list = ls())

library(caret)
library(dbarts)
library(embarcadero)
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(rnaturalearth)

# Define and resave final model for ease
final_model_path <- "models/ROSE/reduced_opt_model_2.rds"

load(file = final_model_path)

saveRDS(sdm, file = "models/final_model.rds")


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

# Generate predictions for each seasonal date
for(i in 1:length(pred_dates)){
  
  # Load in each covariate set
  pred_df <- read.csv(paste0("data/pred/pred_grid_",i,".csv")) %>% 
    filter(complete.cases(.)) %>%
    dplyr::select(-X)
  
  pred_vect <- vect(pred_df, geom=c("x","y"), crs="EPSG:4326")
  
  # NOTE: MUST be a RasterStack type object for predictions, not a SpatRast!
  pred_rast <- rasterize(pred_vect, egy_rast, field = names(pred_vect)) %>% stack()
  
  # Generate risk map with percentiles
  pred_layers <- predict2.bart(object = sdm,
                               x.layers = pred_rast,           
                               quantiles = c(0.025, 0.975),
  )
  
  names(pred_layers) <- c("mean",
                          "lower_cb",
                          "upper_cb")
  
  pred_layers %>% save(file = paste0("models/predictions_", i, ".rds"))
  pred_layers %>% as.data.frame(xy = TRUE) %>% mutate(date = pred_dates[i]) %>% write.csv(file = paste0("models/predictions_df_", i, ".csv"))
  # load(file = "models/predictions.rds")
}

pred_all <- list.files(path="models/", pattern = "\\.csv", full.names=TRUE) %>% 
  map_dfr(read.csv)


# Plot mean predictions
p1 <- ggplot() +
  geom_raster(data = pred_all, aes(x = x, y = y, fill = mean)) +
  geom_sf(data = egy_map[1], color = "black", alpha = 0)  +
  theme_bw() +
  scale_fill_viridis_c(option = "plasma", na.value = "transparent") +
  facet_wrap(date ~ ., nrow = 1) +
  labs(x = "Longitude",
       y = "Latitude",
       fill = "prob.")

ggsave("figures/mean_predictions.png", plot = p1, width = 14, height = 4)

p1_red <- ggplot() +
  geom_raster(data = pred_all, aes(x = x, y = y, fill = mean)) +
  geom_sf(data = egy_map[1], color = "black", alpha = 0)  +
  theme_bw() +
  scale_fill_viridis_c(option = "plasma", na.value = "transparent") +
  facet_wrap(date ~ ., nrow = 1) +
  coord_sf(xlim = c(30,34), ylim = c(29,32)) +
  labs(x = "Longitude",
       y = "Latitude",
       fill = "prob.")

ggsave("figures/mean_predictions_nile_delta.png", plot = p1_red, width = 14, height = 4)

# Plot uncertainty
for(i in 1:length(pred_dates)){

u <- ggplot() +
  geom_raster(data = pred_all %>% 
                filter(date == pred_dates[i]) %>% 
                pivot_longer(cols = mean:upper_cb, names_to = "var"), 
              aes(x = x, y = y, fill = value)) +
  geom_sf(data = egy_map[1], color = "black", alpha = 0)  +
  theme_bw() +
  scale_fill_viridis_c(option = "plasma", na.value = "transparent") +
  facet_wrap(var ~ ., nrow = 1) +
  labs(x = "Longitude",
       y = "Latitude",
       fill = "prob.")

ggsave(paste0("figures/uncertainty_seas_",i,".png"), plot = u, width = 10, height = 4)

}
