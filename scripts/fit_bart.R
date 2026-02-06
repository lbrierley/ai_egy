# In this script we train a BART model to classify sites for avian flu
# presence/absence based on environmental and species abundance factors.
# This script includes cross-validation of BART parameters.

rm(list=ls())

# Global storing optimised k value for BART - which may be updated
K_OPT <- 2

# Number of chains/threads for MCMC run. The code is set up to do one thread per
# chains
N_CHAINS <- 4L

# Ratio of upweighting positive samples to negative
pos_weight <- 2

library(caret)
library(dbarts)
library(dplyr)
library(embarcadero)
library(pROC)
library(raster)
library(reshape)
library(rsample)
library(terra)
library(tidyverse)

# Helper/adapted functions
source("scripts/bart_scripts/bart_functions.R")

# Read in data
df <- read.csv("data/full_env_data.csv") %>% 
  dplyr::select(-X) %>%
  #mutate(across(starts_with("dem"), as.factor)) %>%
  mutate(site = paste(x,y,sep="_")) %>% 
  mutate(Result = as.numeric(as.factor(Result))-1) # convert to 0/1 binary for dbarts

# Determine predictor columns to include

preds <- c("wc2.1_30s_elev",
           "dist_to_coast", 
           "dist_to_wetland",
           "X6_Ch_2015_Aw",
           "X6_Dk_2015_Aw",
           # "dem1", 
           # "dem2", 
           # "dem4", 
           # "dem5", 
           # "dem6", 
           # "dem7", 
           # "dem8", 
           "dem9",                   # only landcover categories 9-12 and 16-17 are represented in the data
           "dem10", 
           "dem11",
           "dem12",
           "dem13",
           # "dem14",
           # "dem15",
           "dem16",
           "dem17", 
           "ndvi",
           "humd",
           "humd_m",
           "ziso",
           "ziso_m",
           "temp",
           "temp_m",
           "diur",
           "diur_m",
           "prec",
           "prec_m"
)

# Predictors which are time-invariant (i.e., constant within each of the 46 geolocated sites)
const_preds <- c("wc2.1_30s_elev",
                 "dist_to_coast", 
                 "dist_to_wetland",
                 "X6_Ch_2015_Aw",
                 "X6_Dk_2015_Aw",
                 # "dem1", 
                 # "dem2", 
                 # "dem4", 
                 # "dem5", 
                 # "dem6", 
                 # "dem7", 
                 # "dem8", 
                 "dem9",                   # only landcover categories 9-12 and 16-17 are represented in the data
                 "dem10", 
                 "dem11",
                 "dem12",
                 "dem13",
                 # "dem14",
                 # "dem15",
                 "dem16",
                 "dem17"
)

daily_preds <- c("humd",
                 "ziso",
                 "temp",
                 "diur",
                 "prec"
)

# preds <- preds[!(preds %in% const_preds)] # remove constants per site
preds <- preds[!(preds %in% daily_preds)] # remove daily covars (and keep monthly covars instead)


# Cluster some sites together based on close proximity - many covariates at the 0.1 deg level, so consider same if within Euclidean distance of 0.15 deg (roughly diagonal distance of a 0.1*0.1 grid cell)

raw_sites <- df %>% dplyr::select(x,y) %>% distinct %>% arrange(x)

dist_sites <- dist(raw_sites, method = "euclidean", diag = T) %>% round(2)
tree_sites <- hclust(dist_sites, method = "complete")
clust_sites <- cutree(tree_sites, h = 0.15)

raw_sites <- raw_sites %>% bind_cols(clustsite = as.factor(clust_sites))

# library(sf)
# library(rnaturalearth)
# egy_map <- ne_countries(country = "Egypt", returnclass = "sf")
# 
# ggplot() +
#   geom_sf(data = egy_map[1], fill = "grey", color = "black", alpha = 0.4)  +
#   geom_text(data = raw_sites, aes(x = x, y = y, label = clustsite, color = as.factor(clustsite))) +
#   guides(color="none") +
#   theme_minimal() +
#   coord_sf(xlim = c(27,34), ylim = c(23.4,32))
# 
# ggplot() +
#   geom_sf(data = egy_map[1], fill = "grey", color = "black", alpha = 0.4)  +
#   geom_text(data = raw_sites, aes(x = x, y = y, label = clustsite, color = as.factor(clustsite))) +
#   guides(color="none") +
#   theme_minimal() +
#   coord_sf(xlim = c(30,34), ylim = c(29,32))

df <- df %>% left_join(raw_sites, by = c("x","y"))

# Define training/test split by random sampling

set.seed(1100)

# # Ordinary training:test split
# split <- rsample::initial_split(df, prop = 3/4)
# Stratified training:test split by clustered site
split <- rsample::group_initial_split(df, prop = 3/4, group = clustsite)

#ROSE to oversample training data
train_temp <- as.data.frame(training(split))
train_resamp <- ROSE::ovun.sample(Result ~ ., data = train_temp,
                                  seed = 1726,
                                  N = 500,
                                  method="over")$data


# # Consider MN data the validation dataset
# MN_sites <- c("32.3061111111111_31.2625", "32.0333333333333_31.1833333333333", "31.8213888888889_31.4166666666667")
# test_MN <- df %>% filter(site %in% MN_sites)
# train_MN <- df %>% filter(!(site %in% MN_sites))

# Define training and test sets
main_train <- train_resamp
main_test <- testing(split)

# Create inner folds for 1 x 5-fold cross-validation
set.seed(1108)
fold_ids <- caret::groupKFold(main_train$clustsite, k = 5)

folds <- lapply(1:length(fold_ids),
                FUN = function(i){
                  main_train[-fold_ids[[i]],] })
# Complements to the folds, i.e. test sets for model trained on that fold
antifolds <- lapply(1:length(fold_ids),
                    FUN = function(i){
                      main_train[fold_ids[[i]],]})

# Now cycle over possible parameter values
k_vals = c(1, 2, 3)
power_vals = c(1.6, 1.8, 2)
base_vals = c(0.75, 0.85, 0.95)
kl <- length(k_vals)
pl <- length(power_vals)
bl <- length(base_vals)
cv_results <- data.frame(k=numeric(),
                         power=numeric(),
                         base=numeric(),
                         auc1=numeric(),
                         auc2=numeric(),
                         auc3=numeric(),
                         auc4=numeric(),
                         auc5=numeric())
cv_results[1:kl*pl*bl, ] <- 0

for (i in 1:length(k_vals)){
  k_val <- k_vals[i]
  for (j in 1:length(power_vals)){
    power_val <- power_vals[j]
    for (m in 1:length(base_vals)){
      base_val <- base_vals[m]
      idx <- (i-1)*pl*bl + (j-1)*bl + m
      # cat("Crossvalidating parameter set ", idx, " of ", kl*pl*bl,"\n") # Uncomment to print where we are in the cycle
      cv_results$k[idx] <- k_val
      cv_results$power[idx] <- power_val
      cv_results$base[idx] <- base_val
      for (fold_no in 1:length(folds)){ # For each fold, fit a model and calculate test AUC
        model <- bart2(formula = reformulate(preds, response = "Result"),
                       data = folds[[fold_no]],
                       test = antifolds[[fold_no]] %>% dplyr::select(all_of(preds)),
                       weights = recode(folds[[fold_no]]$Result, `0` = 1, `1` = pos_weight),       # Upweight positives
                       k = k_val,
                       power = power_val,
                       base = base_val,
                       n.trees = 200,    # small number of trees for tuning parameters
                       n.chains = 1L,
                       n.threads = 1L,
                       keepTrees = TRUE,
                       verbose = FALSE)
        antifold_x <- subset(antifolds[[fold_no]], select=-c(Result))
        antifold_y <- antifolds[[fold_no]]$Result
        cutoff <- get_threshold(model)
        auc <- get_sens_and_spec(model, antifold_x, antifold_y, cutoff)$auc
        cv_results[idx, 3+fold_no] <- auc
        rm(model)
      }
    }
  }
}

# Calculate fold-wise mean AUC associated with each parameter set and choose optimum parameters based on maximum AUC
cv_results <- cv_results %>%
  rowwise() %>%
  mutate(mean_auc = mean(c(auc1,
                           auc2,
                           auc3,
                           auc4,
                           auc5)))
argmax <- which.max(cv_results$mean_auc)
K_OPT <- cv_results$k[argmax]
power_opt <- cv_results$power[argmax]
base_opt <- cv_results$base[argmax]


# Retrain most optimal model, without variable selection
basic_model <- bart( 
  x.train = main_train %>% dplyr::select(all_of(preds)),
  y.train = main_train %>% pull(Result),
  x.test = main_test %>% dplyr::select(all_of(preds)),
  weights = recode(main_train$Result, `0` = 1, `1` = pos_weight),       # Upweight positives
  k = K_OPT,
  ntree = 500,
  power = power_opt,
  base = base_opt,
  nchain = N_CHAINS,
  nthread = N_CHAINS,
  keeptrees = TRUE)
invisible(basic_model$fit$state)

save(basic_model,  file = paste0("models/base_opt_model_", pos_weight, ".rds"))


# Perform variable selection
sdm <- bart.step(
  x.data = main_train %>% dplyr::select(all_of(preds)),
  y.data = main_train %>% pull(Result),
  x.weight = recode(main_train$Result, `0` = 1, `1` = pos_weight),       # Upweight positives
  k = K_OPT,
  power = power_opt,
  base = base_opt)
invisible(sdm$fit$state)

# # Interrupt point to save/load if needed
save(sdm, file = paste0("models/reduced_opt_model_", pos_weight, ".rds"))
# load(file = paste0("models/reduced_opt_model_", pos_weight, ".rds"))

preds_p <- main_test %>% 
  dplyr::select(all_of(preds)) %>%
  stats::predict(object=sdm, type = "bart") %>% 
  pnorm() %>%
  colMeans()

ROC = roc(response = main_test %>% pull(Result), 
          predictor = preds_p,
          direction = "<")

ROC

cutoff <- coords(ROC, "best", best.method="closest.topleft")$threshold

matrix_test <- confusionMatrix(data = if_else(preds_p > cutoff, 1, 0) %>% as.factor(), 
                               reference = main_test %>% pull(Result) %>% as.factor(), 
                               positive = "1")
matrix_test

test_plot <- main_test %>%
  mutate(preds_p = preds_p,
         preds_c = if_else(preds_p > cutoff, 1, 0) %>% as.factor(),
         plotkey = paste(Result, preds_c, sep = "_"))

# Performance seems to vary site-by-site!
test_plot %>% with(., table(Result, preds_c, clustsite))

line <- bind_cols(pos_weight = pos_weight,
                  threshold = cutoff,
                  matrix_test$overall %>% t(),
                  AUC = ROC$auc %>% as.numeric(),
                  matrix_test$byClass %>% t()) %>%
  mutate(across(where(is.numeric), round, 3))

line %>% write.csv(file = paste0("models/metrics_", pos_weight, ".csv"))

# metric_cis <- get_sens_and_spec_ci(sdm, 
#                                      xtest = main_test %>% dplyr::select(all_of(preds)), 
#                                      ytest = main_test %>% pull(Result),
#                                      cutoff)
#
# save(metrics, metric_cis, file = "models/metrics.rds")

# Variable importance

varimp_raw <- sdm$varcount %>% 
  as.data.frame %>%
  mutate(./rowSums(.)) %>%
  mutate(./max(.))

varimp_summ <-  bind_rows(varimp_raw %>% summarise(across(everything(), mean)),
                          varimp_raw %>% summarise(across(everything(), sd))) %>% # IQR, 95% CI?
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "var") %>%
  dplyr::rename("mean" = "V1", "sd" = "V2") %>%
  arrange(-mean)

v <- varimp_summ %>% 
  mutate(var = case_when(var == "humd_m" ~ "mean humidity (prior 30 days)",
                         var == "temp_m" ~ "mean temp. (prior 30 days)",
                         var == "dist_to_wetland" ~ "distance to wetland",
                         var == "ndvi" ~ "NDVI (16-day interval)",
                         var == "ziso_m" ~ "minimum isotherm (prior 30 days)",
                         var == "diur_m" ~ "mean diurnal range (prior 30 days)",
                         var == "prec_m" ~ "total rainfall (prior 30 days)"
                           )) %>%
  ggplot(aes(x = reorder(var, mean), y = mean)) + 
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,0.6)) +
  coord_flip() +
  theme_bw() +
  labs(x = "variable",
       y = "mean variable importance")

ggsave("figures/varimp.png", plot = v, width = 6, height = 2)
ggsave("figures/varimp.pdf", plot = v, width = 6, height = 2)
