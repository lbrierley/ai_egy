rm(list = ls())

library(readxl)
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(rnaturalearth)

## Bring in the different data sets and combine to one large data set 

# MN data
mn_data <- read_excel("data/raw_data/wild bird sampling 2019-2023 collective.xlsx") %>%
  select(Bird, Date, Coordinates, Result, H5, H9) %>%
  rename(Species = Bird) %>%
  mutate(Source = "MN") %>% 
  mutate(Coordinates = gsub("[N|E|S|W|″]", "",Coordinates)) %>%   # Convert character coordinates to decimal
  separate(Coordinates, c("Lat", "Long"), sep=" ") %>%
  mutate(Lat = gsub("[°|′]", " ", Lat)) %>%
  mutate(Long = gsub("[°|′]", " ", Long)) %>%
  separate(Lat, c("lat_d", "lat_m", "lat_s"), sep = " ") %>%
  separate(Long, c("long_d", "long_m", "long_s"), sep = " ") %>%
  mutate(across(lat_d:long_s, ~as.numeric(.))) %>%
  mutate(across(lat_m:lat_s, ~coalesce(., 0))) %>%
  mutate(across(long_m:long_s, ~coalesce(., 0))) %>%
  mutate(Latitude = lat_d + lat_m/60 + lat_s/60^2,
         Longitude = long_d + long_m/60 + long_s/60^2) %>%
  select(Species, Date, Result, H5, H9, Source, Latitude, Longitude)

# FAO 
fao_data  <- read.csv("data/raw_data/epidemiology-raw-data_wild_apr_2024.csv") %>%
  rename(observation.date = "Observation.date..dd.mm.yyyy.") %>%
  rename(report.date = "Report.date..dd.mm.yyyy.") %>%
  filter(Country == "Egypt") %>%
  filter(observation.date != "") %>% # There are rows with no entries for observation date so remove these
  mutate(H5 = case_when(grepl("H5", Serotype) ~ "positive",
                        is.na(Serotype) ~ NA,
                        Serotype %in% c("HPAI", "LPAI", "") ~ NA,
                        TRUE ~ "negative"),
         H9 = case_when(grepl("H9", Serotype) ~ "positive",
                        is.na(Serotype) ~ NA,
                        Serotype %in% c("HPAI", "LPAI", "") ~ NA,
                        TRUE ~ "negative")) %>%
  select(all_of(c("Latitude", "Longitude", "observation.date", "Species", "H5", "H9"))) %>%
  mutate(Source = "FAO") %>%
  mutate(Result = "positive") %>% # only positive records kept in this dataset
  mutate(observation.date = as.Date(observation.date, "%d/%m/%Y")) %>%
  rename(Date = observation.date)

# WAHIS
wahis <- read_excel("data/raw_data/infur_20250120.xlsx", sheet = 2, guess_max = 160000) %>% 
  filter(disease_eng %in% 
           c("Low pathogenic avian influenza (poultry) (2006-2021)",
             "High pathogenicity avian influenza viruses (poultry) (Inf. with)",
             "Influenza A virus (Inf. with)",
             "Influenza A viruses of high pathogenicity (Inf. with) (non-poultry including wild birds) (2017-)"))

wahis_data <- wahis %>% 
  filter(is_wild == TRUE & (wild_type != "captive"|is.na(wild_type))) %>%
  mutate(H5 = case_when(grepl("H5", sero_sub_genotype_eng) ~ "positive",
                        is.na(sero_sub_genotype_eng) ~ NA,
                        sero_sub_genotype_eng %in% c("HPAI", "LPAI", "") ~ NA,
                        TRUE ~ "negative"),
         H9 = case_when(grepl("H9", sero_sub_genotype_eng) ~ "positive",
                        is.na(sero_sub_genotype_eng) ~ NA,
                        sero_sub_genotype_eng %in% c("HPAI", "LPAI", "") ~ NA,
                        TRUE ~ "negative")) %>%
  filter(country == "Egypt") %>%
  select(all_of(c("Latitude", "Longitude", "Outbreak_start_date", "Species", "H5", "H9"))) %>%
  mutate(Source = "WOAH") %>%
  mutate(Result = "positive") %>% # only positive records kept in this dataset
  rename(Date = Outbreak_start_date)

rm(wahis)
gc()

# BVBRC
bvbrc_data <- read.csv("data/raw_data/BVBRC_surveillance.csv") %>%
  filter(Host.Natural.State == "Wild") %>%
  filter(!(is.na(Collection.Year))) %>%
  filter(Host.Species != "Env") %>% # There are samples labelled "env" which I assume are environmental so remove these
  filter(Collection.Country == "Egypt") %>%
  mutate(H5 = case_when(grepl("H5", Subtype) ~ "positive",
                        is.na(Subtype) ~ NA,
                        Subtype %in% c("HPAI", "LPAI", "", "Mixed,mixed(Mixed", "Mixed", "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10") ~ NA,
                        grepl("Hx", Subtype) ~ NA,
                        grepl("HX", Subtype) ~ NA,
                        TRUE ~ "negative"),
         H9 = case_when(grepl("H9", Subtype) ~ "positive",
                        is.na(Subtype) ~ NA,
                        Subtype %in% c("HPAI", "LPAI", "", "Mixed,mixed(Mixed", "Mixed", "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10") ~ NA,
                        grepl("Hx", Subtype) ~ NA,
                        grepl("HX", Subtype) ~ NA,
                        TRUE ~ "negative")) %>%
  mutate(Collection.Date =  as.Date(Collection.Date, "%Y-%m-%d")) %>%
  select(all_of(c("Collection.Date", "Collection.Latitude", "Collection.Longitude", "Pathogen.Test.Result",
                  "Host.Species", "H5", "H9"))) %>%
  mutate(Source = "BVBRC") %>%
  mutate(Pathogen.Test.Result = tolower(Pathogen.Test.Result)) %>%
  rename(Latitude = Collection.Latitude,
         Longitude = Collection.Longitude,
         Date = Collection.Date,
         Species  = Host.Species,
         Result = Pathogen.Test.Result)

# Combine data
egy_data <- bind_rows(mn_data, fao_data, wahis_data, bvbrc_data) %>% as.data.frame

# check for any NA values
egy_data %>% summarise_all(~ sum(is.na(.)))

# check spp
egy_data %>% pull(Species) %>% unique %>% sort

# Check subtypes
table(egy_data$H5, useNA = "always")
table(egy_data$H9, useNA = "always")