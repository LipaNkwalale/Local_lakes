library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(ggpubr)
library(lubridate)
library(data.table)
library(dplyr)
library(CFtime)
library(zoo)
library(tidyverse)
library(purrr)

setwd("C:/Users/nkwalale/Documents/Download_Test")


# List all files
file_list <- list.files()

# Filter files based on specific patterns
b_temp_files <- grep("bottemp", file_list, value = TRUE)
s_temp_files <- grep("surftemp", file_list, value = TRUE)
strat_files <- grep("_strat", file_list, value = TRUE)

# Extract the model names, gcm, scenario, and lake names
model <- sapply(strsplit(b_temp_files, split = "_", fixed = TRUE), `[`, 1)
gcm <- sapply(strsplit(b_temp_files, split = "_", fixed = TRUE), `[`, 2)
scen <- sapply(strsplit(b_temp_files, split = "_", fixed = TRUE), `[`, 4)
lake <- sapply(strsplit(b_temp_files, split = "_", fixed = TRUE), `[`, 8)

# Create an empty data frame to store the results
ISIMIP_data <- data.frame()

# Unique combinations of gcm, scen, and lake
unique_gcm <- unique(gcm)
unique_scen <- unique(scen)
unique_lake <- unique(lake)

# Nested for loops for each combination of gcm, scen, and lake
for (g in unique_gcm) {
  for (s in unique_scen) {
    for (l in unique_lake) {
      
      # Filter files for the specific gcm, scen, and lake
      b_temp_files_filtered <- b_temp_files[gcm == g & scen == s & lake == l]
      s_temp_files_filtered <- s_temp_files[gcm == g & scen == s & lake == l]
      strat_files_filtered <- strat_files[gcm == g & scen == s & lake == l]
      
      # Loop through each filtered file set
      for (i in seq_along(b_temp_files_filtered)) {
        
        # Open and read bottemp file
        nc_data <- nc_open(b_temp_files_filtered[i])
        dates <- as.Date(as_timestamp(CFtime(nc_data$dim$time$units, nc_data$dim$time$calendar, nc_data$dim$time$vals)))
        bt <- ncvar_get(nc_data, "bottemp")
        nc_close(nc_data)
        
        # Open and read surftemp file
        nc_data <- nc_open(s_temp_files_filtered[i])
        st <- ncvar_get(nc_data, "surftemp")
        nc_close(nc_data)
        
        # Open and read strat file
        nc_data <- nc_open(strat_files_filtered[i])
        strat <- ncvar_get(nc_data, "strat")
        nc_close(nc_data)
        
        # Create a data frame for the current set of data
        data_nc <- data.frame(
          dates = dates,
          bt = bt,
          st = st,
          strat = strat,
          model = model[gcm == g & scen == s & lake == l][i],
          gcm = g,
          scen = s,
          lake = l
        )
        
        # Append the data to the main data frame
        ISIMIP_data <- rbind(ISIMIP_data, data_nc)
        
      }
    }
  }
}

ISIMIP_data$year <- year(ISIMIP_data$dates)
ISIMIP_data$dates <- as.Date(dates)

#LER Stratification#####
ISIMIP <- ISIMIP_data |>
  filter(st > bt) #removing winter stratification
ISIMIP$strat <- ifelse(ISIMIP$strat == 0, NA, ISIMIP$strat)

longest_stretch_ISIMIP <- ISIMIP |>
  mutate(year = year(dates)) |>
  group_by(lake, gcm, scen, lake, year) |>
  group_modify(~ {
    ts_zoo <- zoo(.x$strat, .x$dates)
    
    if (all(is.na(ts_zoo))) {
      return(data.frame()) # Return an empty data frame if all values are NA
    }
    
    longest_zoo <- na.contiguous(ts_zoo)
    longest_dates <- index(longest_zoo)
    
    .x %>% filter(dates %in% longest_dates)
  }) %>%
  ungroup()


longest_stretch_ISIMIP <- ISIMIP %>%
  mutate(year = year(dates)) %>%
  group_by(lake, gcm, scen, year) %>%
  group_modify(~ {
    ts_zoo <- zoo(.x$strat, .x$dates)
    
    if (all(is.na(ts_zoo))) {
      return(tibble()) # Return an empty tibble if all values are NA
    }
    
    longest_zoo <- na.contiguous(ts_zoo)
    longest_dates <- index(longest_zoo)
    
    .x %>% filter(dates %in% longest_dates)
  }) %>%
  ungroup()

longer_strat_LER <- data.frame()
longer_strat_LER <- longest_stretch_ISIMIP |> drop_na(strat)

scenarios <- unique(longer_strat_LER$scen)
models <- unique(longer_strat_LER$model)
gcm <- unique(longer_strat_LER$gcm)
lake <- unique(longer_strat_LER$lake)

# Define the function to get yearly average
get_yearly_averages <- function(gcm_name, scen_name, model_name, lake_name, data) {
  strat_data <- data |>
    filter(gcm == gcm_name) |>
    filter(scen == scen_name) |>
    filter(model == model_name) |>
    filter(lake == lake_name) |>
    group_by(year) |>
    summarize(strat.dur = n(), ave.temp = mean(bt, na.rm = TRUE), .groups = 'drop') |>
    mutate(GCM = gcm_name, Scen = scen_name, Model = model_name, lake = lake_name)
  
  if (nrow(strat_data) == 0) {
    return(data.frame(
      GCM = character(0),
      Year = integer(0),
      Strat.dur = integer(0),
      Temp_K = numeric(0),
      Model = character(0),
      Scenario = character(0),
      Lake = character(0)
    ))
  }
  
  data.frame(
    GCM = strat_data$GCM,
    Year = strat_data$year,
    Strat.dur = strat_data$strat.dur,
    Temp_K = strat_data$ave.temp,
    Model = strat_data$Model,
    Scenario = strat_data$Scen,
    Lake = strat_data$lake
  )
}

combinations <- expand.grid(gcm_name = gcm, scen_name = scenarios, model_name = models, lake_name = lake, stringsAsFactors = FALSE)
ISI_lakes <- pmap_dfr(combinations, ~get_yearly_averages(..1, ..2, ..3, ..4, longer_strat_LER))


write.table(ISI_lakes, file = "three_lakes.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
