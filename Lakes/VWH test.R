library(ncdf4) # package for netcdf manipulation
library(ggplot2) # package for plotting
library(ggpubr)
library(lubridate)
library(dplyr)
library(CFtime)
library(zoo)
library(plot3D)
library(tidyr)
library(rLakeAnalyzer)

setwd('C:/Users/nkwalale/Documents/Download_Test/Lakes')

# list files
file_list <- list.files()

# Read the hypso data frames
rappbode <- read.csv("hypsograph.csv", header = TRUE)
bosumtwi <- read.csv("Bosumtwi_hypsometry.csv", header = TRUE)
kivu <- read.csv("Kivu_hypsometry.csv", header = TRUE)

# List all files
file_list <- list.files()
strat_files <- grep("_watertemp", file_list, value = TRUE)
wattemp_files <- grep("_watertemp", file_list, value = TRUE)

# Extract the model names, gcm, scenario, and lake names
model <- sapply(strsplit(wattemp_files, split = "_", fixed = TRUE), `[`, 1)
gcm <- sapply(strsplit(wattemp_files, split = "_", fixed = TRUE), `[`, 2)
scen <- sapply(strsplit(wattemp_files, split = "_", fixed = TRUE), `[`, 4)
lake <- sapply(strsplit(wattemp_files, split = "_", fixed = TRUE), `[`, 8)

# Create an empty data frame to store the results
ISIMIP_data <- data.frame()

# Unique combinations of gcm, scen, and lake
gcm_name <- unique(gcm)
scen_name <- unique(scen)
lake_name <- unique(lake)
model_name <- unique(model)

# Define the function to select the correct hypso data frame based on the lake name
get_hypso <- function(lake) {
  switch(lake,
         "bosumtwi" = bosumtwi,
         "kivu" = kivu,
         "rappbode" = rappbode,
         stop("Unknown lake"))
}

# Nested for loops for each combination of gcm, scen, and lake
for (g in gcm_name) {
  for (s in scen_name) {
    for (l in lake_name) {
      for (m in model_name) {
      
      # Filter files for the specific gcm, scen, lake, and model
      wattemp_files_filtered <- wattemp_files[gcm == g & scen == s & lake == l & model == m]
      
      # Select the appropriate hypso data frame based on the lake name
      hypso <- get_hypso(l)
      
      # Loop through each filtered file set
      for (i in seq_along(wattemp_files_filtered)) {
        
        # Open and read watertemp file
        nc_data <- nc_open(wattemp_files_filtered[i])
        
        t <- ncvar_get(nc_data, "time")
        cf <- CFtime(nc_data$dim$time$units, nc_data$dim$time$calendar, nc_data$dim$time$vals)
        dates <- as_timestamp(cf)
        dates <- as.POSIXct(dates, tz = "UTC")
        
        depth <- ncvar_get(nc_data, "depth")
        watertemp <- t(ncvar_get(nc_data, "watertemp"))
        
        nc_close(nc_data)
        
        watertemp <- data.frame(watertemp)
        watertemp$date <- dates
        colnames(watertemp) <- c(as.character(depth), "date")
        
        if (grepl('flake', wattemp_files_filtered[i], fixed = TRUE)) {
          depth <- c(depth, seq(max(depth), max(hypso$Depth_meter), 2)) # get same structure as other models
          idx <- which(!as.character(depth) %in% colnames(watertemp))
          
          add_m <- matrix(rep(watertemp[, ncol(watertemp)-1], length(idx)), ncol = length(idx))
          df_m <- as.data.frame(add_m)
          colnames(df_m) <- as.character(depth[idx])
          
          watertemp <- cbind(watertemp, df_m)
        }
        
        # Create a data frame for the current set of data
        data_nc <- pivot_longer(watertemp, cols = -date,
                                names_to = "depth",
                                values_to = "watertemp",
                                names_transform = as.numeric) |>
          mutate(scen = s, model = m, gcm = g, lake = l)
        
        hypo_temp <- data_nc |> mutate(watertemp = watertemp - 273.15) |>
          group_by(date, gcm, scen, model, lake) |>
          summarize(hyp_temp = hypo.temperature(watertemp, depth, hypso$Area_meterSquared, hypso$Depth_meter), .groups = 'drop') |>
          mutate(year = year(date)) |>
          group_by(gcm, scen, model, lake, year) |>
          summarize(av_hyp_temp = mean(hyp_temp, na.rm = TRUE), .groups = 'drop')
        
        ISIMIP_data <- rbind(ISIMIP_data, as.data.frame(hypo_temp))
      }
    }
  }
}
}

write.table(ISIMIP_data, file = "VWH_temps.txt", sep = "\t", row.names = FALSE, col.names = TRUE)