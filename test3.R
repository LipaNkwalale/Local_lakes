file_list <- list.files()

b_temp_files <- grep("bottemp", file_list, value = TRUE)
s_temp_files <- grep("surftemp", file_list, value = TRUE)
strat_files <- grep("_strat", file_list, value = TRUE)

model <- sapply(strsplit(b_temp_files, split = "_", fixed = TRUE), `[`, 1) # extracting the model names
gcm <- sapply(strsplit(b_temp_files, split = "_", fixed = TRUE), `[`, 2)
scen <- sapply(strsplit(b_temp_files, split = "_", fixed = TRUE), `[`, 4)
name <- sapply(strsplit(b_temp_files, split = "_", fixed = TRUE), `[`, 8)

working_data <- data.frame()

for (i in seq_along(model)) {
  for (j in seq_along(gcm)) {
    
  }
  for (j in seq_along(s_temp_files)) {
    if (gsub("bottemp", "surftemp", b_temp_files[i]) == s_temp_files[j]) {
      for (k in seq_along(strat_files)) {
        if (gsub("bottemp", "strat", b_temp_files[i]) == strat_files[k]) {
          
          nc_data <- nc_open(b_temp_files[i]) # Open bottemp file
          dates <- as.Date(as_timestamp(CFtime(nc_data$dim$time$units, nc_data$dim$time$calendar, nc_data$dim$time$vals)))
          bt <- ncvar_get(nc_data, "bottemp")
          nc_close(nc_data)
          
          nc_data <- nc_open(s_temp_files[j]) # Open surftemp file
          st <- ncvar_get(nc_data, "surftemp")
          nc_close(nc_data)
          
          nc_data <- nc_open(strat_files[k]) # Open strat file
          strat <- ncvar_get(nc_data, "strat")
          nc_close(nc_data)
          
          data_nc <- data.frame(
            dates = dates,
            bt = bt,
            st = st,
            strat = strat,
            model = model[i],
            gcm = gcm[i], 
            scen = scen[i],
            name = name[i]
          )
          
          working_data <- rbind(working_data, data_nc)
          working_data$year <- year(working_data$dates)
          working_data$dates <- as.Date(dates)
        }
      }
    }
  }
}



#LER Stratification#####
ISIMIP_LER <- working_data|>
  filter(st > bt) #removing winter stratification
ISIMIP_LER$strat <- ifelse(ISIMIP_LER$strat == 0, NA, ISIMIP_LER$strat)

# Find longest contiguous non-NA values for each year
longest_stretch_LER <- ISIMIP_LER |>
  mutate(year = year(dates)) |>
  group_by(model, gcm, name, year) |>
  group_modify(~ {
    ts_zoo <- zoo(.x$strat, .x$dates)
    longest_zoo <- na.contiguous(ts_zoo)
    longest_dates <- index(longest_zoo)
    .x[.x$dates %in% longest_dates, ]
  }) |>
  ungroup()

longer_strat_LER <- data.frame()
longer_strat_LER <- longest_stretch_LER |> drop_na(strat)

scenarios <- unique(longer_strat_LER$scen)
models <- unique(longer_strat_LER$model)

# Define the function to get yearly average
get_yearly_averages <- function(scen_name, model_name, data) {
  strat_data <- data |>
    filter(scen == scen_name) |>
    filter(model == model_name) |>
    group_by(year) |>
    summarize(strat.dur = n(), ave.temp = mean(bt, na.rm = TRUE)) |>
    mutate(Scen = scen_name, Model = model_name)
  
  data.frame(
    Year = strat_data$year,
    Strat.dur = strat_data$strat.dur,
    Temp_K = strat_data$ave.temp,
    Model = strat_data$Model,
    Scenario = strat_data$Scen
  )
}


combinations <- expand.grid(scen_name = scenarios, model_name = models, stringsAsFactors = FALSE)
Bosumtwi_LER <- pmap_dfr(combinations, get_yearly_averages, data = longer_strat_LER)
