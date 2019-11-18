# load all data and save in pre-compiled data objects

# need some packages
library(lubridate)
library(dplyr)

# load some helper functions
source("code/helpers.R")

# load the data
source("code/load-survey-data.R")

# add some rows to data_matrix to get years prior
nyear <- 6
max_year <- max(year(data_matrix$date_formatted))
data_sub <- data_matrix[data_matrix$scientific_name == "Maccullochella peelii" & month(data_matrix$date_formatted) < 7, ]
min_years <- tapply(year(data_matrix$date_formatted), data_matrix$system, min, na.rm = TRUE)
idx <- match(unique(data_sub$system), data_sub$system)
year_seqs <- lapply(min_years, function(x) (x - nyear):max_year)
padded_years <- data_sub[rep(idx, times = sapply(year_seqs, length)), ]
padded_years <- padded_years[order(padded_years$system), ]
padded_years$date_formatted <- c(ymd(unlist(sapply(year_seqs, paste0, "-05-05"))))
padded_years$dataset <- rep("PADDED_TO_GET_FLOW_YEARS", nrow(padded_years))
data_matrix <- rbind(data_matrix, padded_years)

# save loaded data
saveRDS(data_matrix, file = "data/survey-data-loaded.rds")

# load the flow data
source("code/load-flow-data-by-year.R")

# save flow data
saveRDS(flow_data, file = "data/flow-data-loaded.rds")

# load otolith data for all species
all_oti <- read.csv("data/MurrayNtribs_otolith_data.csv", stringsAsFactors = FALSE)

# clean up species codes
all_oti$X.SPECIES[all_oti$X.SPECIES == "MC "] <- "MC"
all_oti$X.SPECIES[all_oti$X.SPECIES == "mc"] <- "MC"
all_oti$X.SPECIES[all_oti$X.SPECIES == "silver perch"] <- "Silver perch"
all_oti$X.SPECIES[all_oti$X.SPECIES == "Silver Perch"] <- "Silver perch"
all_oti$X.SPECIES[all_oti$X.SPECIES == "Sp"] <- "SP"
all_oti$X.SPECIES[all_oti$X.SPECIES == "Tc"] <- "TC"
all_oti <- all_oti[all_oti$X.SPECIES != "BT", ]
all_oti$SPECIES <- switch_names(all_oti$X.SPECIES)

# ages and lengths are character valued for some reason
all_oti$T_Length..mm. <- as.numeric(all_oti$T_Length..mm.)
all_oti$AGE <- as.numeric(all_oti$AGE)

# save otolith data object
saveRDS(all_oti, file = "data/otolith-data-loaded.rds")
