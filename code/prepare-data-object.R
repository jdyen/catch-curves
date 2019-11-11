# prepare all data
setwd("~/Dropbox/research/catch-curves/")

# need some packages
library(lubridate)
library(dplyr)

# load some helper functions
source("code/data-helpers.R")
source("code/helpers.R")

# load the data
source("code/load-and-check-corrected-survey-data.R")

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
saveRDS(data_matrix, file = "data/data-loaded-Nov19.rds")

# load the flow data
source("code/load-flow-data-by-site.R")

# save flow data
saveRDS(flow_data, file = "data/flow-data-loaded-Nov19.rds")
