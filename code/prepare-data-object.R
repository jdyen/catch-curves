# prepare all data
setwd("~/Dropbox/research/catch-curves/")

# need some packages
library(lubridate)
library(dplyr)

# load some helper functions
source("code/length_to_mass_calculations.R")
source("code/length_to_mass_calculations_sra.R")
source("code/data-helpers.R")
source("code/helpers.R")

# load the data
source("code/load-survey-data.R")

# save loaded data
saveRDS(alldat, file = "data/data-loaded-Jul19.rds")

# load the flow data
source("code/load-flow-data-by-site.R")

# save flow data
saveRDS(flow_data, file = "data/flow-data-loaded-Jul19.rds")

# filter the data to major systems
# source("code/filter-survey-data.R")


