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

# save loaded data
saveRDS(data_matrix, file = "data/data-loaded-Aug19.rds")

# load the flow data
source("code/load-flow-data-by-site.R")

# save flow data
saveRDS(flow_data, file = "data/flow-data-loaded-Aug19.rds")
