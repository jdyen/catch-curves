# load some helper functions
source("code/plot-helpers.R")

# load fitted models and data
mod <- readRDS("outputs/fitted/full-model.rds")
data_set <- readRDS("outputs/fitted/data-set.rds")
additional_data <- readRDS("outputs/fitted/additional-data.rds")

# unpack additional data
flow_scales <- additional_data$flow_scales
age_mat <- additional_data$age_mat
nyear <- additional_data$nyear
nsystem <- additional_data$nsystem
system_info <- additional_data$system_info
year_info <- additional_data$year_info
cohort_mat <- additional_data$cohort_mat

# plot some flow effects
system_names <- c("Broken", "Goulburn",
                  "King", "Murray", "Ovens")

i <- 4
plot_associations(mod, variable = "psprw_vec", data = data_set,
                  rescale = flow_scales, xlab = "Spring flow proportional to long-term winter average",
                  system = i, cohort = data_set$cohort_vec[data_set$system_vec == i][1],
                  age_set = 1)
