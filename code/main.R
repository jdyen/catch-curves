# set file paths
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(lubridate)

# source helper functions
source("./code/length_to_mass_calculations.R")
source("./code/length_to_mass_calculations_sra.R")
source("./code/length_age_conversion.R")
source("./code/helpers.R")

# load data
source("./code/load-flow-data.R")
source("./code/load-data.R")

# check corrs in lagged flow data
# remove_correlated(mc_catch_curve$flow)

# OR: use PCA to calculate a few flow components
mc_pc <- calc_flow_pc(mc_catch_curve$flow, scale = FALSE)
mc_catch_curve$flow_pc <- mc_pc$scores[, 1:3]

# set up mvn flow model (only have MC, TC, GP, SP) {could use size classes and add rainbows}
# size x species correlation matrix?
#  size_corr %x% spp_corr


