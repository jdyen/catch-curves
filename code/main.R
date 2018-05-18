# set file paths
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(lubridate)
library(greta)

# source helper functions
source("./code/length_to_mass_calculations.R")
source("./code/length_to_mass_calculations_sra.R")
source("./code/length_age_conversion.R")
source("./code/helpers.R")

# load data
source("./code/load-flow-data.R")
source("./code/load-data.R")

# check corrs in lagged flow data
# mc_catch_curve$flow_ind <- remove_correlated(mc_catch_curve$flow)

# OR: use PCA to calculate a few flow components
## NOTE: this gives different PCs to each species -- very similar but
##       might affect interpretability.
##       Alternative is to pre-calc PCs and then assign (but this won't
##       preserve uncorrelated variables)
mc_pc <- calc_flow_pc(mc_catch_curve$flow, scale = FALSE)
mc_catch_curve$flow_pc <- mc_pc$scores[, seq_len(3)]
tc_pc <- calc_flow_pc(tc_catch_curve$flow, scale = FALSE)
tc_catch_curve$flow_pc <- tc_pc$scores[, seq_len(3)]
gp_pc <- calc_flow_pc(gp_catch_curve$flow, scale = FALSE)
gp_catch_curve$flow_pc <- gp_pc$scores[, seq_len(3)]
sp_pc <- calc_flow_pc(sp_catch_curve$flow, scale = FALSE)
sp_catch_curve$flow_pc <- sp_pc$scores[, seq_len(3)]

# set up mvn flow model (only have MC, TC, GP, SP) {could use size classes and add rainbows}
# age x species correlation matrix?
#  age_corr %x% spp_corr
num_ages <- 6
all_data <- list(mc_catch_curve,
                 tc_catch_curve,
                 gp_catch_curve,
                 sp_catch_curve)
age_data <- do.call("rbind", sapply(all_data, function(x) x$age_dist[, seq_len(num_ages)]))
flow_data <- do.call("rbind", sapply(all_data, function(x) x$flow_pc))
flow_data <- scale(flow_data)
info_data <- do.call("rbind", sapply(all_data, function(x) as.matrix(x$info)))
info_data <- as.data.frame(info_data)
info_data$spp <- rep(c("MC", "TC", "GP", "SP"),
                     times = sapply(all_data, function(x) nrow(x$age_dist)))

# create a separable covariance matrix with spp and age components
covar_mat <- create_greta_covar(n_sp = length(unique(info_data$spp)),
                                n_age = ncol(age_data))

