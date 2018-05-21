# greta scripts for catch curve analysis

create_greta_covar <- function(n_sp, n_age,
                               eta_sp = 1, eta_age = 4,
                               sigma_mean = 0.0, sigma_sd = 1.0) {
  
  if (length(sigma_mean) == 1)
    sigma_mean <- rep(sigma_mean, 2)
  if (length(sigma_sd) == 1)
    sigma_sd <- rep(sigma_sd, 2)
  
  species_corr_est <- lkj_correlation(eta = eta_sp, dim = n_sp)
  age_corr_est <- lkj_correlation(eta = eta_age, dim = n_age)
  species_sigma <- lognormal(mean = sigma_mean[1], sd = sigma_sd[1], dim = n_sp)
  age_sigma <- lognormal(mean = sigma_mean[2], sd = sigma_sd[2], dim = n_age)
  species_covar_est <- (species_sigma %*% t(species_sigma)) * species_corr_est
  age_covar_est <- age_sigma %*% t(age_sigma) * age_corr_est
  covar_est <- kronecker(species_covar_est, age_covar_est)
  
  covar_est
  
}
# check corrs in lagged flow data
# mc_catch_curve$flow_ind <- remove_correlated(mc_catch_curve$flow)

# OR: use PCA to calculate a few flow components
## NOTE: this gives different PCs to each species -- very similar but
##       might affect interpretability.
##       Alternative is to pre-calc PCs and then assign (but this won't
##       preserve uncorrelated variables)
## Could define flow_pca on the rbind() data set defined below (this is probably correct)
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
n_age <- 6
all_data <- list(mc_catch_curve,
                 tc_catch_curve,
                 gp_catch_curve,
                 sp_catch_curve)
age_data <- do.call("rbind", sapply(all_data, function(x) x$age_dist[, seq_len(n_age)]))
flow_data <- do.call("rbind", sapply(all_data, function(x) x$flow_pc))
flow_data <- scale(flow_data)
info_data <- do.call("rbind", sapply(all_data, function(x) as.matrix(x$info)))
info_data <- as.data.frame(info_data)
info_data$spp <- rep(c("MC", "TC", "GP", "SP"),
                     times = sapply(all_data, function(x) nrow(x$age_dist)))

# extract counters
n_sp <- length(unique(info_data$spp))
n_obs <- nrow(age_data)

# create a separable covariance matrix with spp and age components
covar_mat <- create_greta_covar(n_sp = n_sp,
                                n_age = n_age)

# set priors for covariate model
alpha_mean <- normal(mean = 0.0, sd = 1.0, dim = 1)
beta_mean <- normal(mean = 0.0, sd = 1.0, dim = 1)
alpha <- normal(mean = alpha_mean, sd = sigma_alpha, dim = n_sp)
beta <- normal(mean = beta_mean, sd = sigma_beta, dim = n_sp)
gamma <- normal(mean = 0.0, sd = sigma_gamma, dim = n_sp)
delta <- normal(mean = 0.0, sd = sigma_delta, dim = n_age)

# define response
y <- c(age_data)
## DON'T DO THIS
# set up a matrix of nsite x (nsp * nage)
# This will be the observations -- need rows to be unique sites and pad with zeros
#   for unobserved species (check this is OK -- they might not have been sampling
#   for all spp in all surveys)

# define indices
sp <- as.integer(as.factor(info_data$spp))
age <- rep(seq_len(n_age), each = n_obs)

# set linear predictor
# need to include year, reach (this must be river*reach, i.e., reach within river),
# river as random effects
# need a size * age random effect (I think)
# probably need a year * reach random effect (or year * river)
#
# there will be a nsp * n_age vector of means
# this will need to be swept through a matrix that is the
# outer product of sp*age regression coefficients with 
# the flow data

# THIS IS WRONG
#mu <- alpha[sp] + flow_data %*% beta[sp] + gamma[sp] + delta[age]

# set likelihoods
# see if can setup a Stan-style normal() call with 
#   work done on the covariances outside of this -- will be
#   much faster and easier.
# otherwise will need to loop through multivariate_normal calls (although
#   check how nick is doing this with JSDMs)