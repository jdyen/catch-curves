# clear workspace
rm(list = ls())

# set working directory
setwd('~/Dropbox/research/catch-curves/')

# load libraries
library(greta.dynamics)

# load helper functions
source('./code/helpers.R')

# load data
cmr_data <- read.csv('data/compiled-mc-cmr-data.csv', row.names = 1, stringsAsFactors = FALSE)
survey_data <- read.csv('data/compiled-survey-data.csv', row.names = 1, stringsAsFactors = FALSE)
size_age_data <- read.csv('data/compiled-size-age-data.csv', row.names = 1, stringsAsFactors = FALSE)
flow_data <- read.csv('data/compiled-flow-predictors.csv', row.names = 1, stringsAsFactors = FALSE)

# model settings
n_age <- 10
dens_type <- 'bh'
size_breaks <- c(0, 50, 100, 200, 500, 1000, 2000, 5000, 60000)
n_size <- length(size_breaks) - 1

# mcmc settings
warmup <- 5000
n_samples <- 5000
chains <- 4

# filter survey and flow data to Murray cod
flow_data <- flow_data[survey_data$Common.Name %in% c('Murray Cod', 'Murray cod'), ]
survey_data <- survey_data[survey_data$Common.Name %in% c('Murray Cod', 'Murray cod'), ]
flow_data <- flow_data[!is.na(survey_data$Total.Sampled), ]
survey_data <- survey_data[!is.na(survey_data$Total.Sampled), ]
flow_data <- flow_data[survey_data$Total.Sampled > 0, ]
survey_data <- survey_data[survey_data$Total.Sampled > 0, ]

# prepare survey data
size_data <- tapply(survey_data$WEIGHT, list(survey_data$YEAR, survey_data$SYSTEM), unique)
flow_data <- list(mannf_mld = tapply(flow_data$mannf_mld, list(survey_data$YEAR, survey_data$SYSTEM), function(x) mean(unique(x))),
                  covaf_mld = tapply(flow_data$covaf_mld, list(survey_data$YEAR, survey_data$SYSTEM), function(x) mean(unique(x))),
                  madpth_m = tapply(flow_data$madpth_m, list(survey_data$YEAR, survey_data$SYSTEM), function(x) mean(unique(x))),
                  cvdpth_m = tapply(flow_data$cvdpth_m, list(survey_data$YEAR, survey_data$SYSTEM), function(x) mean(unique(x))))
size_out <- vector('list', length = ncol(size_data))
flow_out <- vector('list', length = ncol(size_data))
for (i in seq_len(ncol(size_data))) {
  size_tmp <- size_data[, i]
  flow_tmp <- do.call(cbind, lapply(flow_data, function(x) x[, i]))
  flow_tmp <- flow_tmp[min(which(!sapply(size_tmp, is.null))):max(which(!sapply(size_tmp, is.null))), ]
  size_tmp <- size_tmp[min(which(!sapply(size_tmp, is.null))):max(which(!sapply(size_tmp, is.null)))]
  size_out[[i]] <- size_tmp
  flow_out[[i]] <- flow_tmp
}
names(size_out) <- names(flow_out) <- colnames(size_data)
flow_out[[7]] <- NULL
size_out[[7]] <- NULL

# compile size info
### IDEALLY flatten this but is impossible(ish)
size_binned <- vector('list', length = length(size_out))
for (i in seq_along(size_out)) {
  size_binned[[i]] <- sapply(size_out[[i]], hist_fn, breaks = size_breaks)
}
size_mat <- t(do.call(cbind, size_binned))

# summarise size-age data
size_age_data$size_class <- cut(size_age_data$weight_g, size_breaks, labels = FALSE)
size_age_data$ones <- rep(1, nrow(size_age_data))
size_age_data$age_trunc <- ifelse(size_age_data$age_round >= n_age, n_age - 1, size_age_data$age_round)
size_age_obs <- matrix(0, nrow = n_size, ncol = n_age)
for (i in seq_len(n_size)) {
  size_age_sub <- size_age_data[size_age_data$size_class == i, ]
  size_age_tmp <- tapply(size_age_sub$ones, size_age_sub$age_trunc, sum)
  size_age_obs[i, match(names(size_age_tmp), seq_len(n_age) - 1)] <- size_age_tmp
}

# convert flow data to matrix
flow_mat <- do.call(rbind, flow_out)
flow_mat[apply(flow_mat, 1, function(x) all(is.na(x))), ] <- rep(apply(flow_mat, 2, mean, na.rm = TRUE), each = sum(is.na(flow_mat[, 1])))
flow_mat <- scale(flow_mat)

# setup model settings
n_time <- sapply(size_out, length)
system_list <- rep(names(size_out), times = n_time)
all_systems <- unique(system_list)
n_site <- length(all_systems)

# initialise matrix model
mat <- leslie_matrix(n_age = n_age,
                     density_dependence = dens_type,
                     predictors = flow_mat,
                     params = list(fec_stages = n_age,
                                   n_site = n_site,
                                   surv_sd = 1.0,
                                   fec_sd = 1.0))

# set up model of dynamics
age_dist <- vector('list', length = n_site)
for (i in seq_len(n_site)) {
  age_dist[[i]] <- iterate_matrix_dynamic(matrix = mat$matrix[system_list == all_systems[i], , ],
                                          initial_state = c(mat$init[, i]),
                                          dens_param = mat$dens_param[i],
                                          dens_form = dens_type)
}

# connect sizes to ages
alpha_set <- ifelse(is.na(size_age_obs), 0, 1)
age_size_dist <- dirichlet(alpha = alpha_set)

# estimate size-age distribution from data
size_sums <- apply(size_age_obs, 1, sum)
distribution(size_age_obs) <- multinomial(size = size_sums, p = age_size_dist)

# convert sizes to ages
modelled_sizes <- vector('list', length = length(age_dist))
for (i in seq_along(age_dist))
  modelled_sizes[[i]] <- age_size_dist %*% age_dist[[i]]

# mark-recapture model to estimate detection and survival probabilities
# calculate size-based catch history for each individual
catch_size <- with(cmr_data, tapply(weight, list(idfish, year), mean))
catch_size <- ifelse(is.na(catch_size), 0, catch_size)

# observed at least once?
observed <- apply(catch_size, 1, sum) > 0
catch_size <- catch_size[observed, ]

# convert to size classes
catch_size_class <- matrix(cut(catch_size, size_breaks, labels = FALSE),
                           ncol = ncol(catch_size))
catch_size_class <- ifelse(is.na(catch_size_class), 0, catch_size_class)

# first and final size classes
first_size_class <- apply(catch_size_class, 1, function(x) x[min(which(x > 0))])
final_size_class <- apply(catch_size_class, 1, function(x) x[max(which(x > 0))])

# has it shrunk *many* classes?
size_errors <- final_size_class < (first_size_class - 1)

# if so, remove these observations
catch_size_class <- catch_size_class[!size_errors, ]
first_size_class <- first_size_class[!size_errors]
final_size_class <- final_size_class[!size_errors]

# calculate first and final observations
first_obs <- apply(catch_size_class, 1, function(x) min(which(x > 0)))
final_obs <- apply(catch_size_class, 1, function(x) max(which(x > 0)))

# are any individuals never recaptured?
single_obs <- first_obs == final_obs

# parameters
survival <- mat$surv_params
recruitment <- mat$fec_params

# convert CMR sizes to ages (use age_size_dist)
size_first_capture <- as_data(matrix(hist(first_size_class,
                                          plot = FALSE,
                                          breaks = c(0:n_size + 0.5))$count,
                                     ncol = n_size))
age_first_capture <- size_first_capture %*% age_size_dist
p_age_first_capture <- ones(nrow(survival), ncol(survival))
for (i in seq_len(ncol(survival))) {
  p_age_first_capture[, i] <- tapply(c(survival[, seq_len(i)]),
                                     rep(seq_len(nrow(survival)), times = n_age),
                                     'prod')
}

nyears <- final_obs - first_obs
nyears <- nyears[!single_obs]
p_mid_obs <- p_age_first_capture[idx, age_init + nyears - 1]

# 
# calculate age at all captures
#
# calculate p(hist) = p(first_obs_age) * p(mid_obs_age) * p(final_obs_age)
#    - calculate probs of all trajectories and then sort with idx?
#

# create vectors of fitted and observed data
mu_vec <- do.call(c, modelled_sizes)
size_vec <- c(size_mat)

# extract parameters
mats <- mat$matrix
surv_coef <- mat$coefs_surv
fec_coef <- mat$coefs_fec
init_abund <- mat$init
gamma_year_surv <- mat$gamma_year_surv
gamma_year_fec <- mat$gamma_year_fec
dens_param <- mat$dens_param

# size obs model
distribution(size_vec) <- poisson(mu_vec)

# create greta model
mod <- model(surv_coef, fec_coef, age_size_dist, mu_vec, dens_param)

# set initial values
inits <- initials(surv_coef = matrix(0.0, nrow(surv_coef), ncol(surv_coef)),
                  fec_coef = matrix(0.0, nrow(fec_coef), ncol(fec_coef)),
                  gamma_year_surv = matrix(0.0, nrow(gamma_year_surv), ncol(gamma_year_surv)),
                  gamma_year_fec = matrix(0.0, nrow(gamma_year_fec), ncol(gamma_year_fec)),
                  init_abund = matrix(1.0, n_age, n_site),
                  dens_param = rep(1e-6, n_site))

# sample from model
samples <- mcmc(mod,
                n_samples = n_samples, warmup = warmup,
                thin = max(1, floor(n_samples / 5000)),
                chains = chains,
                initial_values = inits)

# summarise fitted model
mod_summary <- summary(samples)

# reconstruct Leslie matrices
median_fec_coef <- mod_summary$quantiles[grep('fec_coef', rownames(mod_summary$quantiles)), '50%']
median_surv_coef <- mod_summary$quantiles[grep('surv_coef', rownames(mod_summary$quantiles)), '50%']
median_surv_coef <- matrix(median_surv_coef, ncol = n_age)

# project/test fit etc.
survival_vals <- plogis(flow_mat %*% median_surv_coef)
fecundity_vals <- exp(flow_mat %*% median_fec_coef)
mats_est <- list()
for (i in seq_len(nrow(survival_vals))) {
  mat_tmp <- matrix(0, n_age, n_age)
  mat_tmp[n_age, n_age] <- survival_vals[i, n_age]
  lower_diag <- row(mat_tmp) - col(mat_tmp) == 1
  mat_tmp[lower_diag] <- survival_vals[i, seq_len(n_age - 1)]
  mat_tmp[1, n_age] <- fecundity_vals[i]
  mats_est[[i]] <- mat_tmp
}
