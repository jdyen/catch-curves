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
n_age <- 4
dens_type <- greta.integrated:::beverton()
size_breaks <- c(0, 50, 100, 200, 500, 1000, 2000, 10000, 60000)
n_size <- length(size_breaks) - 1

# mcmc settings
warmup <- 1000
n_samples <- 1000
chains <- 2

# filter survey and flow data to Murray cod
flow_data <- flow_data[survey_data$Common.Name %in% c('Murray Cod', 'Murray cod'), ]
survey_data <- survey_data[survey_data$Common.Name %in% c('murraycod', 'Murray Cod', 'Murray cod'), ]
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

## TEMPORARY: remove OVENS data due to missing weight data
flow_out[[6]] <- NULL
size_out[[6]] <- NULL

# compile size info
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
                                          density = dens_type)
}

# connect sizes to ages
alpha_set <- ifelse(is.na(size_age_obs), 0, 1)
age_size_dist <- dirichlet(alpha = alpha_set)

# estimate size-age distribution from data
size_sums <- apply(size_age_obs, 1, sum)

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

# number of years alive
n_alive <- final_obs - first_obs

# are any individuals never recaptured?
single_obs <- first_obs == final_obs

# focus on those with >1 observation
n_alive <- n_alive[!single_obs]
first_sizes <- first_size_class[!single_obs]
first_seen <- first_obs[!single_obs]
last_seen <- final_obs[!single_obs]

# probabiity of age j given size class k (should be n_obs x n_age)
p_age_first_capture <- age_size_dist[first_sizes, ]

# create Phi = P(surv_n_alive | age_at_time_1) (n_alive x n_age)
max_alive <- max(n_alive)
survival_lmr <- mat$surv_params[system_list == 'LOWERMURRAY', ]
survival_lmr <- cbind(survival_lmr,
                      do.call(cbind,
                              ## INDEXING ERRORS MIGHT BE DUE TO THIS
                              lapply(seq_len(max_alive),
                                     function(x) survival_lmr[, n_age])))

# pre-calculate all possible survival trajectories
first_seen_lifespan <- matrix(as.numeric(unlist(strsplit(unique(paste(first_seen, n_alive, sep = '_')), '_'))), ncol = 2, byrow = TRUE)
id_match <- match(paste(first_seen, n_alive, sep = '_'), unique(paste(first_seen, n_alive, sep = '_')))
idx <- apply(first_seen_lifespan, 1, function(x) x[1]:(x[1] + x[2]))
surv_mat <- ones(nrow(first_seen_lifespan), n_age)
for (i in seq_len(nrow(first_seen_lifespan))) {
  mat_tmp <- survival_lmr[idx[[i]], ]
  out <- NULL
  for (j in seq_len(n_age)) {
    out <- c(out, seq(1, length(idx[[i]]) ^ 2, by = (length(idx[[i]]) + 1)) + (j - 1) * length(idx[[i]]))
  }
  idy <- rep(seq_len(n_age), each = length(idx[[i]]))
  surv_mat[i, ] <- tapply(mat_tmp[out], idy, 'prod')
}

# what are survival probabilities for each age at first sight?
p_survival_hist_age <- surv_mat[id_match, ]

# calculate P(survival_history | size = k)
p_survival_hist_size <- rowSums(p_survival_hist_age * p_age_first_capture)

# probs of binary capture history, assume detection is constant
detection <- beta(1, 1)
binary_capture_history <- apply(catch_size_class, 1, function(x) x[min(which(x > 0)):max(which(x > 0))])
binary_capture_history <- ifelse(do.call(c, binary_capture_history) > 0, 1, 0)

# probs of final detection (for each age??)
# p(never_observed_again | final_size, final_obs) = \Sum_ages p(final_age | final_size) p(not_observed | final_obs, final_age) 
# p_final_obs <- SOMETHING
# not_detected is independent of age/size
# survival depends on year and age
# set up recursively (but still needs one entry for each individual)

# parameters
surv_total <- mat$survival
survival <- mat$surv_params
recruitment <- mat$fec_params
growth <- mat$growth_params

# convert sizes to ages
modelled_sizes <- vector('list', length = length(age_dist))
for (i in seq_along(age_dist))
  modelled_sizes[[i]] <- age_size_dist %*% age_dist[[i]]

# create vectors of fitted and observed data
mu <- do.call(c, modelled_sizes)
size_vec <- do.call(c, size_binned)

# extract parameters
mats <- mat$matrix
surv_coef <- mat$coefs_surv
fec_coef <- mat$coefs_fec
init_abund <- mat$init
gamma_year_surv <- mat$gamma_year_surv
gamma_year_fec <- mat$gamma_year_fec
dens_param <- mat$dens_param

# size obs model
distribution(size_vec) <- poisson(mu)

# size-age model
distribution(size_age_obs) <- multinomial(size = size_sums, p = age_size_dist)

# cmr model
# distribution(binary_capture_history) <- binomial(size = 1, p = detection)
# survival_hist_ones <- ones(length(p_survival_hist_size))
# distribution(survival_hist_ones) <- binomial(size = 1, p = p_survival_hist_size)
##
# final_obs_ones <- ones(length(p_final_obs))
# distribution(final_obs_ones) <- binomial(size = 1, p = p_final_obs)

# create greta model
mod <- model(survival, recruitment, dens_param,
             age_size_dist, mu)

# set initial values
inits <- initials(surv_coef = matrix(0.0, nrow(surv_coef), ncol(surv_coef)),
                  fec_coef = matrix(0.0, nrow(fec_coef), ncol(fec_coef)),
                  init_abund = matrix(1.0, n_age, n_site),
                  dens_param = rep(1e-6, n_site))

# sample from model
samples <- mcmc(mod,
                sampler = hmc(Lmin = 10, Lmax = 200, epsilon = 0.1),
                n_samples = n_samples, warmup = warmup,
                thin = max(1, floor(n_samples / 5000)),
                chains = chains,
                initial_values = inits)

# summarise fitted modelthe
mod_summary <- summary(samples)

# pull out fitted vals
fitted_vals <- mod_summary$quantiles[grep('mu\\[', rownames(mod_summary$quantiles)), ]

# pull out age-size distribution
age_size_est <- mod_summary$quantiles[grep("age_size_dist\\[", rownames(mod_summary$quantiles)), ]

# reconstruct Leslie matrices
recruitment_est <- mod_summary$quantiles[grep('^recruitment\\[', rownames(mod_summary$quantiles)), ]
survival_est <- mod_summary$quantiles[grep('^survival\\[', rownames(mod_summary$quantiles)), ]
recruitment_median <- recruitment_est[, "50%"]
survival_median <- matrix(survival_est[, "50%"], ncol = n_age)

# compile fitted matrices
mats_est <- list()
for (i in seq_len(nrow(survival_median))) {
  mat_tmp <- matrix(0, n_age, n_age)
  mat_tmp[n_age, n_age] <- survival_median[i, n_age]
  lower_diag <- row(mat_tmp) - col(mat_tmp) == 1
  mat_tmp[lower_diag] <- survival_median[i, seq_len(n_age - 1)]
  mat_tmp[1, n_age] <- recruitment_median[i]
  mats_est[[i]] <- mat_tmp
}
