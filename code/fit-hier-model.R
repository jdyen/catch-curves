# set working directory
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(greta)

# load some helper functions
source("code/helpers.R")

# load compiled survey data
alldat <- readRDS("data/data-loaded-Jul19.rds")

# filter survey data to MC
to_keep <- alldat$Scientific.Name == "Maccullochella peelii"
alldat <- alldat[to_keep, ]
alldat$Common.Name <- NULL

# load otolith data
oti_data <- readRDS("data/otolith-data-loaded.rds")

# filter to MC
oti_data <- oti_data[oti_data$SPECIES == "Maccullochella peelii", ]

# need to load flow data
flow_data <- readRDS("data/flow-data-loaded-Jul19.rds")

# filter to MC
flow_data <- flow_data[to_keep, ]

# optional: subset to systems of interest
systems_to_keep <- c("BROKEN", "GOULBURN", "KING", "LOWERMURRAY", "OVENS")
flow_data <- flow_data[alldat$SYSTEM %in% systems_to_keep, ]
alldat <- alldat[alldat$SYSTEM %in% systems_to_keep, ]

# prepare survey data
survey_data <- data.frame(length_mm = alldat$totallength,
                          system = alldat$SYSTEM,
                          system_id = rebase_index(alldat$SYSTEM),
                          site = alldat$SITE_CODE,
                          site_id = rebase_index(alldat$SITE_CODE),
                          year = alldat$YEAR,
                          year_id = rebase_index(alldat$YEAR),
                          dataset = alldat$dataset,
                          dataset_id = rebase_index(alldat$dataset),
                          effort = alldat$total_no_passes * alldat$seconds)
flow_data <- flow_data[apply(survey_data, 1, function(x) !any(is.na(x))), ]
survey_data <- survey_data[apply(survey_data, 1, function(x) !any(is.na(x))), ]

# bin otolith data by age and size: offset by -0.4 (0-0.6 = YOY, 0.6-1.6 = 1YO, etc.)
oti_data$age_class <- cut(oti_data$AGE, breaks = c(-0.4:ceiling(max(oti_data$AGE, na.rm = TRUE))),
                          labels = FALSE)
size_breaks <- c(0, 150, 276, 469, 544, 610, 667, 718, 760, max(survey_data$length_mm, na.rm = TRUE))
oti_data$len_class <- cut(oti_data$T_Length..mm., breaks = size_breaks, labels = FALSE)
length_age_matrix <- classify(oti_data$len_class, oti_data$age_class)

# how many length classes do we need to keep to give all 0-4 year olds?
n_len <- max(which(length_age_matrix[, 5] > 0))

# how many age classes do we need to keep to give all 0:n_len length individuals?
n_age <- max(which(length_age_matrix[n_len, ] > 0))

# let's just keep those from the age_length_matrix
length_age_matrix <- length_age_matrix[seq_len(n_len), seq_len(n_age)]

# set otolith likelihood
# can set a 1:1 prior; could update to match a VB curve more closely
length_age_prior <- zeros(n_len, n_age) + 0.1
length_age_prior[row(length_age_prior) == col(length_age_prior)] <- 1
length_to_age <- dirichlet(alpha = length_age_prior)
distribution(length_age_matrix) <- multinomial(size = apply(length_age_matrix, 1, sum),
                                               prob = length_to_age)

# need to bin the survey data by lengths
response_matrix <- do.call(
  rbind, tapply(survey_data$length_mm,
                list(survey_data$system, survey_data$year),
                hist_fn, breaks = size_breaks))

# filter length class survey data to the length classes we care about (1:n_len)
response_matrix <- response_matrix[, seq_len(n_len)]

# define system/year predictors
system <- c(tapply(survey_data$system_id,
                   list(survey_data$system, survey_data$year),
                   unique))
system <- system[!is.na(system)]
year <- c(tapply(survey_data$year_id,
                   list(survey_data$system, survey_data$year),
                   unique))
year <- year[!is.na(year)]

# add a quadratic summer effect to the flow predictors
flow_data$prop_sum_win_sq <- flow_data$prop_sum_lt_win ^ 2
flow_data$prop_sum_win_sq_ym1 <- flow_data$prop_sum_lt_win_ym1 ^ 2

# compile flow predictors
vars_to_include <- c("rrang_spwn_mld", "rrang_spwn_mld_ym1",
                     "prop_spr_lt_win", "prop_spr_lt_win_ym1",
                     "prop_sum_lt_win", "prop_sum_lt_win_ym1",
                     "prop_sum_win_sq", "prop_sum_win_sq_ym1",
                     "maxan_mld", "maxan_mld_ym1",
                     "spwntmp_c")
flow_compiled <- sapply(vars_to_include,
                        function(x) tapply(get(x, flow_data),
                                           list(survey_data$system_id, survey_data$year_id),
                                           mean, na.rm = TRUE))
flow_compiled <- flow_compiled[!is.na(flow_compiled[, 1]), ]

# replace missing temperature data:
#   - King = Ovens (2011, 2012, 2017 filled with mean of adjacent years)
#   - Murray 1999-2002 filled with average of Murray 2003-2008
flow_compiled[system == 4 & year %in% c(1:4), "spwntmp_c"] <- mean(flow_compiled[system == 4 & year %in% c(5:10), "spwntmp_c"])
flow_compiled[system == 3, "spwntmp_c"] <- flow_compiled[match(paste0("5", year[system == 3]), paste0(system, year)), "spwntmp_c"]
flow_compiled[system == 3 & year %in% c(12:13), "spwntmp_c"] <- mean(flow_compiled[system == 3 & year %in% c(11, 14), "spwntmp_c"])
flow_compiled[system == 3 & year %in% c(18), "spwntmp_c"] <- mean(flow_compiled[system == 5 & year %in% c(17, 19), "spwntmp_c"])

# standardised flow data
flow_std <- apply(flow_compiled, 2, scale)
flow_scales <- apply(flow_compiled, 2, extract_standards)

# how many systems/years are we including?
n_system <- max(system)
n_year <- max(year)
n_flow <- ncol(flow_std)

# we need a system x age vector to flatten everything out
sys_vec <- rep(system, times = n_age)
age_vec <- rep(seq_len(n_age), each = length(system))
year_vec <- rep(year, times = n_age)
survey_vec <- rep(seq_along(system), times = n_age)
n_survey <- max(survey_vec)
sys_age_vec <- n_age * (sys_vec - 1) + age_vec
n_sys_age <- max(sys_age_vec)

# tricky bit: create a matrix of indices identifying cohorts
cohort_mat <- matrix(NA, nrow = length(system), ncol = n_age)
current_max <- 0
for (i in seq_len(n_system)) {
  sys_sub <- system == i
  year_sort <- year[sys_sub]
  cohort_tmp <- matrix(NA, nrow = sum(sys_sub), ncol = n_age)
  cohort_tmp[1, ] <- rev(seq_len(n_age))
  year_diffs <- diff(year_sort)
  for (j in seq_len(sum(sys_sub))[-1])
    cohort_tmp[j, ] <- cohort_tmp[j - 1, ] + year_diffs[j - 1]
  cohort_tmp <- cohort_tmp + current_max
  current_max <- max(cohort_tmp)
  cohort_mat[which(sys_sub)[order(year_sort)], ] <- cohort_tmp
}
cohort_vec <- c(cohort_mat)
n_cohort <- max(cohort_mat)

# need to set some variance priors for hierarchical coefficients
sigma_alpha <- normal(0, 10, truncation = c(0, Inf))
sigma_beta <- normal(0, 10, truncation = c(0, Inf))
sigma_year <- normal(0, 10, truncation = c(0, Inf))
sigma_cohort <- normal(0, 10, truncation = c(0, Inf))

# need priors on the regression coefs
alpha <- normal(0, sigma_alpha, dim = n_system)
beta <- normal(0, sigma_beta, dim = n_system, truncation = c(-Inf, 0))
gamma_year <- normal(0, sigma_year, dim = n_year)
gamma_cohort <- normal(0, sigma_cohort, dim = n_cohort)

# define flow priors (trickier than the others because we have crossed, shared variances)
sigma_flow_sub <- normal(0, 10, dim = c(n_age, n_flow), truncation = c(0, Inf))
sigma_flow <- do.call(rbind, lapply(seq_len(n_system), function(i) sigma_flow_sub))
flow_effects <- normal(0, sigma_flow, dim = c(n_sys_age, n_flow))

# define linear predictor: includes system and age specific flow effects, with random intercepts for year and cohort
mu <- alpha[sys_vec] + beta[sys_vec] * age_vec +
  gamma_cohort[cohort_vec] + gamma_year[year_vec] +
  rowSums(flow_effects[sys_age_vec, ] * flow_std[survey_vec, ])

# put back on original scale
modelled_ages <- exp(mu)
dim(modelled_ages) <- c(n_survey, n_age)

# calculate modelled sizes from ages
age_to_length <- t(length_to_age)
modelled_lengths <- modelled_ages %*% age_to_length

# flatten the modelled and observed lengths
length_vec <- c(modelled_lengths)
response_vec <- c(response_matrix)

# set likelihood for modelled sizes
distribution(response_vec) <- poisson(length_vec)

# set mcmc settings
nkeep <- 2500
nthin <- ifelse(nkeep > 1000, floor(nkeep / 1000), 1)

# compile and sample from model
mod <- model(alpha, beta, flow_effects, 
             length_to_age, modelled_ages,
             sigma_alpha, sigma_beta, sigma_flow,
             sigma_year, sigma_cohort,
             gamma_cohort, gamma_year)
draws <- mcmc(mod, n_samples = (2 * nkeep), warmup = nkeep, thin = nthin)

# summarise fitted model
mod_summary <- summary(draws)

# length-age conversions
a <- mod_summary$quantiles[grep("length_to_age", rownames(mod_summary$quantiles)), "50%"]
a <- matrix(a, nrow = n_len)
image(a, col = viridis::inferno(256))

# modelled ages and lengths
b <- mod_summary$quantiles[grep("modelled_ages", rownames(mod_summary$quantiles)), "50%"]
b <- matrix(b, nrow = nrow(modelled_ages))
d <- b %*% t(a)

# flow effects
beta_flow <- mod_summary$quantiles[grep("flow_effects", rownames(mod_summary$quantiles)), "50%"]
beta_flow <- array(beta_flow, dim = c(n_sys_age, n_flow))
beta_flow_lower <- mod_summary$quantiles[grep("flow_effects", rownames(mod_summary$quantiles)), "2.5%"]
beta_flow_lower <- array(beta_flow_lower, dim = c(n_sys_age, n_flow))
beta_flow_upper <- mod_summary$quantiles[grep("flow_effects", rownames(mod_summary$quantiles)), "97.5%"]
beta_flow_upper <- array(beta_flow_upper, dim = c(n_sys_age, n_flow))
colnames(beta_flow) <- vars_to_include
flow_system <- systems_to_keep[rep(seq_len(n_system), each = n_age)]
flow_age <- rep(seq_len(n_age), times = n_system)
