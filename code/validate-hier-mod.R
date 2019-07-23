# validate full hier model

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
             length_to_age,
             gamma_cohort, gamma_year)
draws <- mcmc(mod, n_samples = (2 * nkeep), warmup = nkeep, thin = nthin)

# summarise fitted model
mod_summary <- summary(draws)

# length-age conversions
convert <- mod_summary$quantiles[grep("length_to_age", rownames(mod_summary$quantiles)), "50%"]
convert <- matrix(convert, nrow = n_len)

# need new_sys, new_cohort, new_year, new_sys_age, flow_new_std
##  can drop cohort/year if we think it's unimportant (or use to predict?)
##  Need to match indices, replace with 0 if missing a level
mu_pred <- alpha_est[new_sys] + beta_est[new_sys] * age_vec +
  cohort_est[new_cohort] + year_est[new_year] +
  rowSums(flow_effects[new_sys_age, ] * flow_new_std)

ages_pred <- exp(mu_pred)
ages_pred <- matrix(ages_pred, nrow = n_pred)

lengths_pred <- ages_pred %*% t(convert)

