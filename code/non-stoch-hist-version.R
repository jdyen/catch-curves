# set working directory
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(greta)

# source helper functions
source("code/greta-helpers.R")

# load compiled survey data
alldat <- readRDS("data/data-loaded.rds")

# filter survey data to MC
to_keep <- alldat$Scientific.Name == "Maccullochella peelii"
alldat <- alldat[to_keep, ]
alldat$Common.Name <- NULL

# need to load flow data
flow_data <- readRDS("data/flow-data-loaded.rds")

# filter to MC
flow_data <- flow_data[to_keep, ]

# define a function to convert size to age based on model in Todd & Koehn 2008
inverse_growth <- function(x, length_inf, time_zero, k_param, c_param) {
  
  par1 <- 1 / k_param
  ratio_param <- 1 - (time_zero / length_inf)
  par2 <- 1 / (1 - c_param * ratio_param)
  par3 <- (length_inf - time_zero) / (length_inf - x)
  par4 <- c_param * ratio_param
  
  par1 * log(par2 * (par3 - par4))
  
}

# data prep
survey_data <- data.frame(length = alldat$totallength / 10,
                          system = alldat$SYSTEM,
                          site = alldat$SITE_CODE,
                          year = alldat$YEAR,
                          dataset = alldat$dataset)
survey_data <- survey_data[apply(survey_data, 1, function(x) !any(is.na(x))), ]
survey_data$system <- as.integer(as.factor(survey_data$system))
survey_data$site <- as.integer(as.factor(survey_data$site))
survey_data$year <- as.integer(as.factor(survey_data$year))
survey_data$dataset <- as.integer(as.factor(survey_data$dataset))

# convert observed lengths to ages
len_par <- 150
time_par <- 6
k_par <- 0.0011
c_par <- -103
age_vec <- inverse_growth(survey_data$length,
                          len_par, time_par, k_par, c_par)
age_vec[age_vec < 0] <- 0

# pull out indices for random effects
nsystem <- max(survey_data$system)
nsite <- max(survey_data$site)
nyear <- max(survey_data$year)
ndataset <- max(survey_data$dataset)

# settings for linear model
n_int <- 100
max_age <- ceiling(max(age_vec))

# priors for PPM
alpha_age <- normal(0, 10)
beta_age <- normal(0, 10)

# variance priors for random effects
sigma_system <- normal(0, 1, truncation = c(0, Inf))
sigma_site <- normal(0, 1, truncation = c(0, Inf))
sigma_year <- normal(0, 1, truncation = c(0, Inf))
sigma_dataset <- normal(0, 1, truncation = c(0, Inf))
# sigma_delta_sys <- normal(0, 1, truncation = c(0, Inf))
# sigma_delta_yr <- normal(0, 1, truncation = c(0, Inf))

# main priors for random effects (expanded with an extra zero for integration points)
gamma_system <- normal(0, sigma_system, dim = nsystem)
gamma_site <- normal(0, sigma_site, dim = nsite)
gamma_year <- normal(0, sigma_year, dim = nyear)
gamma_dataset <- normal(0, sigma_dataset, dim = ndataset)
# delta_sys <- zeros(nsystem + 1)
# delta_sys[seq_len(nsystem)] <- normal(0, sigma_delta_sys, dim = nsystem)
# delta_yr <- zeros(nyear + 1)
# delta_yr[seq_len(nyear)] <- normal(0, sigma_delta_yr, dim = nyear)

# we need binned data by site and year
age_seq <- seq(0, max_age, by = 1)
hist_fn <- function(x, breaks) {
  hist(x, breaks = breaks, plot = FALSE)$counts
}
age_counts <- tapply(age_vec, list(survey_data$system, survey_data$year), hist_fn, breaks = age_seq)
age_mat <- do.call(rbind, c(age_counts))
site_info <- rep(rownames(age_counts), times = nyear)
year_info <- rep(colnames(age_counts), each = nsystem)
to_keep <- !sapply(c(age_counts), is.null)
system_info <- as.numeric(site_info[to_keep])
year_info <- as.numeric(year_info[to_keep])
# system_info <- survey_data$system[match(site_info, survey_data$site)]
# dataset_info <- survey_data$dataset[match(site_info, survey_data$site)]

# now we need to define a linear predictor
response_vec <- c(age_mat)
mu <- alpha_age + beta_age * age_seq[-length(age_seq)]
mu <- rep(mu, times = nrow(age_mat))
mu <- mu + rep(gamma_system[system_info], each = max_age) +
  rep(gamma_year[year_info], each = max_age) # +
  # rep(gamma_site[site_info], each = max_age) +
  # rep(gamma_dataset[dataset_info], each = max_age)

# need to add offset and exponentiate linear predictor
lambda <- exp(mu)

# set likelihood
distribution(response_vec) <- poisson(lambda)

# compile model
mod <- model(alpha_age, beta_age,
             gamma_system, gamma_year, # gamma_site, gamma_dataset,
             sigma_system, sigma_year) #sigma_site, sigma_dataset)

# sample from modelht
draws <- mcmc(mod, n_samples = 5000, warmup = 5000)

# summarise model outputs
mod_summary <- summary(draws)

# pull out values of interest
age_vec_est <- age_vec
alpha_est <- mod_summary$quantiles[grep("alpha_age", rownames(mod_summary$quantiles)), ]
beta_est <- mod_summary$quantiles[grep("beta_age", rownames(mod_summary$quantiles)), ]
delta_est <- mod_summary$quantiles[grep("delta_age", rownames(mod_summary$quantiles)), ]
system_est <- mod_summary$quantiles[grep("gamma_system", rownames(mod_summary$quantiles)), ]
# site_est <- mod_summary$quantiles[grep("gamma_site", rownames(mod_summary$quantiles)), ]
year_est <- mod_summary$quantiles[grep("gamma_year", rownames(mod_summary$quantiles)), ]
# dataset_est <- mod_summary$quantiles[grep("gamma_dataset", rownames(mod_summary$quantiles)), ]
# delta_sys_est <- mod_summary$quantiles[grep("delta_sys", rownames(mod_summary$quantiles)), ]
# delta_yr_est <- mod_summary$quantiles[grep("delta_yr", rownames(mod_summary$quantiles)), ]

# predict fitted curves at integration points
p <- exp(alpha_est[3] + (beta_est[3]) * integration_ages + delta_est[3] * integration_ages * integration_ages)
# site_est[nsite, 3] +
# dataset_est[ndataset, 3])
hist(age_vec_est, breaks = c(-0.5, integration_ages),
     col = grey(0.6), border = NA, las = 1,
     xlab = "Age", ylab = "Density", freq = FALSE)
lines(c(p * binsize / sum(integration_ages)) ~ integration_ages, lwd = 3)
