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

# need to load all otolith data
oti_data <- readRDS("data/otolith-data-loaded.rds")

# filter to MC
oti_data <- oti_data[oti_data$SPECIES == "Maccullochella peelii", ]

# need to load flow data
flow_data <- readRDS("data/flow-data-loaded.rds")

# filter to MC
flow_data <- flow_data[to_keep, ]

# define model
oti_analysis_data <- data.frame(length = oti_data$T_Length..mm.,
                                age = oti_data$AGE,
                                site = as.integer(as.factor(oti_data$Site)),
                                year = as.integer(as.factor(oti_data$Year)))
oti_analysis_data <- oti_analysis_data[apply(oti_analysis_data, 1, function(x) !any(is.na(x))), ]
oti_analysis_data$site <- as.integer(as.factor(oti_analysis_data$site))
oti_analysis_data$year <- as.integer(as.factor(oti_analysis_data$year))

# length to age (vals for CT MC mod)
len_par <- normal(150, 50, truncation = c(0, Inf))
time_par <- normal(6, 50, truncation = c(0, Inf))
k_par <- exponential(1 / 0.0011)
c_par <- normal(-103, 50)

# define a function to convert size to age based on model in Todd & Koehn 2008
inverse_growth <- function(x, length_inf, time_zero, k_param, c_param) {
  
  par1 <- 1 / k_param
  ratio_param <- 1 - (time_zero / length_inf)
  par2 <- 1 / (1 - c_param * ratio_param)
  par3 <- (length_inf - time_zero) / (length_inf - x)
  par4 <- c_param * ratio_param
  
  par1 * log(par2 * (par3 - par4))

}

# random intercepts
sigma_site_oti <- normal(0, sd(oti_analysis_data$age), truncation = c(0, Inf))
sigma_year_oti <- normal(0, sd(oti_analysis_data$age), truncation = c(0, Inf))
gamma_site_oti <- normal(0, sigma_site_oti, dim = length(unique(oti_analysis_data$site)))
gamma_year_oti <- normal(0, sigma_year_oti, dim = length(unique(oti_analysis_data$year)))

# linear predictor
age_est <- inverse_growth(oti_analysis_data$length / 10,
                          len_par, time_par, k_par, c_par) +
  gamma_site_oti[oti_analysis_data$site] + gamma_year_oti[oti_analysis_data$year]

# add likelihood for age model
sigma_main <- normal(0, 10, truncation = c(0, Inf))
distribution(oti_analysis_data$age) <- normal(age_est, sigma_main)

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
age_vec <- inverse_growth(survey_data$length,
                          len_par, time_par, k_par, c_par)

# setup PPM as a GLM
n_int <- 100
max_age <- 50

# expand random effect data to include integration points
system_expanded <- c(survey_data$system, rep(max(survey_data$system) + 1, n_int))
site_expanded <- c(survey_data$site, rep(max(survey_data$site) + 1, n_int))
year_expanded <- c(survey_data$year, rep(max(survey_data$year) + 1, n_int))
dataset_expanded <- c(survey_data$dataset, rep(max(survey_data$dataset) + 1, n_int))

# priors for PPM
alpha_age <- normal(0, 1)
beta_age <- normal(0, 1)
sigma_system <- normal(0, 1, truncation = c(0, Inf))
sigma_site <- normal(0, 1, truncation = c(0, Inf))
sigma_year <- normal(0, 1, truncation = c(0, Inf))
sigma_dataset <- normal(0, 1, truncation = c(0, Inf))
gamma_system <- normal(0, sigma_system, dim = max(system_expanded))
gamma_site <- normal(0, sigma_site, dim = max(site_expanded))
gamma_year <- normal(0, sigma_year, dim = max(year_expanded))
gamma_dataset <- normal(0, sigma_dataset, dim = max(dataset_expanded))

# setup integration points for PPM
integration_ages <- exp(seq(log(1e-2), log(max_age), length = n_int))
eps <- .Machine$double.eps
binsize <- diff(c(0, integration_ages))

# expand age data to include integration points
age_expanded <- c(age_vec, integration_ages)

# setup linear predictor and response variable
response_vec <- rep(1:0, c(length(age_vec), n_int))
offset <- c(rep(eps, length(age_vec)), binsize)
lambda <- exp(log(offset) + alpha_age + beta_age * age_expanded) # +
                # gamma_system[system_expanded] * age_expanded +
                # gamma_site[site_expanded] +
                # gamma_year[year_expanded] * age_expanded + 
                # gamma_dataset[dataset_expanded])

# set likelihood
distribution(response_vec) <- poisson(lambda)

# compile model
mod <- model(len_par, time_par, k_par, c_par,
             age_vec,
             alpha_age, beta_age,
             gamma_system, gamma_site, gamma_year, gamma_dataset,
             sigma_system, sigma_site, sigma_year, sigma_dataset)

# sample from model
init_set <- initials(len_par = 150, time_par = 6, k_par = 0.001, c_par = -100)
draws <- mcmc(mod, initial_values = init_set,
              n_samples = 1000, warmup = 1000)

# summarise model outputs
mod_summary <- summary(draws)

# pull out values of interest
age_vec_est <- mod_summary$quantiles[grep("age_vec", rownames(mod_summary$quantiles)), ]
alpha_est <- mod_summary$quantiles[grep("alpha_age", rownames(mod_summary$quantiles)), ]
beta_est <- mod_summary$quantiles[grep("beta_age", rownames(mod_summary$quantiles)), ]
sys_est <- mod_summary$quantiles[grep("gamma_system", rownames(mod_summary$quantiles)), ]
site_est <- mod_summary$quantiles[grep("gamma_system", rownames(mod_summary$quantiles)), ]
yr_est <- mod_summary$quantiles[grep("gamma_system", rownames(mod_summary$quantiles)), ]
ds_est <- mod_summary$quantiles[grep("gamma_system", rownames(mod_summary$quantiles)), ]

# predict fitted curves at integration points
p <- exp(log(binsize) + alpha_est[3] + beta_est[3] * integration_ages +
           sys_est[1, 3] + site_est[1, 3] + yr_est[1, 3] + ds_est[1, 3])
hist(age_vec_est[, 3], breaks = c(-0.5, integration_ages),
     col = grey(0.6), border = NA, las = 1,
     xlab = "Age", ylab = "Density")
lines(p ~ integration_ages, lwd = 3)
