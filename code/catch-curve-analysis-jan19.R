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
sigma_site <- normal(0, sd(oti_analysis_data$age), truncation = c(0, Inf))
sigma_year <- normal(0, sd(oti_analysis_data$age), truncation = c(0, Inf))
gamma_site <- normal(0, sigma_site, dim = length(unique(oti_analysis_data$site)))
gamma_year <- normal(0, sigma_year, dim = length(unique(oti_analysis_data$year)))

# linear predictor
age_est <- inverse_growth(oti_analysis_data$length / 10,
                          len_par, time_par, k_par, c_par) +
  gamma_site[oti_analysis_data$site] + gamma_year[oti_analysis_data$year]

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

# linear predictor
age_vec <- inverse_growth(survey_data$length,
                          len_par, time_par, k_par, c_par)
      # + REs gamma_system + gamma_site + gamma_year + gamma_dataset

# priors for PPM
# alpha_age <- normal(0, 100)
# beta_age <- normal(0, 100)
# 
# # setup PPM as a GLM
# n_int <- 200
# max_age <- 20
# integration_ages <- seq(0, max_age, len = n_int)
# binsize <- abs(diff(range(integration_ages))) / n_int
# age_expanded <- c(age_vec, integration_ages)
# response_vec <- rep(1:0, c(length(age_vec), n_int))
# offset <- rep(c(1, binsize), c(length(age_vec), n_int))
# lambda <- exp(log(offset) + alpha_age + beta_age * age_expanded)
# distribution(response_vec) <- poisson(lambda)

# compile model
mod <- model(len_par, time_par, k_par, c_par,
             age_vec)#,
             # lambda,
             # sigma_main, sigma_site, sigma_year)

# sample from model
init_set <- initials(len_par = 150, time_par = 6, k_par = 0.001, c_par = -100)
draws <- mcmc(mod, initial_values = init_set)

# summarise model outputs
mod_summary <- summary(draws)

