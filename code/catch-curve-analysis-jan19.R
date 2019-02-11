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
age_alpha <- normal(0, 10)
age_beta <- normal(0, 10)

# define a function to convert size to age based on model in Todd & Koehn 2008
inverse_growth <- function(x, int, slope) {
  
  int + slope * x

}

# random intercepts
sigma_site_oti <- normal(0, sd(oti_analysis_data$age), truncation = c(0, Inf))
sigma_year_oti <- normal(0, sd(oti_analysis_data$age), truncation = c(0, Inf))
gamma_site_oti <- normal(0, sigma_site_oti, dim = length(unique(oti_analysis_data$site)))
gamma_year_oti <- normal(0, sigma_year_oti, dim = length(unique(oti_analysis_data$year)))

# linear predictor
length_scaled <- scale(oti_analysis_data$length)
age_est <- inverse_growth(length_scaled,
                          age_alpha, age_beta)

# add likelihood for age model
sigma_oti <- normal(0, 10, truncation = c(0, Inf))
distribution(oti_analysis_data$age) <- normal(age_est, sigma_oti, truncation = c(0, Inf))

# data prep
survey_data <- data.frame(length = alldat$totallength,
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
survey_length_scaled <- (survey_data$length - attributes(length_scaled)$`scaled:center`) /
  attributes(length_scaled)$`scaled:scale`
age_vec <- inverse_growth(survey_length_scaled,
                          age_alpha, age_beta)

# setup PPM as a GLM
n_int <- 20
max_age <- 60

# pull out indices fo rrandom effects
nsystem <- max(survey_data$system)
nsite <- max(survey_data$site)
nyear <- max(survey_data$year)
ndataset <- max(survey_data$dataset)

# priors for PPM
alpha_age <- normal(0, 10)
beta_age <- normal(0, 10)

# variance priors for random effects
sigma_system <- normal(0, 1, truncation = c(0, Inf))
sigma_site <- normal(0, 1, truncation = c(0, Inf))
sigma_year <- normal(0, 1, truncation = c(0, Inf))
sigma_dataset <- normal(0, 1, truncation = c(0, Inf))

# main priors for random effects (expanded with an extra zero for integration points)
gamma_system <- zeros(nsystem + 1)
gamma_system[seq_len(nsystem)] <- normal(0, sigma_system, dim = nsystem)
gamma_site <- zeros(nsite + 1)
gamma_site[seq_len(nsite)] <- normal(0, sigma_site, dim = nsite)
gamma_year <- zeros(nyear + 1)
gamma_year[seq_len(nyear)] <- normal(0, sigma_year, dim = nyear)
gamma_dataset <- zeros(ndataset + 1)
gamma_dataset[seq_len(ndataset)] <- normal(0, sigma_dataset, dim = ndataset)

# setup integration points for PPM
integration_ages <- seq(0, max_age, length = n_int)
eps <- .Machine$double.eps
binsize <- diff(range(integration_ages)) / n_int

# expand age data to include integration points
age_expanded <- c(age_vec, integration_ages)
system_expanded <- c(survey_data$system, rep(nsystem + 1, n_int))
site_expanded <- c(survey_data$site, rep(nsite + 1, n_int))
year_expanded <- c(survey_data$year, rep(nyear + 1, n_int))
dataset_expanded <- c(survey_data$dataset, rep(ndataset + 1, n_int))

# setup linear predictor and response variable
response_vec <- rep(1:0, c(length(age_vec), n_int))
offset <- c(rep(eps, length(age_vec)), rep(binsize, n_int))
mu <- alpha_age + beta_age * age_expanded # +
  # gamma_system[system_expanded] * age_expanded +
  # gamma_site[site_expanded] +
  # gamma_year[year_expanded] * age_expanded +
  # gamma_dataset[dataset_expanded]

## WHY CAN'T WE HAVE hist(age_vec) = f(age)? (there was a reason)

# need to add offset and exponentiate linear predictor
lambda <- exp(log(offset) + mu)

# set likelihood
distribution(response_vec) <- poisson(lambda)

# compile model
mod <- model(len_par, time_par, k_par, c_par,
             age_vec,
             sigma_oti,
             alpha_age, beta_age)#,
             # gamma_system, gamma_site, gamma_year, gamma_dataset,
             # sigma_system, sigma_site, sigma_year, sigma_dataset)

# sample from modell
# init_set <- initials(len_par = 150, time_par = 6, k_par = 0.001, c_par = -100)
init_set <- initials(age_alpha = 6.0, age_beta = 3.0)
# opt_est <- opt(mod, max_iterations = 500, tolerance = 1e-8,
#                initial_values = init_set)
# init_set <- initials(len_par = opt_est$par$len_par,
#                      time_par = opt_est$par$time_par,
#                      k_par = opt_est$par$k_par,
#                      c_par = opt_est$par$c_par,
#                      alpha_age = opt_est$par$alpha_age,
#                      beta_age = opt_est$par$beta_age)
draws <- mcmc(mod, initial_values = init_set,
              n_samples = 4000, warmup = 4000)

# summarise model outputs
mod_summary <- summary(draws)

# pull out values of interest
age_vec_est <- mod_summary$quantiles[grep("age_vec", rownames(mod_summary$quantiles)), ]
alpha_est <- mod_summary$quantiles[grep("alpha_age", rownames(mod_summary$quantiles)), ]
beta_est <- mod_summary$quantiles[grep("beta_age", rownames(mod_summary$quantiles)), ]
alpha_oti_est <- mod_summary$quantiles[grep("age_alpha", rownames(mod_summary$quantiles)), ]
beta_oti_est <- mod_summary$quantiles[grep("age_beta", rownames(mod_summary$quantiles)), ]
# sys_est <- mod_summary$quantiles[grep("gamma_system", rownames(mod_summary$quantiles)), ]
# site_est <- mod_summary$quantiles[grep("gamma_site", rownames(mod_summary$quantiles)), ]
# yr_est <- mod_summary$quantiles[grep("gamma_year", rownames(mod_summary$quantiles)), ]
# ds_est <- mod_summary$quantiles[grep("gamma_dataset", rownames(mod_summary$quantiles)), ]

# predict fitted curves at integration points
p <- exp(alpha_est[3] + beta_est[3] * integration_ages)
hist(age_vec_est[, 3], #breaks = c(-0.5, integration_ages),
     col = grey(0.6), border = NA, las = 1,
     xlab = "Age", ylab = "Density", freq = FALSE)
lines(c(p * binsize / sum(integration_ages)) ~ integration_ages, lwd = 3)
