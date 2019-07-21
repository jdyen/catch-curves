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

# pull out necessary indices
n_system <- max(survey_data$system_id)
system <- c(tapply(survey_data$system_id,
                   list(survey_data$system, survey_data$year),
                   unique))
system <- system[!is.na(system)]
year <- c(tapply(survey_data$year_id,
                   list(survey_data$system, survey_data$year),
                   unique))
year <- year[!is.na(year)]

# create a linear predictor for modelled ages
sigma_system <- normal(0, 10, truncation = c(0, Inf))
sigma_system2 <- normal(0, 10, truncation = c(0, Inf))
# sigma_system <- t(do.call(cbind, lapply(seq_len(n_system),
#                                         function (x) rep(normal(0, 1, truncation = c(0, Inf)), n_age))))
# sigma_age <- do.call(cbind, lapply(seq_len(n_age),
#                                    function (x) rep(normal(0, 1, truncation = c(0, Inf)), n_system)))
# sigma_alpha <- sigma_system + sigma_age

# alpha <- normal(zeros(n_system, n_age), sigma_alpha, dim = c(n_system, n_age))
alpha <- normal(0, sigma_system, dim = n_system)
beta_sub <- normal(0, sigma_system2, dim = n_system, truncation = c(-Inf, 0))
beta <- t(do.call(cbind, lapply(system, function(i) beta_sub[i] * seq_len(n_age))))
mu <- sweep(beta[system, ], 1, alpha[system], "+") # + X %*% beta_matrix + gamma_year[year]

## COUlD INFER COHORT AS YEAR - EST_COHORT??

## COULD TRY POINT PROCESS VERSION WITH SIMILAR NEGATIVE CONSTRAINT ON BETA?

modelled_ages <- exp(mu)

# calculate modelled sizes from ages
age_to_length <- t(length_to_age)
modelled_lengths <- modelled_ages %*% age_to_length

# flatten the modelled and observed lengths
length_vec <- c(modelled_lengths)
response_vec <- c(response_matrix)

# set likelihood for modelled sizes
distribution(response_vec) <- poisson(length_vec)

# compile and sample from model
mod <- model(alpha, sigma_system, sigma_system2, length_to_age, modelled_ages)
draws <- mcmc(mod, n_samples = 5000, warmup = 2500, thin = 1)

# summarise fitted model
mod_summary <- summary(draws)
a <- mod_summary$quantiles[grep("length_to_age", rownames(mod_summary$quantiles)), "50%"]
a <- matrix(a, nrow = n_len)

image(a, col = viridis::inferno(256))


b <- mod_summary$quantiles[grep("modelled_ages", rownames(mod_summary$quantiles)), "50%"]
b <- matrix(b, nrow = nrow(mu))
  
  
