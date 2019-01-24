
# possibly zero-inflated

# sample RE is for unequal sample effort among years (could be year??)
## use a sample size offset in a Po model

# input indices
site ## site id of each obs
sample ## sample id of each obs
system ## system id of each obs

# length to age (vals for CT MC mod)
length_inf <- normal(150, 100, truncation = c(0, Inf))
time_zero <- normal(6, 100, truncation = c(0, Inf))
k_param <- exponential(1 / 0.0011)
c_param <- normal(-103, 100)

# define a function to convert size to age based on model in Todd & Koehn 2008
inverse_growth <- function(x, length_inf, time_zero, k_param, c_param) {
  
  par1 <- 1 / k_param
  ratio_param <- 1 - (time_zero / length_inf)
  par2 <- 1 / (1 - c_param * ratio_param)
  par3 <- (length_inf - time_zero) / (length_inf - x)
  par4 <- c_param * ratio_param
  
  par1 * log(par2 * (par3 - par4))
  
}
## ADD REs + priors
age_est <- inverse_growth(length_obs) +
  gamma_year + gamma_site + gamma_sample + gamma_system

# add likelihood for age model
distribution(age_obs) <- poisson(age_est)

# variance priors
### prob don't need in this form for NB mod
###sigma_main <- normal(0.0, 10.0, truncation = c(0, Inf))
sigma_offset <- normal(0.0, 10.0, truncation = c(0, Inf))
sigma_system <- normal(0.0, 10.0, truncation = c(0, Inf))
sigma_site <- normal(0.0, 10.0, truncation = c(0, Inf))

# priors for random effects
gamma_system <- normal(0.0, sigma_system)
gamma_site <- normal(0.0, sigma_site)

# calculate sample sizes
sample_size_offset <- normal(0.0, sigma_offset)

# catch ~ age + (1 | sample) + flow
age_pred <- inverse_growth(full_length_data)
age_count <- tf_hist(age_pred)
mu_catch_link <- beta_age[system] * age_count +
  sample_size_offset[sample] +
  gamma_system[system] +
  gamma_site[site]
mu_catch <- exp(mu_catch_link)

## alternative: check indexing
# multinom to convert modelled age to size
age_size <- dirichlet(alpha = matrix(1, n_age, n_obs))
# check dims and array index
distribution(age_size_obs) <- multinomial(size = apply(age_size_obs, 1, sum),
                                          prob = age_size)

## CHECK matmuls; should be able to keep as n_age x n_obs matrix
age_count <- alpha + beta[system] * age_index + gamma
size_count <- age_size %*% age_count
distribution(size_obs) <- poisson(size_count)
## THIS IS THE CLEANEST APPROACH; truncate ages/sizes.

## ALTERNATIVE2: try a poisson process model
n_int <- 200
max_age <- 40
## need to add linear predictor to rates;
age_vec <- inverse_growth(size_vec)
linear_predictor <- alpha + beta[system] * age_vec + gammas
p_est <- ilogit(linear_predictor)
distribution(ones(length(size_vec))) <- binomial(size = 1, prob = p_est)
## EXAMPLE CODE FOR THIS
integration_sizes <- seq(0, max_age, len = n_int)
binsize <- abs(diff(range(integration_sizes))) / n_int
weights <- rep(binsize, n_int)
df_ppm <- data.frame(sizes = c(sizes, integration_sizes),
                     response = rep(1:0, c(n, n_int)),
                     offset = rep(c(1, binsize), c(n, n_int)))
lambda <- exp(alpha + beta * df_ppm$sizes)  # add offset?
lambda_p <- ilogit(alpha + beta * df_ppm$sizes)  # add offset?
distribution(df_ppm$response) <- poisson(lambda)
distribution(df_ppm$response) <- binomial(size = 1, prob = lambda_p)
## PREDICTOR HERE IS SIZE OCNVERTED TO AGE, not sizes. Can't have age as
## a response because it's unobserved discrete.

# possibly need to account for zero-inflation but probably not: check residuals
# if so, look into mixture() distribution, with a binomial and poisson component

### add CJS for detection probab
# recapture data
catch_size <- catch_size[apply(catch_size, 1, function(x) sum(x > 0)) > 1, ]
catch_binary <- ifelse(catch_size > 0, 1, 0)

# extract summary info from data
first_obs <- apply(catch_binary, 1, function(x) min(which(x > 0)))
final_obs <- apply(catch_binary, 1, function(x) max(which(x > 0)))
obs_id <- apply(catch_binary, 1, function(x) seq(min(which(x > 0)), max(which(x > 0)), by = 1)[-1])
obs_id <- unlist(obs_id)
capture_vec <- apply(catch_binary, 1, function(x) x[min(which(x > 0)):max(which(x > 0))][-1])
capture_vec <- unlist(capture_vec)
n_time <- ncol(catch_binary)

# priors
detection_mean <- beta(1, 1, dim = nyear)
phi <- beta(1, 1, dim = nyear)
### MAKE THESE f(flow)
# phi <- ilogit(alpha_phi + flow_surv %*% beta_phi)

# derived parameter
chi <- ones(n_time)
for (i in seq_len(n_time - 1)) {
  tn <- n_time - i
  chi[tn] <- (1 - phi[tn]) + phi[tn] * (1 - detection_mean[tn + 1]) * chi[tn + 1]
}

# dummy variables
alive_data <- ones(length(obs_id))            # definitely alive
not_seen_last <- final_obs != ncol(catch_binary) # ignore observations in last timestep
final_observation <- ones(sum(not_seen_last)) # final observation

# set CJS likelihoods
distribution(alive_data) <- bernoulli(phi[obs_id - 1])
distribution(capture_vec) <- bernoulli(detection_mean[obs_id])
distribution(final_observation) <- bernoulli(chi[final_obs[not_seen_last]])  

# set other likelihoods
catch_data <- hist(full_length_data)$counts
distribution(catch_data) <- poisson(mu_catch)

# OR
# transform mu + sigma to NB params
nb_size <- (mu_catch * mu_catch) / (sigma_main * sigma_main - mu_catch)
nb_prob <- 1 - (mu_catch / (sigma_main * sigma_main))
# distribution(catch_data) <- negative_binomial(nb_size, nb_prob)

# abundance estimates
abund <- mu_catch / detection_mean[year_id]

# trying to fit this:
# class_abunds ~ age_classes + covars + offset(sum(class_abunds)) + REs
# Now have multiple systems

# compile model
mod <- model()

# mcmc settings
n_warmup <- 1000
n_samples <- 1000
n_chains <- 4
n_thin <- 1

# initialise?
init_list <- vector("list", length = n_chains)
for (i in seq_len(n_chains)) {
  opt_est <- opt(mod, optimiser = adam(), max_iterations = 1000)
  init_list[[i]] <- initials()
}

# sample from model
draws <- mcmc(mod,
              sampler = hmc(Lmax = 50),
              n_samples = n_samples, warmup = n_warmup,
              chains = n_chains, thin = n_thin,
              initial_values = init_list)
  
# summarise fitted
mod_summary <- summary(draws)