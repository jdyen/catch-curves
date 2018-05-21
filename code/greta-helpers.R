# greta scripts for catch curve analysis

# create separable covariance function from two component covariances
create_greta_covar <- function(n_sp, n_age,
                               eta_sp = 1, eta_age = 4,
                               sigma_mean = 0.0, sigma_sd = 1.0) {
  
  if (length(sigma_mean) == 1)
    sigma_mean <- rep(sigma_mean, 2)
  if (length(sigma_sd) == 1)
    sigma_sd <- rep(sigma_sd, 2)
  
  species_corr_est <- lkj_correlation(eta = eta_sp, dim = n_sp)
  age_corr_est <- lkj_correlation(eta = eta_age, dim = n_age)
  species_sigma <- lognormal(mean = sigma_mean[1], sd = sigma_sd[1], dim = n_sp)
  age_sigma <- lognormal(mean = sigma_mean[2], sd = sigma_sd[2], dim = n_age)
  species_covar_est <- (species_sigma %*% t(species_sigma)) * species_corr_est
  age_covar_est <- age_sigma %*% t(age_sigma) * age_corr_est
  covar_est <- kronecker(species_covar_est, age_covar_est)
  
  covar_est
  
}

# set priors for covariate model
prepare_greta_model <- function(stage_data, flow_data, info_data, n_sp, n_class) {

  # calculate length of full response vector
  n_class_sp <- (n_sp * n_class)
  
  # pull out n_obs
  n_obs <- nrow(stage_data)
  
  # create a separable covariance matrix with spp and age components
  covar_mat <- create_greta_covar(n_sp = n_sp,
                                  n_age = n_class)

  # set priors
  alpha_mean <- normal(mean = 0.0, sd = 1.0, dim = 1)
  beta_mean <- normal(mean = 0.0, sd = 1.0, dim = 1)
  alpha <- normal(mean = alpha_mean, sd = sigma_alpha, dim = n_sp)
  beta <- normal(mean = beta_mean, sd = sigma_beta, dim = n_sp)
  gamma <- normal(mean = 0.0, sd = sigma_gamma, dim = n_sp)
  delta <- normal(mean = 0.0, sd = sigma_delta, dim = n_class)
  kappa <- normal(mean = 0.0, sd = sigma_kappa, dim = n_class_sp)
  # set linear predictor
  # need to include year, reach (this must be river*reach, i.e., reach within river),
  # river as random effects
  # need a size * age random effect (I think)
  # probably need a year * reach random effect (or year * river)
  
  #
  # there will be a nsp * n_age vector of means
  # this will need to be swept through a matrix that is the
  # outer product of sp*age regression coefficients with 
  # the flow data

  # use centered MVN model 
  # set likelihoods
  # see if can setup a Stan-style normal() call with 
  #   work done on the covariances outside of this -- will be
  #   much faster and easier.
  # otherwise will need to loop through multivariate_normal calls (although
  #   check how nick is doing this with JSDMs)
  resid <- multivariate_normal(mean = rep(0, n_class_sp), Sigma = covar_mat, dim = n_obs)

  #mod <- model()
  
  mod
    
}
