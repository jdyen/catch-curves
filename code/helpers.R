
# define leslie matrix function
leslie_matrix <- function(n_age, density_dependence, predictors, params = list()) {
  
  # initialise model parameters
  param_list <- list(dens_lower = 0,
                     dens_upper = 10,
                     fec_stages = n_age,
                     n_site = 1,
                     surv_mean = 0.0,
                     surv_sd = 10.0,
                     fec_mean = 0.0,
                     fec_sd = 10.0,
                     init_mean = 0.0,
                     init_sd = 2.0)
  param_list[names(params)] <- params
  
  # how many observations and predictor variables?
  n_obs <- nrow(predictors)
  n_pred <- ncol(predictors)
  n_fec <- length(param_list$fec_stages)
  
  # initialise effects of predictors on survival
    ## SET UP hierarchical priors on these
  coefs_surv <- normal(mean = param_list$surv_mean, sd = param_list$surv_sd, dim = c(n_pred, n_age))
  
  # initialise effects of predictors on fecundity
   ## SET UP hierarchical priors on these
  coefs_fec <- normal(mean = param_list$fec_mean, sd = param_list$fec_sd, dim = c(n_pred, n_fec)) 
  
  # add additional noise to survival and fecundity parameters
  sigma_year_surv <- normal(0.0, 10.0, truncation = c(0, Inf))
  gamma_year_surv <- normal(0.0, sigma_year_surv, dim = c(n_obs, n_age))
  sigma_year_fec <- normal(0.0, 10.0, truncation = c(0, Inf))
  gamma_year_fec <- normal(0.0, sigma_year_fec, dim = c(n_obs, n_fec))
  
  # setup linear predictors
  surv_params <- ilogit(predictors %*% coefs_surv + gamma_year_surv)
  fec_params <- exp(predictors %*% coefs_fec + gamma_year_fec)

  # construct leslie matrices
  leslie_matrix <- zeros(n_obs, n_age, n_age)
  for (i in seq_len(n_age - 1)) {
    leslie_matrix[, i + 1, i] <- surv_params[, i]
  }
  leslie_matrix[, n_age, n_age] <- surv_params[, n_age]
  leslie_matrix[, param_list$fec_stages, 1] <- fec_params
  
  # initialise abundances
  init <- lognormal(meanlog = param_list$init_mean, sdlog = param_list$init_sd,
                    dim = c(n_age, n_site))
  
  # create a parameter to capture density dependence
  density <- uniform(min = param_list$dens_lower,
                     max = param_list$dens_upper,
                     dim = n_site)
  
  # collate and return outputs  
  list(matrix = leslie_matrix,
       init = init,
       dens_param = density,
       coefs_surv = coefs_surv,
       coefs_fec = coefs_fec,
       surv_params = surv_params,
       fec_params = fec_params,
       gamma_year_surv = gamma_year_surv,
       gamma_year_fec = gamma_year_fec)
  
}

# define leslie matrix function
lefkovitch_matrix <- function(n_stage, density_dependence, predictors, params = list()) {
  
  # initialise model parameters
  param_list <- list(dens_lower = 0,
                     dens_upper = 10,
                     fec_stages = n_stage,
                     n_site = 1,
                     surv_mean = 0.0,
                     surv_sd = 10.0,
                     growth_mean = 0.0,
                     growth_sd = 10.0,
                     fec_mean = 0.0,
                     fec_sd = 10.0,
                     init_mean = 0.0,
                     init_sd = 2.0)
  param_list[names(params)] <- params
  
  # how many observations and predictor variables?
  n_obs <- nrow(predictors)
  n_pred <- ncol(predictors)
  n_fec <- length(param_list$fec_stages)
  
  # initialise effects of predictors on vital rates
  coefs_surv <- normal(mean = param_list$surv_mean, sd = param_list$surv_sd, dim = c(n_pred, n_stage))
  coefs_growth <- normal(mean = param_list$growth_mean, sd = param_list$growth_sd, dim = c(n_pred, n_stage - 1))
  coefs_fec <- normal(mean = param_list$fec_mean, sd = param_list$fec_sd, dim = c(n_pred, n_fec)) 
  
  # add additional noise to survival and fecundity parameters
  sigma_year_surv <- normal(0.0, 10.0, truncation = c(0, Inf))
  gamma_year_surv <- normal(0.0, sigma_year_surv, dim = c(n_obs, n_stage))
  sigma_year_growth <- normal(0.0, 10.0, truncation = c(0, Inf))
  gamma_year_growth <- normal(0.0, sigma_year_surv, dim = c(n_obs, n_stage - 1))
  sigma_year_fec <- normal(0.0, 10.0, truncation = c(0, Inf))
  gamma_year_fec <- normal(0.0, sigma_year_fec, dim = c(n_obs, n_fec))
  
  # setup linear predictors
  surv_params <- ilogit(predictors %*% coefs_surv + gamma_year_surv)
  growth_params <- ilogit(predictors %*% coefs_growth + gamma_year_growth)
  fec_params <- exp(predictors %*% coefs_fec + gamma_year_fec)
  
  # construct leslie matrices
  lefkovitch_matrix <- zeros(n_obs, n_stage, n_stage)
  for (i in seq_len(n_stage - 1)) {
    lefkovitch_matrix[, i, i] <- surv_params[, i]
    lefkovitch_matrix[, i + 1, i] <- growth_params[, i]
  }
  lefkovitch_matrix[, n_stage, n_stage] <- surv_params[, n_stage]
  lefkovitch_matrix[, param_list$fec_stages, 1] <- fec_params
  
  # initialise abundances
  init <- lognormal(meanlog = param_list$init_mean, sdlog = param_list$init_sd,
                    dim = c(n_stage, n_site))
  
  # create a parameter to capture density dependence
  density <- uniform(min = param_list$dens_lower,
                     max = param_list$dens_upper,
                     dim = n_site)
  
  # collate and return outputs  
  list(matrix = lefkovitch_matrix,
       init = init,
       dens_param = density,
       coefs_surv = coefs_surv,
       coefs_growth = coefs_growth,
       coefs_fec = coefs_fec,
       surv_params = surv_params,
       fec_params = fec_params,
       growth_params = growth_params,
       gamma_year_surv = gamma_year_surv,
       gamma_year_growth = gamma_year_growth,
       gamma_year_fec = gamma_year_fec)
  
}

# switch out system codes for names
system_switch_fun <- function(x) {
  
  if (!is.character(x)) {
    x <- as.character(x)
  }
  
  out <- rep(NA, length(x))
  for (i in seq_along(x)) {
    out[i] <- switch(substr(x[i], 1, 2),
                     'GO' = 'GOULBURN',
                     'BR' = 'BROKEN',
                     'PC' = 'PYRAMID CK',
                     'LO' = 'LODDON',
                     'CA' = 'CAMPASPE')
  }  
  
  out
  
}

# hist function to handle missing data and return counts only
hist_fn <- function(x, breaks) {
  out <- rep(0, length(breaks) - 1)
  if (length(x) & !all(is.na(x))) {
    out <- hist(x, plot = FALSE, breaks = breaks)$counts
  }
  out
}
