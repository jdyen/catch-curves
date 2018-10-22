# greta scripts for catch curve analysis

# create separable covariance function from two component covariances
create_greta_covar <- function(n_sp, n_class,
                               eta_sp = 1, eta_class = 4,
                               sigma_mean = 0.0, sigma_sd = 1.0) {
  
  if (length(sigma_mean) == 1)
    sigma_mean <- rep(sigma_mean, 2)
  if (length(sigma_sd) == 1)
    sigma_sd <- rep(sigma_sd, 2)
  
  species_corr_est <- lkj_correlation(eta = eta_sp, dim = n_sp)
  class_corr_est <- lkj_correlation(eta = eta_class, dim = n_class)
  species_sigma <- lognormal(mean = sigma_mean[1], sd = sigma_sd[1], dim = n_sp)
  class_sigma <- lognormal(mean = sigma_mean[2], sd = sigma_sd[2], dim = n_class)
  species_covar_est <- (species_sigma %*% t(species_sigma)) * species_corr_est
  class_covar_est <- class_sigma %*% t(class_sigma) * class_corr_est
  covar_est <- kronecker(species_covar_est, class_covar_est)
  
  covar_est
  
}

# set priors for covariate model
prepare_greta_model <- function(data, flow_data = NULL, n_pc = 3, ...) {
  
  # extract data
  stage_data <- data$age
  if (is.null(flow_data)) {
    flow_data <- data$flow
  } else {
    flow_data <- flow_data$scores[, seq_len(n_pc)]
  }
  info_data <- data$info
  n_sp <- data$n_sp
  n_class <- data$n_class
  n_obs <- data$n_obs
  
  # calculate indices used repeatedly
  n_class_sp <- (n_sp * n_class)
  sp <- rep(seq_len(n_sp), each = n_class)
  class <- rep(seq_len(n_class), times = n_sp)
  n_preds <- ncol(flow_data)
  
  # create a separable covariance matrix with spp and age components
  covar_mat <- create_greta_covar(n_sp = n_sp,
                                  n_class = n_class)
  
  # define mean vector
  alpha <- normal(mean = 0.0, sd = 1.0, dim = 1)
  alpha_sp <- normal(mean = 0.0, sd = 1.0, dim = n_sp)
  alpha_class <- normal(mean = 0.0, sd = 1.0, dim = n_class)
  alpha_sp_by_class <- normal(mean = 0.0, sd = 1.0, dim = n_class_sp)
  Alpha <- alpha + alpha_sp[sp] + alpha_class[class] + alpha_sp_by_class

  # define coefficient matrix  
  beta <- normal(mean = 0.0, sd = 1.0, dim = n_preds)
  beta_sp <- normal(mean = 0.0, sd = 1.0, dim = c(n_preds, n_sp))
  beta_class <- normal(mean = 0.0, sd = 1.0, dim = c(n_preds, n_class))
  beta_sp_by_class <- normal(mean = 0.0, sd = 1.0, dim = c(n_preds, n_class_sp))
  Beta <- zeros(n_preds, n_class_sp)
  for (i in seq_len(ncol(Beta)))
    Beta[, i] <- beta + beta_sp[, sp[i]] + beta_class[, class[i]] + beta_sp_by_class[, i]

  # define exchangeable terms (pseudo random effects)
  river <- as.integer(as.factor(info_data$system))
  year <- as.integer(as.factor(info_data$year))
  river_by_year <- as.integer(as.factor(apply(cbind(info_data$system, info_data$year), 1, paste, collapse = "_")))
  river_by_reach <- as.integer(as.factor(apply(cbind(info_data$system, info_data$reach), 1, paste, collapse = "_")))
  nriver <- length(unique(river))
  nyear <- length(unique(year))
  nriver_by_year <- length(unique(river_by_year))
  nriver_by_reach <- length(unique(river_by_reach))
  gamma_river <- normal(mean = 0.0, sd = 1.0, dim = nriver)
  gamma_year <- normal(mean = 0.0, sd = 1.0, dim = nyear)
  gamma_river_by_year <- normal(mean = 0.0, sd = 1.0, dim = (nriver_by_year))
  gamma_river_by_reach <- normal(mean = 0.0, sd = 1.0, dim = (nriver_by_reach))
  gamma_river_by_sp <- normal(mean = 0.0, sd = 1.0, dim = c(nriver, n_sp))
  gamma_year_by_sp <- normal(mean = 0.0, sd = 1.0, dim = c(nyear, n_sp))
  gamma_river_by_class <- normal(mean = 0.0, sd = 1.0, dim = c(nriver, n_class))
  gamma_year_by_class <- normal(mean = 0.0, sd = 1.0, dim = c(nyear, n_class))
  gamma_river_by_sp_by_class <- normal(mean = 0.0, sd = 1.0, dim = c(nriver, n_class_sp))
  gamma_year_by_sp_by_class <- normal(mean = 0.0, sd = 1.0, dim = c(nyear, n_class_sp))
  Gamma <- zeros(n_obs, n_class_sp)
  for (i in seq_len(ncol(Gamma))) {
    Gamma[, i] <- gamma_river_by_sp[river, sp[i]] +
      gamma_year_by_sp[year, sp[i]] +
      gamma_river_by_class[river, class[i]] +
      gamma_year_by_class[river, class[i]] +
      gamma_river_by_sp_by_class[river, i] +
      gamma_year_by_sp_by_class[year, i]
  }
  Gamma <- sweep(Gamma, 1,
                 gamma_river[river] + gamma_year[year] +
                   gamma_river_by_year[river_by_year] +
                   gamma_river_by_reach[river_by_reach],
                 "+")
    
  # calculate linear predictor from its components
  mu <- sweep((flow_data %*% Beta + Gamma), 2, Alpha, "+")

  # add MVN errors to IID mu values
  resid <- multivariate_normal(mean = rep(0.0, n_class_sp), Sigma = covar_mat, dim = n_obs)
  mu_mvn <- mu + resid
  ## SHOULD WORK (but is probably slower?):
  # mvn_helper <- normal(mean = 0.0, sd = 1.0, dim = dim(mu))
  # mu_mvn <- mu + mvn_helper %*% chol(covar_mat)
     
  # set likelihood
  distribution(stage_data) <- poisson(lambda = exp(mu_mvn))
  
  # compile model
  mod <- model(Alpha, Beta, Gamma, mu, mu_mvn, covar_mat)
  
  mod
  
}
