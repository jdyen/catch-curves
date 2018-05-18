# helpers for catch curve analysis
remove_correlated <- function(flow) {
  
  # calculate correlations > 0.7 and remove one at a time, taking left-most of the columns
  #   correlated with the most other variables
  cor_vals <- cor(flow, use = "complete")
  sum_correlated <- apply(cor_vals, 2, function(x) sum(abs(x) > 0.7) - 1)
  while(any(sum_correlated > 0)) {
    cor_vals <- cor_vals[-max(which(sum_correlated == max(sum_correlated))), -max(which(sum_correlated == max(sum_correlated)))]
    sum_correlated <- apply(cor_vals, 2, function(x) sum(abs(x) > 0.7) - 1)
  }
  flow[, match(names(sum_correlated), colnames(flow))]
  
}

na_replace_fun <- function(x) {
  
  if (any(is.na(x))) 
    x[is.na(x)] <- mean(x, na.rm = TRUE)

  x
  
}

calc_flow_pc <- function(flow, scale = FALSE) {
  
  flow <- apply(flow, 2, na_replace_fun)
  if (scale)
    flow <- scale(flow)
  pc_out <- princomp(flow)
  
  pc_out
  
}

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