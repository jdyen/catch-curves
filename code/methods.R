# helper methods for a CCR model

# get a parameter from a fitted MCMC object
get_param <- function(samples, regex_par)
  samples[, grep(regex_par, colnames(samples))]

# predict to newdata from a fitted MCMC object 
predict.ccr_model <- function(obj, newdata = NULL, thin = 1, lengths = TRUE,
                              year = FALSE, cohort = FALSE,
                              effort = NULL) {
  
  # are newdata provided?
  if (is.null(newdata))
    newdata <- obj$data
  
  # what if effort data are not provided?
  if (is.null(newdata$effort))
    newdata$effort <- rep(1, nrow(newdata$predictors))
  
  # create data inputs from provided data
  data <- prepare_data(newdata$length_age_matrix,
                       newdata$system,
                       newdata$year,
                       newdata$predictors,
                       newdata$effort)
  
  # do we want to standardise efforts?
  if (!is.null(effort)) {
    if (length(effort) == 1)
      effort <- rep(effort, length(data$sys_vec))
    data$effort_vec <- effort
  }
  
  # extract MCMC samples from fitted model
  samples <- do.call(rbind, obj$draws)
  
  # unpack newdata
  new_sys <- data$sys_vec
  new_cohort <- c(newdata$cohort_mat)
  new_year <- data$year_vec
  new_predictors <- newdata$predictors
  new_age_vec <- data$age_vec
  new_survey <- data$survey_vec
  new_effort <- data$effort_vec
  
  # calculate sys_age combo from system and age
  n_age <- ncol(newdata$length_age_matrix)
  new_sys_age <- n_age * (new_sys - 1) + new_age_vec
  
  # combine the draws and thin if needed
  if (thin > 1)
    samples <- samples[seq(1, nrow(samples), by = thin), ]
  
  # how many samples and predictions do we have?
  n_samples <- nrow(samples)
  n_pred <- length(newdata$system)
  
  # pull out params we need
  convert <- get_param(samples, "length_to_age")
  alpha_est <- get_param(samples, "alpha\\[")
  beta_est <- get_param(samples, "beta\\[")
  cohort_est <- get_param(samples, "gamma_cohort")
  year_est <- get_param(samples, "gamma_year")
  pred_est <- get_param(samples, "pred_effects")
  
  # we need some indices to reformat arrays
  n_age <- as.numeric(strsplit(strsplit(colnames(convert)[ncol(convert)], ",")[[1]][2], "\\]")[[1]])
  n_len <- ncol(convert) / n_age
  n_system <- ncol(alpha_est)
  n_cohort <- ncol(cohort_est)
  n_year <- ncol(year_est)
  n_sys_age <- n_age * n_system
  n_predictors <- ncol(pred_est) / n_sys_age
  
  # reformat into correct dimensions
  pred_array <- array(pred_est, dim = c(n_samples, n_sys_age, n_predictors))
  convert_array <- array(convert, dim = c(n_samples, n_len, n_age))
  
  # need to check that new systems, cohorts, and years are within previous bounds
  new_sys[new_sys > n_system] <- n_system + 1
  new_cohort[new_cohort > n_cohort] <- n_cohort + 1
  new_year[new_year > n_year] <- n_year + 1
  new_sys_age[new_sys_age > n_sys_age] <- n_sys_age + 1
  
  # add zeros to set up marginal predictions for new levels
  alpha_est <- cbind(alpha_est, rep(0, n_samples))
  beta_est <- cbind(beta_est, rep(0, n_samples))
  cohort_est <- cbind(cohort_est, rep(0, n_samples))
  year_est <- cbind(year_est, rep(0, n_samples))
  pred_array <- abind::abind(pred_array, along = 2, array(0, dim = c(n_samples, n_predictors)))
  
  # calculate linear predictor
  mu_pred <- alpha_est[, new_sys] +
    sweep(beta_est[, new_sys], 2, new_age_vec, "*") +
    apply(sweep(pred_array[, new_sys_age, ], c(2, 3), new_predictors[new_survey, ], "*"), c(1, 2), sum) +
    log(new_effort)

  # should we add random effects?
  if (year)
    mu_pred <- mu_pred + year_est[, new_year]
  if (cohort)
    mu_pred <- mu_pred + cohort_est[, new_cohort] 
    
  # put back onto the observation scale
  ages_pred <- mu_pred
  # ages_pred <- exp(mu_pred)
  
  # reformat predictions
  unique_ages <- unique(new_age_vec)
  ages_formatted <- array(0, dim = c(n_samples, n_pred, n_age))
  for (i in seq_along(unique_ages)) {
    idx <- new_age_vec == unique_ages[i]
    ages_formatted[, new_survey[idx], unique_ages[i]] <- ages_pred[, idx]
  }
  
  # convert to length (would be good to vectorise this but lapply options are slower)
  if (lengths) {
    pred_out <- array(NA, dim = c(n_samples, n_pred, n_len))
    for (i in seq_len(n_samples))
      pred_out[i, , ] <- ages_formatted[i, , ] %*% t(convert_array[i, , ])
  } else {
    pred_out <- ages_formatted
  }
  
  # return all predictions
  pred_out
  
}

# create a function to define a set of predictor combos from some simple settings
create_newdata <- function(obj, var,
                           nplot = 100,
                           system = 1, year = 1,
                           effort = 1) {
  
  # is a model object provided?
  if (missing(obj))
    stop("obj must be a ccr_model object", call. = FALSE)
  
  # is a model object provided?
  if (!"ccr_model" %in% class(obj))
    stop("obj must be a ccr_model object", call. = FALSE)
  
  # are one or more variables listed?
  if (missing(var))
    stop("var must include a character name of one variable", call. = FALSE)

  # are one or more variables listed?
  if (length(var) > 1)
    stop("var must include one variable only", call. = FALSE)
  
  # how many observations do we need to predict?
  nsystem <- length(unique(system))
  nyear <- length(unique(year))

  # what range do the variables have?
  var_range <- range(obj$data$predictors[, var])
  
  # make up newdata objects for each variable
  sys_vec <- rep(rep(system, each = nplot), each = nyear)
  year_vec <- rep(rep(year, each = nplot), times = nsystem)
  predictors <- matrix(0, nrow = length(year_vec), ncol = ncol(obj$data$predictors))
  colnames(predictors) <- colnames(obj$data$predictors)
  predictors[, var] <- rep(seq(var_range[1], var_range[2], length = nplot), times = nyear * nsystem)
  effort <- rep(effort, times = nrow(predictors))
  
  # what if we are changing a variable that's also included quadratically?
  if (length(grep("prop_sum", var))) {
    predictors[, "prop_sum_win_sq"] <- rep(seq(min(obj$data$predictors[, "prop_sum_win_sq"]),
                                               max(obj$data$predictors[, "prop_sum_win_sq"]),
                                               length = nplot),
                                           times = nyear * nsystem)
    predictors[, "prop_sum_lt_win"] <- rep(seq(min(obj$data$predictors[, "prop_sum_lt_win"]),
                                               max(obj$data$predictors[, "prop_sum_lt_win"]),
                                               length = nplot),
                                           times = nyear * nsystem)
  }

  # return the newdata object
  list(length_age_matrix = obj$data$length_age_matrix,
       system = sys_vec,
       year = year_vec,
       predictors = predictors,
       effort = effort)
  
}

# calculate fit metrics for a fitted CCR model
calculate_metrics <- function(obj) {

  # pull out observed data
  observed <- obj$data$response
  
  # calculate fitted values
  fitted <- predict(obj, random = TRUE)

  # pull out a few different r2 values
  dev_null <- NULL
  dev_fitted <- NULL
  
  
}

# calculate deviance from a fitted model
calc_deviance <- function(x, y) {
  
  loglik <- NULL
  
  -2 * loglik
  
}
