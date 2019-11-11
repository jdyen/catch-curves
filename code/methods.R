# helper methods for a CCR model

# get a parameter from a fitted MCMC object
get_param <- function(samples, regex_par)
  samples[, grep(regex_par, colnames(samples))]

# predict to newdata from a fitted MCMC object 
predict.ccr_model <- function(
  obj,
  newdata = NULL,
  thin = 1,
  lengths = TRUE,
  survey = FALSE,
  effort = NULL
) {
  
  # are newdata provided?
  if (is.null(newdata))
    newdata <- obj$data

  # turn off survey random effects if not included in fitted model
  if (survey)
    survey <- ifelse(obj$include$survey, TRUE, FALSE)
  
  # what if effort data are not provided?
  if (is.null(newdata$effort))
    newdata$effort <- rep(1, times = length(newdata$year))

  # create data inputs from provided data
  data <- prepare_data(newdata$length_age_matrix,
                       newdata$system,
                       newdata$predictors,
                       newdata$effort)

  # match up surveys if needed
  if (survey) {
    survey_ids <- data.frame(survey = seq_along(obj$data$system),
                             system = obj$data$system,
                             year = obj$data$year)
    idx <- paste(survey_ids$system, survey_ids$year, sep = "_")
    idy <- paste(data$sys_vec, rep(newdata$year, ncol(newdata$length_age_matrix)), sep = "_")
    data$survey_vec <- survey_ids$survey[match(idy, idx)]
    if (any(is.na(data$survey_vec)))
      data$survey_vec[is.na(data$survey_vec)] <- max(data$survey_vec, na.rm = TRUE) + 1
    
  }
  
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
  new_predictors <- newdata$predictors
  new_age_vec <- data$age_vec
  survey_index <- data$survey_vec
  new_survey <- data$survey_vec
  new_effort <- data$effort_vec
  
  # are predictors formatted correctly?
  if (nrow(new_predictors) != length(new_survey))
    new_predictors <- new_predictors[new_survey, ]
  
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
  survey_est <- get_param(samples, "gamma_survey")
  pred_est <- get_param(samples, "pred_effects")
  
  # we need some indices to reformat arrays
  n_age <- as.numeric(strsplit(strsplit(colnames(convert)[ncol(convert)], ",")[[1]][2], "\\]")[[1]])
  n_len <- ncol(convert) / n_age
  n_system <- ncol(alpha_est)
  n_survey <- ncol(survey_est)
  n_sys_age <- n_age * n_system
  if (obj$include$sys_flow)
    n_predictors <- ncol(pred_est) / n_system
  else
    n_predictors <- ncol(pred_est)
    
  # reformat into correct dimensions
  if (obj$include$sys_flow)
    pred_array <- array(pred_est, dim = c(n_samples, n_system, n_predictors))
  else 
    pred_array <- array(pred_est, dim = c(n_samples, n_predictors))
  convert_array <- array(convert, dim = c(n_samples, n_len, n_age))
  
  # need to check that new systems, cohorts, and years are within previous bounds
  new_sys[new_sys > n_system] <- n_system + 1
  new_survey[new_survey > n_survey] <- n_survey + 1
  new_sys_age[new_sys_age > n_sys_age] <- n_sys_age + 1
  
  # add zeros to set up marginal predictions for new levels
  alpha_est <- cbind(alpha_est, rep(0, n_samples))
  beta_est <- cbind(beta_est, rep(0, n_samples))
  survey_est <- cbind(survey_est, rep(0, n_samples))
  if (obj$include$sys_flow)
    pred_array <- abind::abind(pred_array, along = 2, array(0, dim = c(n_samples, n_predictors)))
  
  # calculate linear predictor
  mu_pred <- alpha_est[, new_sys] +
    sweep(beta_est[, new_sys], 2, new_age_vec, "*") +
    log(new_effort)

  # is flow included by system?
  if (obj$include$sys_flow)
    mu_pred <- mu_pred + apply(sweep(pred_array[, new_sys, ], c(2, 3), new_predictors, "*"), c(1, 2), sum)
  else
    mu_pred <- mu_pred + pred_array %*% t(new_predictors)
  
  # should we add random effects?
  if (survey)
    mu_pred <- mu_pred + survey_est[, new_survey]

  # put back onto the observation scale
  ages_pred <- mu_pred

  # reformat predictions
  unique_ages <- unique(new_age_vec)
  ages_formatted <- array(0, dim = c(n_samples, n_pred, n_age))
  for (i in seq_along(unique_ages)) {
    idx <- new_age_vec == unique_ages[i]
    ages_formatted[, survey_index[idx], unique_ages[i]] <- ages_pred[, idx]
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
create_newdata <- function(obj, 
                           var,
                           nplot = 100,
                           system = 1,
                           year = 1,
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

  # make up newdata objects for each variable
  sys_vec <- rep(rep(system, each = nplot), each = nyear)
  year_vec <- rep(rep(year, each = nplot), times = nsystem)
  predictors <- matrix(0, nrow = length(year_vec), ncol = ncol(obj$data$predictors))
  colnames(predictors) <- colnames(obj$data$predictors)
  for (i in seq_len(nsystem)) {

    # which rows do we want to change in the new data?    
    idx <- sys_vec == i
    
    # what about the old data?
    idy <- rep(obj$data$system, times = ncol(obj$data$length_age_matrix)) == i
    
    # what range do the variables have in system i?
    var_range <- range(obj$data$predictors[idy, var])
    
    # create a sequence spanning that range
    predictors[idx, var] <- rep(seq(var_range[1], var_range[2], length = nplot), times = nyear)
    
    # fill the other variables with their system-level mean
    excluded <- !colnames(predictors) %in% var
    predictors[idx, excluded] <- matrix(
      rep(
        apply(obj$data$predictors[idy, excluded], 2, mean),
        each = sum(idx)
      ),
      ncol = ncol(predictors) - 1)
  }
  effort <- rep(effort, times = nrow(predictors))

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
