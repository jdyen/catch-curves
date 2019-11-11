# set validate as an S3 method
validate <- function(obj, folds, ...) {
  UseMethod("validate")
}

# run k-fold cross validation on a fitted CCR model
validate.ccr_model <- function(obj, folds = 10, ...) {

  # how many observations?
  n <- nrow(obj$data$response)
  
  # pull the settings from the fitted model
  priors <- obj$priors
  mcmc_settings <- obj$mcmc_settings

  # create folds
  folds <- create_folds(n, folds)
  
  # create test and train data sets from folds
  data_folds <- lapply(folds, extract_folds, obj$data)

  # fit a model for each fold
  mod_cv <- lapply(
    data_folds,
    validate_internal,
    priors,
    mcmc_settings,
    ...
  )

  # unpack list
  out <- list()
  out$predicted <- do.call(rbind, lapply(mod_cv, function(x) x$predicted))
  out$observed <- do.call(rbind, lapply(mod_cv, function(x) x$observed))
  
  # return predicted and observed
  out
  
}

# internal function to run cross validation
validate_internal <- function(data, priors, mcmc_settings, ...) {
  
  # unpack data
  response <- data$data_train$response
  length_age_matrix <- data$data_train$length_age_matrix
  predictors <- data$data_train$predictors
  effort <- data$data_train$effort
  system <- data$data_train$system
  year <- data$data_train$year
  n_age <- ncol(response)
  sys_year <- data.frame(system = rep(system, n_age),
                         year = c(sapply(seq_len(n_age), function(x) year - x + 1)))

  # fit a model
  mod_tmp <- fit_ccr(response = response,
                     length_age_matrix = length_age_matrix,
                     predictors = predictors,
                     effort = effort,
                     system = system,
                     year = year,
                     sys_year_flow = sys_year,
                     include = list(sys_flow = TRUE,
                                    survey = TRUE),
                     priors = priors,
                     mcmc_settings = mcmc_settings)
  
  # predict to holdout data
  predictions <- predict(mod_tmp, data$data_test, ...)
  
  # reduce predictions to a mean value
  predictions <- apply(predictions, c(2, 3), mean)
  
  # return predictions and observed values
  list(predicted = predictions, observed = data$data_test$response)
  
}

# check folds for cross validation and create if needed
create_folds <- function(n, folds) {
  
  # have folds been pre-prepared?
  already_list <- is.list(folds)
  out <- folds
  
  # if not, create them from a numeric value
  if (!already_list) {
    
    # is folds numeric?
    if (!is.numeric(folds) & !is.integer(folds))
      stop("folds must be provided as a list or scalar value", call. = FALSE)
    
    # how big is each fold?
    cv_size <- floor(n / folds)
    
    # create a vector random indices to sample from (shuffling the data prior to CV)
    idx <- sample(seq_len(n), size = n, replace = FALSE)
    
    # fill all
    out <- list()
    for (i in seq_len(folds))
      out[[i]] <- idx[((i - 1) * cv_size + 1):(i * cv_size)]
    
    # replace final fold with one that includes the final obs
    out[[folds]] <- idx[((folds - 1) * cv_size + 1):n]
    
  }
  
  # return outputs
  out
  
}

# construct data lists with correct test/train splits
extract_folds <- function(folds, data) {
  
  # pull out training data
  data_train <- list(response = data$response[-folds, ],
                     length_age_matrix = data$length_age_matrix,
                     system = data$system[-folds],
                     year = data$year[-folds],
                     predictors = data$predictors[-folds, ],
                     effort = data$effort[-folds],
                     cohort_mat = data$cohort_mat[-folds, ])
  
  # pull out testing data
  data_test <- list(response = data$response[folds, ],
                    length_age_matrix = data$length_age_matrix,
                    system = data$system[folds],
                    year = data$year[folds],
                    predictors = data$predictors[folds, ],
                    effort = data$effort[folds],
                    cohort_mat = data$cohort_mat[folds, ])
  
  # return outputs
  list(data_train = data_train, data_test = data_test)
  
}
