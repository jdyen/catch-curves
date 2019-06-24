validate_glmer <- function(obj, folds, settings = list()) {
  
  # unpack settings
  sets <- list(iter = obj$call$iter,
               chains = obj$call$chains)
  sets[names(settings)] <- settings
  
  # how many obs?
  all_vars <- get_all_vars(obj$formula, obj$data)
  n_obs <- nrow(all_vars)
  
  # define folds
  if (is.numeric(folds))
    folds <- define_cv_folds(folds, n_obs)

  # run cv function
  cv_vals <- lapply(folds, cv_fun, obj, all_vars, sets)
  
  # return validation metrics
  list(r2_naive = cor(obj$fitted.values, all_vars[, 1]),
       r2_cv = cor(do.call(c, sapply(cv_vals, function(x) x$fitted)),
                   do.call(c, sapply(cv_vals, function(x) x$observed))))
  
}

cv_fun <- function(idx, obj, all_vars, settings) {
  
  # subset data
  data_train <- all_vars[-idx, ]
  data_test <- all_vars[idx, ]

  # unpack settings
  iter <- settings$iter
  chains <- settings$chains
  
  # fit model
  out <- rstanarm::stan_glmer(obj$formula, data = data_train,
                              iter = iter, chains = chains, cores = 1)

  # predict
  out <- posterior_predict(out, newdata = data_test)
  
  # return mean predictions
  list(fitted = apply(out, 2, mean), observed = data_test[, 1])
  
}

define_cv_folds <- function(folds, n_obs) {
  
  # how big is each fold?
  cv_size <- floor(n_obs / folds)
  
  # create a vector random indices to sample from (shuffling the data prior to CV)
  idx <- sample(seq_len(n_obs), size = n_obs, replace = FALSE)
  
  # fill all
  out <- list()
  for (i in seq_len(folds))
    out[[i]] <- idx[((i - 1) * cv_size + 1):(i * cv_size)]
  
  # replace final fold with one that includes the final obs
  out[[folds]] <- idx[((folds - 1) * cv_size + 1):n_obs]
  
  # return outputs
  out

}
