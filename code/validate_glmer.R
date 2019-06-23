validate_glmer <- function(obj, folds, settings = list(), future = future::sequential()) {
  
  # unpack settings
  sets <- list()  ## TAKE FROM obj
  sets[names(settings)] <- settings
  
  # how many obs?
  n_obs <- length(obj$data$y)
  
  # define folds
  if (is.numeric(folds))
    folds <- define_cv_folds(folds, n_obs)
  
  # set eval type
  future::future(future)
  
  # run function
  future.apply::future_sapply(folds, cv_fun, obj, settings)
  
}

cv_fun <- function(idx, obj, settings) {
  
  # subset data
  data_train
  data_test

  # unpack settings
  
  
  # fit model
  out <- rstanarm::stan_glmer(form)

  # predict
  out <- posterior_predict(out, newdata = data_test)
  
  # pull out mean predictions if needed
  #out <- apply(out, 2, mean)
  
  out
    
}

define_cv_folds <- function(folds, n_obs) {
  
  # how big is each fold?
  cv_size <- floor(n_obs / folds)
  
  # fill all
  out <- list()
  for (i in seq_len(folds))
    out[[i]] <- ((i - 1) * cv_size + 1):(i * cv_size)
  
  # replace final fold with one that includes the final obs
  out[[folds]] <- ((folds - 1) * cv_size + 1):n_obs
  
  # return outputs
  out

}
