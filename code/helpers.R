# histogram function to return counts only
hist_fn <- function(x, breaks) {
  hist(x, breaks = breaks, plot = FALSE)$counts
}

# create a pairwise matrix matching classes from two variables
classify <- function(x, y) {
  out <- matrix(0, nrow = max(x, na.rm = TRUE), ncol = max(y, na.rm = TRUE))
  for (i in seq_along(x)) {
    if (!is.na(x[i]) & !is.na(y[i]))
      out[x[i], y[i]] <- out[x[i], y[i]] + 1
  }
  out
}

# define a function to convert size to age based on model in Todd & Koehn 2008
inverse_growth <- function(x, length_inf, time_zero, k_param, c_param) {
  
  par1 <- 1 / k_param
  ratio_param <- 1 - (time_zero / length_inf)
  par2 <- 1 / (1 - c_param * ratio_param)
  par3 <- (length_inf - time_zero) / (length_inf - x)
  par4 <- c_param * ratio_param
  
  par1 * log(par2 * (par3 - par4))
  
}

# function from Todd & Koehn 2008 to estimate length from age
length_from_age <- function(x, length_inf, time_zero, k_param, c_param) {
  
  par1 <- length_inf
  par2 <- par1 - time_zero
  ratio_param <- 1 - (time_zero / length_inf)
  exp_param <- exp(k_param * x)
  par3 <- exp_param + (1 - exp_param) * c_param * ratio_param
  
  par1 - (par2 / par3)
  
}

# rebase a vector intended as an index to start at 1 and count up with no gaps
rebase_index <- function(x)
  as.integer(as.factor(x))

# rescale data and store the standardisations
extract_standards <- function(x) {
  c("mean" = mean(x, na.rm = TRUE), "sd" = sd(x, na.rm = TRUE))
}

# subset a data object from a fitted CCR model
subset_data <- function(x, idx) {
  
  list(response = x$response[idx, ],
       length_age_matrix = x$length_age_matrix,
       system = x$system[idx],
       year = x$year[idx],
       predictors = x$predictors[idx, ],
       cohort_mat = x$cohort_mat[idx, ])
  
}

# helper function to match cols and fill empties
col_match <- function(x, names, classes = NULL) {
  
  idx <- match(names, colnames(x))
  if (any(is.na(idx)))
    idx[is.na(idx)] <- ncol(x) + 1
  x <- cbind(x, rep(NA, nrow(x)))
  
  out <- x[, idx]
  colnames(out) <- names
  
  out  
  
}

add_scientific_names <- function(x) {
  
  sci_names <- c(maccullochellapeelii = "Maccullochella peelii",
                 maccullochellamacquariensis = "Maccullochella macquariensis",
                 cyprinuscarpio = "Cyprinus carpio",
                 macquariaambigua = "Macquaria ambigua",
                 macquariaaustralasica = "Macquaria australasica",
                 melanotaeniafluviatilis = "Melanotaenia fluviatilis",
                 bidyanusbidyanus = "Bidyanus bidyanus")
  
  sci_names[x]
  
}

add_common_names <- function(x) {
  
  common_names <- c(maccullochellapeelii = "Murray cod",
                    maccullochellamacquariensis = "Trout cod",
                    cyprinuscarpio = "Common carp",
                    macquariaambigua = "Golden perch",
                    macquariaaustralasica = "Macquarie perch",
                    melanotaeniafluviatilis = "Murray river rainbowfish",
                    bidyanusbidyanus = "Silver perch") 
  
  common_names[x]
  
}

# filter a matrix down by column to remove correlated columns
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

# define some helper functions
max_fun <- function(x) {
  
  out <- NA
  
  if (any(!is.na(x)))
    out <- max(x, na.rm = TRUE)
  
  out
  
}

na_replace_fun <- function(x) {
  
  if (any(is.na(x))) 
    x[is.na(x)] <- mean(x, na.rm = TRUE)

  x
  
}

calculate_length_conversions <- function(data) {
  
  # pull out species names
  sp_names <- unique(data$Common.Name)
  
  # set up output list
  length_weight_conversion <- list()
  
  # loop through each species
  for (i in seq_along(sp_names)) {
    
    # subset for species i
    dat <- data[which(data$Common.Name == sp_names[i]), ]
    
    # log transform length and weight
    rows_to_rm <- NULL
    if (any(is.na(dat$totallength)))
      rows_to_rm <- c(rows_to_rm, which(is.na(dat$totallength)))
    if (any(is.na(dat$WEIGHT)))
      rows_to_rm <- c(rows_to_rm, which(is.na(dat$WEIGHT)))
    if (any(dat$totallength <= 0, na.rm = TRUE))
      rows_to_rm <- c(rows_to_rm, which(dat$totallength <= 0))
    if (any(dat$WEIGHT <= 0, na.rm = TRUE))
      rows_to_rm <- c(rows_to_rm, which(dat$WEIGHT <= 0))
    if (length(rows_to_rm))
      dat <- dat[-rows_to_rm, ]
    log_length <- log(dat$totallength)
    log_weight <- log(dat$WEIGHT)
    
    # fit linear model to log-log transformed data
    if (length(log_weight) > 10)
      mod_tmp <- lm(log_weight ~ log_length)
    
    # pull out coefs and other stats
    if (length(log_weight) > 10) {
      coefs_tmp <- coef(mod_tmp)
      n_obs <- length(log_length)
      r2_vals <- summary(mod_tmp)$r.squared
      resid_tmp <- residuals(mod_tmp)
    } else {
      coefs_tmp <- NA
      n_obs <- 0
      r2_vals <- NA
      resid_tmp <- NA
    }
    
    # save outputs
    length_weight_conversion[[sp_names[i]]] <-
      list(n = n_obs,
           coef = coefs_tmp,
           length = dat$totallength,
           weight = dat$WEIGHT,
           r2 = r2_vals,
           resid = resid_tmp)
    
  }
  
  coefs_all <- lapply(length_weight_conversion, function(x) x$coef)
  coefs_all <- coefs_all[-which(is.na(coefs_all))]
  
  length_weight_conversion$generic <- apply(matrix(unlist(coefs_all), ncol = 2, byrow = TRUE),
                                            2, mean)
  
  length_weight_conversion
  
}

impute_weights <- function(data, length_conversions) {
  
  # fill missing weights
  sp_tmp <- unique(data$Common.Name)
  for (i in seq_along(sp_tmp)) {
    
    # subset to single species
    dat_tmp <- data[which(data$Common.Name == sp_tmp[i]), ]
    
    # check to make sure some weights are missing      
    if (any(is.na(dat_tmp$WEIGHT))) {
      
      # use species-specific equation if it exists, generic otherwise
      if (length_conversions[[sp_tmp[i]]]$n) {
        coefs <- length_conversions[[sp_tmp[i]]]$coef
      } else {
        coefs <- length_conversions$generic
      }
      
      # subset NA observations
      na_sub <- which(is.na(dat_tmp$WEIGHT))
      
      # estimate weight from length
      dat_tmp$WEIGHT[na_sub] <- exp(coefs[1] + coefs[2] * log(dat_tmp$totallength[na_sub]))
      
      # return estimated data to full data set
      data[which(data$Common.Name == sp_tmp[i]), ] <- dat_tmp
      
    }
  }
  
  data
  
}

clean_reaches <- function(data) {
  
  # arrange alldat to have more consistent set of variables
  data <- data.frame(date = data$Date,
                     year = data$YEAR,
                     site = data$SITE_CODE,
                     system = data$SYSTEM,
                     reach = data$Reach,
                     species = data$Common.Name,
                     length = data$totallength,
                     weight = data$WEIGHT,
                     abundance = data$Total.Sampled,
                     intensity = (data$total_no_passes * data$seconds))
  
  data$reach_alt <- data$reach
  data$reach_alt <- ifelse(data$system == "BROKEN",
                           ifelse(data$reach_alt == 4, 5, data$reach_alt), # downstream 
                           data$reach_alt)
  data$reach_alt <- ifelse(data$system == "THOMSON",
                           ifelse(data$reach_alt == 2, 3,             # downstream 
                                  ifelse(data$reach_alt == 6, 5,      # upstream
                                         data$reach_alt)),
                           data$reach_alt)
  data$reach_alt <- ifelse(data$system == "LODDON",
                           ifelse(data$reach_alt == 2, 3,             # downstream
                                  ifelse(data$reach_alt == 5, 4,      # upstream 
                                         data$reach_alt)),
                           data$reach_alt)
  
  data
  
}

# internal functon for rolling_range (below)
get_range <- function(x) {
  
  out <- NA
  
  if (any(!is.na(x)))
    out <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    
  out
  
}

# alternative internal functon for rolling_range (below)
get_range2 <- function(x) {
  
  out <- NA
  
  if (any(!is.na(x))) {
    out <- 0
    if (min(x, na.rm = TRUE) > 0)
      out <- max(x, na.rm = TRUE) / min(x, na.rm = TRUE)
  }
  
  out
  
}

# calculate maximum rolling range in a variable over a set lag or lead time
rolling_range <- function(data, lag, variable = NULL){ 
  
  if (length(dim(data)) < 2) {
    data <- as.matrix(data, ncol = 1)
    variable <- 1
  }
  
  if (ncol(data) == 1)
    variable <- 1
  
  if (is.null(variable))
    stop("`variable` must be provided if `data` has more than one column", call. = FALSE)
  
  nrows <- nrow(data)
  
  idx <- sapply(rev(seq_len(lag)) - 1, function(x) rep(x, nrows))
  idx <- sweep(idx, 1, seq_len(nrows), "+")
  idx <- ifelse(idx > nrows, NA, idx)
  
  df <- matrix(data[c(idx), variable], nrow = nrows)
  
  diff <- apply(df, 1, get_range2)
  
  max(diff, na.rm = TRUE)
  
}   

# function to rename river systems from sites
system_switch_fun <- function(x) {
  
  x <- as.character(x)
  
  out <- rep(NA, length(x))
  for (i in seq_along(x)) {
    
    out[i] <- switch(substr(x[i], 1, 2),
                     "LO" = "LODDON",
                     "GO" = "GOULBURN",
                     "BR" = "BROKEN",
                     "PC" = "PYRAMIDCK",
                     "CA" = "CAMPASPE")
    
  }
  
  out
  
}

# set month definitions
define_months <- function() {
  list(
    all_months = c(1:12),
    summer_months = c(12, 1:3),
    cooler_months = c(5:7),
    winter_months = c(6:8),
    spawning_months = c(10:12),
    spring_months = c(9:11)
  )
}

switch_names <- function(x) {
  
  out <- rep(NA, length(x))
  for (i in seq_along(x)) {
    out[i] <- 
      switch(x[i],
             "MC" = "Maccullochella peelii",
             "Murray Cod" = "Maccullochella peelii",
             "SP" = "Bidyanus bidyanus",
             "Silver perch" = "Bidyanus bidyanus",
             "TC" = "Maccullochella macquariensis",
             "YB" = "Macquaria ambigua",
             "GP" = "Macquaria ambigua",
             "Golden Perch" = "Macquaria ambigua",
             "RF" = "Melanotaenia fluviatilis")
  }
  
  out
  
}

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
