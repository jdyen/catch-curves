# fit a catch-curve regression model
fit_ccr <- function(response, length_age_matrix,
                    predictors, effort,
                    system, year,
                    priors = list(),
                    mcmc_settings = list()) {
  
  # unpack priors
  prior_set <- list(sigma = 2,
                    sigma_random = 1)
  prior_set[names(priors)] <- priors
  
  # unpack mcmc settings
  mcmc_set <- list(n_samples = 2000,
                   warmup = 2000,
                   chains = 4)
  mcmc_set[names(mcmc_settings)] <- mcmc_settings
  
  # ideally want to keep stored samples down
  if (is.null(mcmc_settings$thin))
    mcmc_set$thin <- ifelse(mcmc_set$n_samples > 1000, floor(mcmc_set$n_samples / 1000), 1)
  
  # are effort data provided?
  if (missing(effort))
    effort <- rep(1, nrow(predictors))
  
  # prepare data from inputs
  data <- prepare_data(length_age_matrix, system, year, predictors, effort)
  
  # unpack data outputs
  n_len <- data$n_len
  n_age <- data$n_age
  n_system <- data$n_system
  n_year <- data$n_year
  n_predictors <- data$n_predictors
  sys_vec <- data$sys_vec
  age_vec <- data$age_vec
  year_vec <- data$year_vec
  survey_vec <- data$survey_vec
  n_survey <- data$n_survey
  sys_age_vec <- data$sys_age_vec
  n_sys_age <- data$n_sys_age
  cohort_vec <- data$cohort_vec
  n_cohort <- data$n_cohort
  effort_vec <- data$effort_vec
  
  # set priors on length-age conversion
  length_age_prior <- zeros(n_len, n_age) + 0.001
  length_age_prior[row(length_age_prior) == col(length_age_prior)] <- 100
  length_to_age <- dirichlet(alpha = length_age_prior)
  
  # set likelihood on length-age conversion data
  distribution(length_age_matrix) <-
    multinomial(size = apply(length_age_matrix, 1, sum), prob = length_to_age)
  
  # need to set some variance priors for hierarchical coefficients
  sigma_alpha <- normal(0, prior_set$sigma, truncation = c(0, Inf))
  sigma_beta <- normal(0, prior_set$sigma, truncation = c(0, Inf))
  sigma_year <- normal(0, prior_set$sigma_random, truncation = c(0, Inf))
  sigma_cohort <- normal(0, prior_set$sigma_random, truncation = c(0, Inf))
  
  # need priors on the regression coefs
  alpha <- normal(0, sigma_alpha, dim = n_system)
  beta <- normal(0, sigma_beta, dim = n_system, truncation = c(-Inf, 0))
  gamma_year <- normal(0, sigma_year, dim = n_year)
  gamma_cohort <- normal(0, sigma_cohort, dim = n_cohort)
  
  # define flow priors (trickier than the others because we have crossed, shared variances)
  sigma_pred_sub <- normal(0, prior_set$sigma, dim = c(n_age, n_predictors), truncation = c(0, Inf))
  sigma_pred <- do.call(rbind, lapply(seq_len(n_system), function(i) sigma_pred_sub))
  pred_effects <- normal(0, sigma_pred, dim = c(n_sys_age, n_predictors))
  
  # define linear predictor: includes system and age specific predictor effects, with random intercepts for year and cohort
  mu <- alpha[sys_vec] + beta[sys_vec] * age_vec +
    gamma_cohort[cohort_vec] + gamma_year[year_vec] +
    rowSums(pred_effects[sys_age_vec, ] * predictors[survey_vec, ]) +
    log(effort_vec)
  
  # put back on original scale
  modelled_ages <- exp(mu)
  dim(modelled_ages) <- c(n_survey, n_age)
  
  # calculate modelled sizes from ages
  # age_to_length <- t(length_to_age)
  modelled_lengths <- t(length_to_age %*% t(modelled_ages))
  # modelled_lengths <- modelled_ages %*% age_to_length
  
  # flatten the modelled and observed lengths
  length_vec <- c(modelled_lengths)
  response_vec <- c(response)
  
  # set likelihood for modelled sizes
  distribution(response_vec) <- poisson(length_vec)
  
  # compile and sample from model
  mod <- model(alpha, beta, pred_effects, 
               length_to_age,
               gamma_cohort, gamma_year)
  draws <- mcmc(mod, n_samples = mcmc_set$n_samples,
                warmup = mcmc_set$warmup,
                thin = mcmc_set$thin,
                chains = mcmc_set$chains)
  
  # compile data into a clean output obj
  data_list <- list(response = response,
                    length_age_matrix = length_age_matrix,
                    system = system,
                    year = year,
                    predictors = predictors,
                    effort = effort,
                    cohort_mat = data$cohort_mat)
  
  # compile fitted
  out <- list(data = data_list,
              draws = do.call(rbind, draws),
              priors = priors,
              mcmc_settings = mcmc_settings)
  
  # set class
  out <- as.ccr_model(out)
  
  # return
  out
  
}

# prepare data for a CCR model
prepare_data <- function(length_age_matrix, system, year, predictors, effort) {
  
  # we need a few indices to keep track of things
  n_len <- nrow(length_age_matrix)
  n_age <- ncol(length_age_matrix)
  
  # how many systems/years are we including?
  n_system <- max(system)
  n_year <- max(year)
  n_predictors <- ncol(predictors)
  
  # we need a system x age vector to flatten everything out
  sys_vec <- rep(system, times = n_age)
  age_vec <- rep(seq_len(n_age), each = length(system))
  year_vec <- rep(year, times = n_age)
  survey_vec <- rep(seq_along(system), times = n_age)
  n_survey <- max(survey_vec)
  sys_age_vec <- n_age * (sys_vec - 1) + age_vec
  n_sys_age <- max(sys_age_vec)
  
  # we also need to flatten out the effort data
  effort_vec <- rep(effort, times = n_age)
  
  # tricky bit: create a matrix of indices identifying cohorts
  cohort_mat <- matrix(NA, nrow = length(system), ncol = n_age)
  current_max <- 0
  for (i in seq_len(n_system)) {
    if (any(system == i)) {
      sys_sub <- system == i
      year_sub <- year[sys_sub]
      idx <- order(year_sub)
      year_sort <- year_sub[idx]
      cohort_tmp <- matrix(NA, nrow = sum(sys_sub), ncol = n_age)
      cohort_tmp[idx[1], ] <- rev(seq_len(n_age))
      year_diffs <- diff(year_sort)
      for (j in seq_len(sum(sys_sub))[-1])
        cohort_tmp[idx[j], ] <- cohort_tmp[idx[j - 1], ] + year_diffs[j - 1]
      cohort_tmp <- cohort_tmp + current_max
      current_max <- max(cohort_tmp)
      cohort_mat[which(sys_sub)[order(year_sort)], ] <- cohort_tmp
    }
  }
  n_cohort <- max(cohort_mat, na.rm = TRUE)
  cohort_vec <- c(cohort_mat)
  
  # compile into a clean object
  list(n_len = n_len,
       n_age = n_age,
       n_system = n_system,
       n_year = n_year,
       n_predictors = n_predictors,
       sys_vec = sys_vec,
       age_vec = age_vec,
       year_vec = year_vec,
       survey_vec = survey_vec,
       n_survey = n_survey,
       sys_age_vec = sys_age_vec,
       n_sys_age = n_sys_age,
       cohort_vec = cohort_vec,
       n_cohort = n_cohort,
       cohort_mat = cohort_mat,
       effort_vec = effort_vec)
  
}

# internal function: create ccr_model object
as.ccr_model <- function(object) {
  class(object) <- c("ccr_model", class(object))
  object
}
