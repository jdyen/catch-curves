# fit a catch-curve regression model
fit_ccr <- function(
  response,
  length_age_matrix,
  predictors,
  effort,
  system,
  year,
  sys_year_flow,
  include = list(),
  priors = list(),
  mcmc_settings = list(),
  optim_start = FALSE
) {
  
  # unpack priors
  prior_set <- list(sigma = sd(log(response + 1)),
                    sigma_random = 0.5 * sd(log(response + 1)))
  prior_set[names(priors)] <- priors
  
  # unpack mcmc settings
  mcmc_set <- list(n_samples = 2000,
                   warmup = 2000,
                   chains = 4)
  mcmc_set[names(mcmc_settings)] <- mcmc_settings

  # unpack inclusions
  inclusions <- list(survey = TRUE,
                     sys_flow = TRUE,
                     predictors = TRUE)
  inclusions[names(include)] <- include
    
  # ideally want to keep stored samples down
  if (is.null(mcmc_settings$thin))
    mcmc_set$thin <- ifelse(mcmc_set$n_samples > 1000, floor(mcmc_set$n_samples / 1000), 1)
  
  # are effort data provided?
  if (missing(effort))
    effort <- rep(1, nrow(predictors))
  
  # prepare data from inputs
  data <- prepare_data(length_age_matrix, system, predictors, effort)
  
  # unpack data outputs
  n_len <- data$n_len
  n_age <- data$n_age
  n_system <- data$n_system
  n_survey <- data$n_survey
  n_predictors <- data$n_predictors
  n_sys_age <- data$n_sys_age
  sys_vec <- data$sys_vec
  age_vec <- data$age_vec
  survey_vec <- data$survey_vec
  sys_age_vec <- data$sys_age_vec
  effort_vec <- data$effort_vec
  
  # set priors on length-age conversion
  length_age_prior <- zeros(n_len, n_age) + 0.001
  length_age_prior[row(length_age_prior) == col(length_age_prior)] <- 1000
  length_to_age <- dirichlet(alpha = length_age_prior)

  # set likelihood on length-age conversion data
  distribution(length_age_matrix) <-
    multinomial(size = apply(length_age_matrix, 1, sum), prob = length_to_age)
  
  # need to set some variance priors for hierarchical coefficients
  sigma_alpha <- normal(0, prior_set$sigma, truncation = c(0, Inf))
  sigma_beta <- normal(0, prior_set$sigma, truncation = c(0, Inf))
  sigma_survey <- normal(0, prior_set$sigma_random, truncation = c(0, Inf))
  
  # need priors on the regression coefs
  alpha <- normal(0, sigma_alpha, dim = n_system)
  beta <- normal(0, sigma_beta, dim = n_system, truncation = c(-Inf, 0))
  gamma_survey <- normal(0, sigma_survey, dim = n_survey)
  
  # define flow priors (trickier than the others because we have crossed, shared variances)
  if (inclusions$sys_flow) {
    sigma_pred_sub <- normal(0, prior_set$sigma, dim = c(1, n_predictors), truncation = c(0, Inf))
    sigma_pred <- do.call(rbind, lapply(seq_len(n_system), function(i) sigma_pred_sub))
    pred_effects <- normal(0, sigma_pred, dim = c(n_system, n_predictors))
  } else {
    sigma_pred <- normal(0, prior_set$sigma, dim = n_predictors, truncation = c(0, Inf))
    pred_effects <- normal(0, sigma_pred, dim = n_predictors)
  }
  
  # define linear predictor: includes system and age specific predictor effects, with random intercepts for year and cohort
  mu <- alpha[sys_vec] + beta[sys_vec] * age_vec + log(effort_vec)

  # include survey ids as random?    
  if (inclusions$survey)
    mu <- mu + gamma_survey[survey_vec]
  
  # pull out rows based on modified sys_year to account for staggered flows by age
  sys_year_observed <- paste(rep(system, n_age), c(sapply(seq_len(n_age), function(x) year - x + 1)), sep = "_")
  expanded_rows <- match(sys_year_observed, paste(sys_year_flow$system, sys_year_flow$year, sep = "_"))
  predictors_expanded <- predictors[expanded_rows, ]
  
  # include flow effects by system?
  if (inclusions$predictors) {
    if (inclusions$sys_flow) {
      mu <- mu + rowSums(pred_effects[sys_vec, ] * predictors_expanded)
    } else {
      mu <- mu + predictors_expanded %*% pred_effects
    }
  }
  
  # put back on original scale
  modelled_ages <- exp(mu)
  dim(modelled_ages) <- c(n_survey, n_age)
  
  # calculate modelled sizes from ages
  age_to_length <- t(length_to_age)
  modelled_lengths <- modelled_ages %*% age_to_length
  
  # flatten the modelled and observed lengths
  length_vec <- c(modelled_lengths)
  response_vec <- c(response)
  
  # set likelihood for modelled sizes
  subset <- rep(seq_len(n_age), each = n_survey) > 0
  distribution(response_vec[subset]) <- poisson(length_vec[subset])
  
  # compile model
  if (inclusions$survey) {
    mod <- model(alpha,
                 beta,
                 pred_effects, 
                 length_to_age,
                 gamma_survey)
  } else {
    mod <- model(alpha,
                 beta,
                 pred_effects,
                 length_to_age)
  }
  
  # set initial values
  if (optim_start) {
    inits <- list()
    for (i in seq_len(mcmc_set$chains)) {
      opt_start <- opt(mod, max_iterations = 1000, tolerance = 1e-4,
                       initial_values = initials(alpha = rnorm(length(alpha)),
                                                 beta = rep(-1, length(beta))))
      inits[[i]] <- initials(alpha = opt_start$par$alpha,
                             beta = opt_start$par$beta,
                             pred_effects = opt_start$par$pred_effects)
    }
  } else {
    inits <- lapply(seq_len(mcmc_set$chains), function(i) initials(alpha = rep(-5, length(alpha))))
  }
  
  # sample from model
  draws <- mcmc(mod,
                sampler = hmc(),
                n_samples = mcmc_set$n_samples,
                warmup = mcmc_set$warmup,
                thin = mcmc_set$thin,
                chains = mcmc_set$chains,
                initial_values = inits)
  
  # compile data into a clean output obj
  data_list <- list(response = response,
                    length_age_matrix = length_age_matrix,
                    system = system,
                    year = year,
                    predictors = predictors_expanded,
                    effort = effort)
  
  # compile fitted
  out <- list(data = data_list,
              draws = draws,
              priors = priors,
              mcmc_settings = mcmc_settings,
              include = inclusions)
  
  # set class
  out <- as.ccr_model(out)
  
  # return
  out
  
}

# prepare data for a CCR model
prepare_data <- function(length_age_matrix, system, predictors, effort) {
  
  # we need a few indices to keep track of things
  n_len <- nrow(length_age_matrix)
  n_age <- ncol(length_age_matrix)
  
  # how many systems/years are we including?
  n_system <- max(system)
  n_predictors <- ncol(predictors)
  
  # we need a system x age vector to flatten everything out
  sys_vec <- rep(system, times = n_age)
  age_vec <- rep(seq_len(n_age), each = length(system))
  survey_vec <- rep(seq_along(system), times = n_age)
  n_survey <- max(survey_vec)
  sys_age_vec <- n_age * (sys_vec - 1) + age_vec
  n_sys_age <- max(sys_age_vec)
  
  # we also need to flatten out the effort data
  effort_vec <- rep(effort, times = n_age)
  
  # compile into a clean object
  list(n_len = n_len,
       n_age = n_age,
       n_system = n_system,
       n_predictors = n_predictors,
       sys_vec = sys_vec,
       age_vec = age_vec,
       survey_vec = survey_vec,
       n_survey = n_survey,
       sys_age_vec = sys_age_vec,
       n_sys_age = n_sys_age,
       effort_vec = effort_vec)
  
}

# internal function: create ccr_model object
as.ccr_model <- function(object) {
  class(object) <- c("ccr_model", class(object))
  object
}
