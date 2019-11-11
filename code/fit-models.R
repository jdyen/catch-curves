# set working directory
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(rstanarm)
library(lubridate)

# load some helper functions
source("code/helpers.R")

# load compiled survey data
alldat <- readRDS("data/data-loaded-Oct19.rds")

# filter survey data to MC
to_keep <- alldat$scientific_name == "Maccullochella peelii"
alldat <- alldat[to_keep, ]

# load otolith data
oti_data <- readRDS("data/otolith-data-loaded.rds")

# filter to MC
oti_data <- oti_data[oti_data$SPECIES == "Maccullochella peelii", ]

# need to load flow data
flow_data <- readRDS("data/flow-data-loaded-Oct19.rds")

# filter to MC
flow_data <- flow_data[to_keep, ]

# optional: subset to systems of interest
systems_to_keep <- c("broken", "goulburn", "king", "murray", "ovens")
flow_data <- flow_data[alldat$system %in% systems_to_keep, ]
alldat <- alldat[alldat$system %in% systems_to_keep, ]

# hack for now because ovens site data incomplete
alldat$site[is.na(alldat$site)] <- 1

# prepare survey data
survey_data <- data.frame(length_mm = alldat$total_length_mm,
                          system = alldat$system,
                          system_id = rebase_index(alldat$system),
                          site = alldat$site,
                          site_id = rebase_index(alldat$site),
                          date = alldat$date_formatted,
                          dataset = alldat$dataset,
                          dataset_id = rebase_index(alldat$dataset),
                          effort = alldat$ef_seconds_total,
                          stringsAsFactors = FALSE)

# add years
survey_data$year <- year(survey_data$date)
survey_data$year_id <- rebase_index(survey_data$year)

# fill NAs in length data for padded rows
survey_data$length_mm[survey_data$dataset == "PADDED_TO_GET_FLOW_YEARS"] <- 0

# remove any rows with missing data
flow_data <- flow_data[apply(survey_data, 1, function(x) !any(is.na(x))), ]
survey_data <- survey_data[apply(survey_data, 1, function(x) !any(is.na(x))), ]

# remove padded rows from survey data to calculate age matrix
survey_data_filtered <- survey_data[survey_data$dataset != "PADDED_TO_GET_FLOW_YEARS", ]

# convert observed lengths to ages
len_par <- 150
time_par <- 6
k_par <- 0.0011
c_par <- -103
age_vec <- inverse_growth(survey_data_filtered$length_mm / 10,
                          len_par, time_par, k_par, c_par)
age_vec[age_vec < 0] <- 0

# pull out indices for random effects
nsystem <- max(survey_data_filtered$system_id)
nyear <- max(survey_data_filtered$year_id)

# settings for linear model
max_age <- ceiling(max(age_vec))

# we need binned data by site and year
age_seq <- seq(-0.4, max_age + 1, by = 1)
age_counts <- tapply(age_vec, list(survey_data_filtered$system_id, survey_data_filtered$year_id), hist_fn, breaks = age_seq)
age_mat <- do.call(rbind, c(age_counts))

# how many age classes to keep (default = 6)?
n_age <- 6

# pull out adult abundances
age_vec_full <- inverse_growth(survey_data$length_mm / 10,
                               len_par, time_par, k_par, c_par)
age_vec_full[age_vec_full < 0] <- 0
age_counts_full <- tapply(age_vec_full, list(survey_data$system_id, survey_data$year_id), hist_fn, breaks = age_seq)
age_mat_full <- do.call(rbind, c(age_counts_full))
adult_catch <- apply(age_mat_full[, 5:ncol(age_mat_full)], 1, sum)

# cut off at 5 year olds based on linearity of count ~ age
age_mat <- age_mat[, seq_len(n_age)]

# how much effort went into each survey? This should be summed over all sits in a given system
#   (but not double-count effort for each individual fish)
sys_site <- paste(survey_data_filtered$system_id, survey_data_filtered$site, sep = "_")
effort_site <- tapply(survey_data_filtered$effort,
                      list(sys_site, survey_data_filtered$year_id),
                      unique)
year_ids <- rep(colnames(effort_site), each = nrow(effort_site))
system_ids <- sapply(strsplit(rownames(effort_site), "_"), function(x) x[1])
system_ids <- rep(system_ids, times = ncol(effort_site))
effort_site <- sapply(c(effort_site), sum)
effort <- tapply(effort_site, list(system_ids, year_ids), sum)
effort <- effort[, order(as.numeric(colnames(effort)))]
sys_site_full <- paste(survey_data$system, survey_data$site, sep = "_")
effort_full <- tapply(survey_data$effort,
                      list(sys_site_full, survey_data$year_id),
                      unique)
year_ids_full <- rep(colnames(effort_full), each = nrow(effort_full))
system_ids_full <- sapply(strsplit(rownames(effort_full), "_"), function(x) x[1])
system_ids_full <- rep(system_ids_full, times = ncol(effort_full))
effort_full <- sapply(c(effort_full), sum)
effort_full <- tapply(effort_full, list(system_ids_full, as.integer(year_ids_full)), sum)

# want to identify years with synthetic data because effort and catch is incorrect
synth_data <- survey_data[survey_data$dataset == "PADDED_TO_GET_FLOW_YEARS", ]
all_years <- tapply(synth_data$year, synth_data$system, function(x) sort(unique(x)))
all_years$murray <- sort(c(all_years$murray, 2011:2018))
observed_years <- tapply(survey_data_filtered$year, survey_data_filtered$system, function(x) sort(unique(x)))
years_of_surveys <- mapply(`%in%`, all_years, observed_years)
years_of_surveys_matrix <- t(effort_full)
years_of_surveys_matrix[years_of_surveys_matrix > 0] <- do.call(c, years_of_surveys)
years_of_surveys_matrix <- t(years_of_surveys_matrix)

# compile flow predictors
var_sets <- list(c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "prop_sum_lt_win_ym1", "prop_sum_win_sq", "maxan_mld", "maxan_mld_ym1", "spwntmp_c"),
                 c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "prop_sum_win_sq", "maxan_mld", "maxan_mld_ym1", "spwntmp_c"),
                 c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "maxan_mld", "maxan_mld_ym1", "spwntmp_c"),
                 c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "maxan_mld", "spwntmp_c"),
                 c("rrang_spwn_mld", "spwntmp_c"),
                 c("prop_spr_lt_win", "spwntmp_c"),
                 c("prop_sum_lt_win", "spwntmp_c"),
                 c("maxan_mld", "spwntmp_c"))

vars_to_include <- var_sets[[3]]
# vars_to_include <- c("prop_spr_lt_win", "prop_sum_lt_win", "maxan_mld", "spwntmp_c")

flow_compiled <- sapply(vars_to_include,
                        function(x) tapply(get(x, flow_data),
                                           list(survey_data$system, survey_data$year),
                                           mean, na.rm = TRUE))
sys_year <- data.frame(system = rep(sort(unique(survey_data$system_id)), times = length(unique(survey_data$year_id))),
                       year = rep(sort(unique(survey_data$year_id)), each = length(unique(survey_data$system_id))))
sys_year <- sys_year[!is.na(flow_compiled[, 1]), ]
flow_compiled <- flow_compiled[!is.na(flow_compiled[, 1]), ]

# expand system and year vectors
system <- rep(as.integer(as.factor(rownames(effort))), times = ncol(effort))
year <- rep(as.numeric(colnames(effort)), each = nrow(effort))
is_missing <- effort == 0
system <- system[!is_missing]
year <- year[!is_missing]
effort <- effort[!is_missing]

#   - King = Ovens (2004-200? filled with average of following 4 years)
#   - Murray 1993-2002 filled with average of Murray 2003-2008
if ("spwntmp_c" %in% vars_to_include) {
  flow_compiled[sys_year$system == 4 & sys_year$year <= 10, "spwntmp_c"] <-
    mean(flow_compiled[sys_year$system == 4 & sys_year$year %in% c(11:16), "spwntmp_c"])
  flow_compiled[sys_year$system == 3, "spwntmp_c"] <- 
    flow_compiled[match(paste0("5", sys_year$year[sys_year$system == 3]), paste0(sys_year$system, sys_year$year)), "spwntmp_c"]
  flow_compiled[sys_year$system == 3 & sys_year$year <= 15, "spwntmp_c"] <-
    mean(flow_compiled[sys_year$system == 3 & sys_year$year %in% c(16:19), "spwntmp_c"])
}

# standardised flow data
flow_std <- apply(flow_compiled, 2, scale)
flow_scales <- apply(flow_compiled, 2, extract_standards)

# add total abundance as a predictor
years_of_surveys <- c(years_of_surveys_matrix)
years_of_surveys <- years_of_surveys[effort_full > 0]
effort_full <- effort_full[effort_full > 0]
for (i in seq_len(max(sys_year$system))) {
  idx <- sys_year$system == i & years_of_surveys == 0
  idy <- sys_year$system == i & years_of_surveys == 1
  min_year <- min(sys_year$year[idy])
  idy <- sys_year$system == i & sys_year$year %in% c(min_year:(min_year + 3))
  effort_full[idx] <- mean(effort_full[idy])
  adult_catch[idx] <- median(adult_catch[idy])
}
flow_std <- cbind(flow_std, "adult_cpue" = c(scale(adult_catch / effort_full)))
flow_scales <- cbind(flow_scales, "adult_cpue" = c(mean(adult_catch / effort_full), sd(adult_catch / effort_full)))

# indices to use later
n_age <- ncol(age_mat)
n_obs <- nrow(age_mat)

# pull out rows based on modified sys_year to account for staggered flows by age
sys_year_observed <- paste(rep(system, n_age), c(sapply(seq_len(n_age), function(x) year - x + 1)), sep = "_")
expanded_rows <- match(sys_year_observed, paste(sys_year$system, sys_year$year, sep = "_"))
flow_expanded <- flow_std[expanded_rows, ]

# now we need to create response and predictor variables
data_set <- data.frame(age_predictor = rep(seq_len(n_age), each = n_obs),
                       system_vec = factor(rep(system, times = n_age)),
                       year_vec = factor(rep(year - min(year) + 1, times = n_age)),
                       flow_expanded,
                       response_vec = c(age_mat),
                       age_factor = factor(rep(seq_len(n_age), each = n_obs)),
                       sampling_effort = effort / 60)
data_set$survey_vec <- rebase_index(paste(data_set$system_vec, data_set$year_vec, sep = "_"))

# remove early years for all systems except Murray
to_keep <- rep(seq_len(n_age), each = n_obs) > 1
# to_keep[data_set$system_vec == 4] <- TRUE
data_set <- data_set[to_keep, ]

# mcmc settings
n_iter <- 10000
n_chains <- 3
n_cores <- n_chains

# mod_tmp <- list()
# for (i in seq_along(vars_to_include)) {
#   
#   formula_tmp <- paste0("response_vec ~ (age_predictor | system_vec) + ",
#                         paste(paste0("(-1 + ", vars_to_include[i], " | system_vec)"), collapse = " + "),
#                         " + (1 | year_vec)")
#   
#   mod_tmp[[i]] <- lme4::glmer(formula_tmp,
#                               offset = log(sampling_effort),
#                               data = data_set,
#                               family = stats::poisson(),
#                               na.action = "na.fail")
# 
# }
# 
# mod_fixed <- paste("response_vec ~ ", 
#                    paste0("(age_predictor | system_vec) + ",
#                           paste(vars_to_include, collapse = " + ")),
#                      " + (1 | year_vec)",
#                      sep = "")
# 
# mod_test <- lme4::glmer(mod_fixed,
#                         offset = log(sampling_effort),
#                         data = data_set,
#                         family = stats::poisson(),
#                         na.action = "na.fail")
# many_mods <- dredge(mod_test)

vars_to_include <- c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "maxan_mld", "spwntmp_c", "adult_cpue")

sys_term <- paste0("(age_predictor | system_vec) + ",
                   paste(paste0("(-1 + ", vars_to_include, " | system_vec)"), collapse = " + "))

mod_formula <- paste("response_vec ~ ", 
                     sys_term,
                     " + (1 | survey_vec)",
                     sep = "")
# 
# mod_test2 <- lme4::glmer(mod_formula,
#                          offset = log(sampling_effort),
#                          data = data_set,
#                          family = stats::poisson(),
#                          na.action = "na.fail")


# fit the full model with all predictors
mod_full <- stan_glmer(mod_formula,
                       iter = n_iter, chains = n_chains,
                       data = data_set, offset = log(sampling_effort),
                       family = neg_binomial_2(), cores = n_cores)

# fit a reduced model without flow predictors
mod_noflow <- stan_glmer(response_vec ~ (age_predictor | system_vec) +
                           (1 | survey_vec),
                         iter = n_iter, chains = n_chains,
                         data = data_set, offset = log(sampling_effort),
                         family = neg_binomial_2(), cores = n_cores)

# save fitted model/s
saveRDS(mod_full, file = "outputs/fitted/full-model.rds")
saveRDS(mod_noflow, file = "outputs/fitted/noflow-model.rds")

# save data and related
saveRDS(data_set, file = "outputs/fitted/data-set.rds")
additional_data <- list(flow_scales = flow_scales,
                        flow_expanded = flow_expanded,
                        age_mat = age_mat,
                        system = system,
                        year = year)
saveRDS(additional_data, file = "outputs/fitted/additional-data.rds")


## check residuals
## NB rather than Poisson model based on mean_PPD and LOO_IC
# resid(mod_full) ~ fitted(mod_full)
## plot residual image plots without fitted model values
## but based on mod_noflow (i.e., what is the variation we are trying
##  to explain with flow?)
## plot flow associations
## run multiple models to check which is "best"
## calc bayes_r2 values for "best" model.
## Plot variance components to highlight variable flow responses
var_sets <- list(
  c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "maxan_mld", "maxan_mld_ym1", "spwntmp_c"),
  c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "maxan_mld", "spwntmp_c"),
  c("prop_spr_lt_win", "prop_sum_lt_win", "maxan_mld", "spwntmp_c"),
  c("rrang_spwn_mld", "prop_sum_lt_win", "maxan_mld", "spwntmp_c"),
  c("rrang_spwn_mld", "prop_spr_lt_win", "maxan_mld", "spwntmp_c"),
  c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "spwntmp_c"),
  c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "maxan_mld"),
  c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "maxan_mld", "spwntmp_c", "adult_cpue")
)

### ADD IN THREE_VAR COMBOS??

loo_out <- list()
for (i in seq_along(var_sets)) {

  vars_to_include <- var_sets[[i]]
  sys_term <- paste0("(age_predictor | system_vec) + ",
                     paste(paste0("(-1 + ", vars_to_include, " | system_vec)"), collapse = " + "))
  mod_formula <- paste("response_vec ~ ", sys_term, sep = "")
  
  # fit the full model with all predictors
  mod_compare <- stan_glmer(mod_formula,
                            iter = n_iter, chains = n_chains,
                            data = data_set, offset = log(sampling_effort),
                            family = neg_binomial_2(), cores = n_cores)
  
  # how does it stack up?
  loo_out[[i]] <- loo(mod_compare)
  
}

