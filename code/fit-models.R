# set working directory
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(rstanarm)

# load some helper functions
source("code/helpers.R")

# load compiled survey data
alldat <- readRDS("data/data-loaded-Jul19.rds")

# filter survey data to MC
to_keep <- alldat$Scientific.Name == "Maccullochella peelii"
alldat <- alldat[to_keep, ]
alldat$Common.Name <- NULL

# need to load flow data
flow_data <- readRDS("data/flow-data-loaded-Jul19.rds")

# filter to MC
flow_data <- flow_data[to_keep, ]

# optional: subset to systems of interest
systems_to_keep <- c("BROKEN", "GOULBURN",
                     "KING", "LOWERMURRAY", "OVENS")
flow_data <- flow_data[alldat$SYSTEM %in% systems_to_keep, ]
alldat <- alldat[alldat$SYSTEM %in% systems_to_keep, ]

# data prep
survey_data <- data.frame(length = alldat$totallength / 10,
                          system = alldat$SYSTEM,
                          site = alldat$SITE_CODE,
                          year = alldat$YEAR,
                          dataset = alldat$dataset,
                          effort = alldat$total_no_passes * alldat$seconds)
flow_data <- flow_data[apply(survey_data, 1, function(x) !any(is.na(x))), ]
survey_data <- survey_data[apply(survey_data, 1, function(x) !any(is.na(x))), ]
survey_data$system <- as.integer(as.factor(survey_data$system))
survey_data$site <- as.integer(as.factor(survey_data$site))
survey_data$year <- as.integer(as.factor(survey_data$year))
survey_data$dataset <- as.integer(as.factor(survey_data$dataset))

# convert observed lengths to ages
len_par <- 150
time_par <- 6
k_par <- 0.0011
c_par <- -103
age_vec <- inverse_growth(survey_data$length,
                          len_par, time_par, k_par, c_par)
age_vec[age_vec < 0] <- 0

# pull out indices for random effects
nsystem <- max(survey_data$system)
nyear <- max(survey_data$year)

# settings for linear model
max_age <- ceiling(max(age_vec))

# we need binned data by site and year
age_seq <- seq(-0.4, max_age + 1, by = 1)
age_counts <- tapply(age_vec, list(survey_data$system, survey_data$year), hist_fn, breaks = age_seq)
age_mat <- do.call(rbind, c(age_counts))
age_mat <- do.call(rbind, c(age_counts))
age_mat <- age_mat[, 1:4]

# need to bin survey effort as well
effort <- tapply(survey_data$effort, list(survey_data$system, survey_data$year), mean)

# include flow in year of survey only, assume cohort effects are captured in
#   survival link among years (flow affects YOY, which carries through to later years)
rrang_compiled <- tapply(flow_data$rrang_spwn_mld, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
rrang_ym1_compiled <- tapply(flow_data$rrang_spwn_mld_ym1, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
psprw_compiled <- tapply(flow_data$prop_spr_lt_win, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
psumw_compiled <- tapply(flow_data$prop_sum_lt_win, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
psprw_ym1_compiled <- tapply(flow_data$prop_spr_lt_win_ym1, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
psumw_ym1_compiled <- tapply(flow_data$prop_sum_lt_win_ym1, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
minwin_compiled <- tapply(flow_data$numlow_days, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
spwn_temp_compiled <- tapply(flow_data$spwntmp_c, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)

# replace King temperature data (all missing) with Ovens data
spwn_temp_compiled[3, ]  <- spwn_temp_compiled[5, ]

## MISSING temperature data for early Murray years
# replace with mean from 2003-2010
spwn_temp_compiled[4, 1:4] <- mean(spwn_temp_compiled[4, 5:10])

# pull out system and year info
system_info <- rep(rownames(age_counts), times = nyear)
year_info <- rep(colnames(age_counts), each = nsystem)

# subset to observed years and systems
to_keep <- !sapply(c(age_counts), is.null)
system_info <- as.numeric(system_info[to_keep])
year_info <- as.numeric(year_info[to_keep])
effort <- effort[to_keep]
rrang_compiled <- rrang_compiled[to_keep]
rrang_ym1_compiled <- rrang_ym1_compiled[to_keep]
psprw_compiled <- psprw_compiled[to_keep]
psumw_compiled <- psumw_compiled[to_keep]
psprw_ym1_compiled <- psprw_ym1_compiled[to_keep]
psumw_ym1_compiled <- psumw_ym1_compiled[to_keep]
minwin_compiled <- minwin_compiled[to_keep]
spwn_temp_compiled <- spwn_temp_compiled[to_keep]

# standardise flow values
rrang_std <- scale(rrang_compiled)
rrang_ym1_std <- scale(rrang_ym1_compiled)
psprw_std <- scale(psprw_compiled)
psumw_std <- scale(psumw_compiled)
psprw_ym1_std <- scale(psprw_ym1_compiled)
psumw_ym1_std <- scale(psumw_ym1_compiled)
minwin_std <- scale(minwin_compiled)
spwntmp_std <- scale(spwn_temp_compiled)

# pull out means and SDs of unscaled flow variables
flow_scales <- list()
flow_scales$rrang_vec$mean <- attributes(rrang_std)$`scaled:center`
flow_scales$rrang_vec$sd <- attributes(rrang_std)$`scaled:scale`
flow_scales$rrang_ym1_vec$mean <- attributes(rrang_ym1_std)$`scaled:center`
flow_scales$rrang_ym1_vec$sd <- attributes(rrang_ym1_std)$`scaled:scale`
flow_scales$psprw_vec$mean <- attributes(psprw_std)$`scaled:center`
flow_scales$psprw_vec$sd <- attributes(psprw_std)$`scaled:scale`
flow_scales$psprw_ym1_vec$mean <- attributes(psprw_ym1_std)$`scaled:center`
flow_scales$psprw_ym1_vec$sd <- attributes(psprw_ym1_std)$`scaled:scale`
flow_scales$psumw_vec$mean <- attributes(psumw_std)$`scaled:center`
flow_scales$psumw_vec$sd <- attributes(psumw_std)$`scaled:scale`
flow_scales$psumw_ym1_vec$mean <- attributes(psumw_ym1_std)$`scaled:center`
flow_scales$psumw_ym1_vec$sd <- attributes(psumw_ym1_std)$`scaled:scale`
flow_scales$minwin_vec$mean <- attributes(minwin_std)$`scaled:center`
flow_scales$minwin_vec$sd <- attributes(minwin_std)$`scaled:scale`
flow_scales$spwntmp_vec$mean <- attributes(spwntmp_std)$`scaled:center`
flow_scales$spwntmp_vec$sd <- attributes(spwntmp_std)$`scaled:scale`

# create a matrix of indices identifying cohorts
cohort_mat <- matrix(NA, nrow = length(system_info), ncol = ncol(age_mat))
current_max <- 0
for (i in seq_len(nsystem)) {
  sys_sub <- system_info == i
  year_sort <- year_info[sys_sub]
  cohort_tmp <- matrix(NA, nrow = sum(sys_sub), ncol = ncol(age_mat))
  cohort_tmp[1, ] <- rev(seq_len(ncol(age_mat)))
  for (j in seq_len(sum(sys_sub))[-1])
    cohort_tmp[j, ] <- cohort_tmp[j - 1, ] + 1
  cohort_tmp <- cohort_tmp + current_max
  current_max <- max(cohort_tmp)
  cohort_mat[which(sys_sub)[order(year_sort)], ] <- cohort_tmp
}

# now we need to create response and predictor variables
data_set <- data.frame(age_predictor = rep(seq_len(ncol(age_mat)), each = nrow(age_mat)),
                       system_vec = rep(system_info, times = ncol(age_mat)),
                       year_vec = rep(year_info, times = ncol(age_mat)),
                       rrang_vec = rep(rrang_std, times = ncol(age_mat)),
                       rrang_ym1_vec = rep(rrang_ym1_std, times = ncol(age_mat)),
                       psprw_vec = rep(psprw_std, times = ncol(age_mat)),
                       psprw_ym1_vec = rep(psprw_ym1_std, times = ncol(age_mat)),
                       psumw_vec = rep(psumw_std, times = ncol(age_mat)),
                       psumw_ym1_vec = rep(psumw_ym1_std, times = ncol(age_mat)),
                       minwin_vec = rep(minwin_std, times = ncol(age_mat)),
                       spwntmp_vec = rep(spwntmp_std, times = ncol(age_mat)),
                       cohort_vec = c(cohort_mat),
                       response_vec = c(age_mat),
                       ncohort = length(unique(c(cohort_mat))),
                       age_factor = factor(rep(seq_len(ncol(age_mat)), each = nrow(age_mat))),
                       sampling_effort = effort)
data_set$system_vec <- factor(data_set$system_vec)

# mcmc settings
n_iter <- 10000
n_chains <- 4
n_cores <- n_chains

# create freq-hist plots

## SWITCH FROM LT WINTER TO LT ANNUAL MEDIAN
# TEMP ONLY FOR SPAWNING -- CAN WE FOCUS ON YOY ONLY?
# CHANGE IN FLOW SPR/SUM IS ABOUT SPAWNING
## WINTER IS AFFECTING 1YO-up because flow is pre-spawning.

## COUDL FIT MULTIPLE MODELS:
## YOY - SPRING FLOWS ,SUMMER FLOWS, SPWN_TMP, YM1_MAX flows.
## SURVIVAL MODEL

## RRANG IS NOW PROPORTIONAL CHANGE (max / min)
## LT_WIN is now LT_MEDIAN over all months (for psprw and psumw)

# COUDL ADD YM2 if fits? (probably not; could use some average of two years??)

# fit the full model with all predictors
mod_full <- stan_glmer(response_vec ~ age_predictor + 
                         (rrang_vec + rrang_ym1_vec +
                            psprw_vec + psprw_ym1_vec +
                            psumw_vec + psumw_ym1_vec + 
                            (psumw_vec ^ 2) + 
                            minwin_vec + spwntmp_vec | system_vec) +
                         (-1 +
                            rrang_vec + rrang_ym1_vec +
                            psprw_vec + psprw_ym1_vec +
                            psumw_vec + psumw_ym1_vec +
                            (psumw_vec ^ 2) +
                            minwin_vec + spwntmp_vec | age_factor) +
                         (1 | year_vec) +
                         (1 | cohort_vec) +
                         offset(sampling_effort),
                       iter = n_iter, chains = n_chains,
                       data = data_set,
                       family = stats::poisson, cores = n_cores)

# fit a reduced model without age-specific flow effects
mod_noage <- stan_glmer(response_vec ~ age_predictor + 
                          (rrang_vec + rrang_ym1_vec +
                             psprw_vec + psprw_ym1_vec +
                             psumw_vec + psumw_ym1_vec + 
                             minwin_vec + spwntmp_vec | system_vec) +
                          (1 | year_vec) +
                          (1 | cohort_vec),
                        iter = n_iter, chains = n_chains,
                        data = data_set,
                        family = stats::poisson, cores = n_cores)

# fit a reduced model without system-specific flow effects
mod_nosys <- stan_glmer(response_vec ~ age_predictor + 
                          (-1 + rrang_vec + rrang_ym1_vec +
                             psprw_vec + psprw_ym1_vec +
                             psumw_vec + psumw_ym1_vec + 
                             minwin_vec + spwntmp_vec | age_factor) +
                          (1 | year_vec) +
                          (1 | cohort_vec),
                        iter = n_iter, chains = n_chains,
                        data = data_set,
                        family = stats::poisson, cores = n_cores)

# fit a reduced model without system- or age-specific flow effects
mod_nosys_noage <- stan_glmer(response_vec ~ age_predictor + 
                                rrang_vec + rrang_ym1_vec +
                                psprw_vec + psprw_ym1_vec +
                                psumw_vec + psumw_ym1_vec + 
                                minwin_vec + spwntmp_vec +
                                (1 | year_vec) +
                                (1 | cohort_vec),
                              iter = n_iter, chains = n_chains,
                              data = data_set,
                              family = stats::poisson, cores = n_cores)

# fit a reduced model without flow predictors
mod_noflow <- stan_glmer(response_vec ~ age_predictor + 
                           (1 | system_vec) +
                           (1 | year_vec) +
                           (1 | cohort_vec),
                         iter = n_iter, chains = n_chains,
                         data = data_set,
                         family = stats::poisson, cores = n_cores)

# save fitted model/s
saveRDS(mod_full, file = "outputs/fitted/full-model.rds")
saveRDS(mod_noflow, file = "outputs/fitted/noflow-model.rds")
saveRDS(mod_noage, file = "outputs/fitted/noage-model.rds")
saveRDS(mod_nosys, file = "outputs/fitted/nosys-model.rds")
saveRDS(mod_nosys_noage, file = "outputs/fitted/nosysnoage-model.rds")

# save data and related
saveRDS(data_set, file = "outputs/fitted/data-set.rds")
additional_data <- list(flow_scales = flow_scales,
                        age_mat = age_mat,
                        nyear = nyear,
                        nsystem = nsystem,
                        system_info = system_info,
                        year_info = year_info,
                        cohort_mat = cohort_mat)
saveRDS(additional_data, file = "outputs/fitted/additional-data.rds")
