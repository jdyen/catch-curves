# set working directory
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(lme4)

# load some helper functions
source("code/plot-helpers.R")

# load compiled survey data
alldat <- readRDS("data/data-loaded-Jun19.rds")

### MISSING TEMP DATA FOR OVENS@WANG AND KING (COULD SUB OVENS DATA FOR KING)

# filter survey data to MC
to_keep <- alldat$Scientific.Name == "Maccullochella peelii"
alldat <- alldat[to_keep, ]
alldat$Common.Name <- NULL

# need to load flow data
flow_data <- readRDS("data/flow-data-loaded-Jun19.rds")

# filter to MC
flow_data <- flow_data[to_keep, ]

# optional: subset to systems of interest
systems_to_keep <- c("BROKEN", "GOULBURN",
                     "KING", "LOWERMURRAY", "OVENS")
flow_data <- flow_data[alldat$SYSTEM %in% systems_to_keep, ]
alldat <- alldat[alldat$SYSTEM %in% systems_to_keep, ]

# define a function to convert size to age based on model in Todd & Koehn 2008
inverse_growth <- function(x, length_inf, time_zero, k_param, c_param) {
  
  par1 <- 1 / k_param
  ratio_param <- 1 - (time_zero / length_inf)
  par2 <- 1 / (1 - c_param * ratio_param)
  par3 <- (length_inf - time_zero) / (length_inf - x)
  par4 <- c_param * ratio_param
  
  par1 * log(par2 * (par3 - par4))
  
}

# data prep
survey_data <- data.frame(length = alldat$totallength / 10,
                          system = alldat$SYSTEM,
                          site = alldat$SITE_CODE,
                          year = alldat$YEAR,
                          dataset = alldat$dataset)
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

# CT breaks: works out to (0, 200, 300, 400, 500)
# CT v2 (with -0.4:n_age breaks): (0, 150, 275, 380, 470, 545, above)
#  - maybe not if rounding rather than flooring
# size_breaks <- c(-10, 120, 175, 270, 350, 1400)
# size_breaks <- c(-10, 115, 165, 250, 340, 1400)
# age_vec <- cut(10 * survey_data$length, breaks = size_breaks,
#                labels = FALSE)
# age_vec <- age_vec - 1

# pull out indices for random effects
nsystem <- max(survey_data$system)
nsite <- max(survey_data$site)
nyear <- max(survey_data$year)
ndataset <- max(survey_data$dataset)

# settings for linear model
max_age <- ceiling(max(age_vec))

# we need binned data by site and year
age_seq <- seq(-0.4, max_age + 1, by = 1)
hist_fn <- function(x, breaks) {
  hist(x, breaks = breaks, plot = FALSE)$counts
}
age_counts <- tapply(age_vec, list(survey_data$system, survey_data$year), hist_fn, breaks = age_seq)
age_mat <- do.call(rbind, c(age_counts))
age_mat <- do.call(rbind, c(age_counts))
age_mat <- age_mat[, 1:4]

# include flow in year of survey only, assume cohort effects are captured in
#   survival link among years (flow affects YOY, which carries through to later years)
rrang_compiled <- tapply(flow_data$rrang_spwn_mld, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
rrang_ym1_compiled <- tapply(flow_data[, "rrang_spwn_mld_ym1"], list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
psprw_compiled <- tapply(flow_data$prop_spr_lt_win, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
psumw_compiled <- tapply(flow_data$prop_sum_lt_win, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
psprw_ym1_compiled <- tapply(flow_data$prop_spr_lt_win_ym1, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
psumw_ym1_compiled <- tapply(flow_data$prop_sum_lt_win_ym1, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
minwin_compiled <- tapply(flow_data$numlow_days, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
spwn_temp_compiled <- tapply(flow_data$spwntmp_c, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)

# if wanting to work with flows relative to current winter rather than relative to long-term winter flows
# psprw_compiled <- tapply(flow_data$prop_spr_win, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
# psumw_compiled <- tapply(flow_data$prop_sum_win, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
# psprw_ym1_compiled <- tapply(flow_data$prop_spr_win_ym1, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
# psumw_ym1_compiled <- tapply(flow_data$prop_sum_win_ym1, list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
# psprw_compiled[is.infinite(psprw_compiled)] <- max(psprw_compiled[is.finite(psprw_compiled)], na.rm = TRUE)
# psumw_compiled[is.infinite(psumw_compiled)] <- max(psumw_compiled[is.finite(psumw_compiled)], na.rm = TRUE)
# psprw_ym1_compiled[is.infinite(psprw_ym1_compiled)] <- max(psprw_ym1_compiled[is.finite(psprw_ym1_compiled)], na.rm = TRUE)
# psumw_ym1_compiled[is.infinite(psumw_ym1_compiled)] <- max(psumw_ym1_compiled[is.finite(psumw_ym1_compiled)], na.rm = TRUE)

# pull out system and year info
system_info <- rep(rownames(age_counts), times = nyear)
year_info <- rep(colnames(age_counts), each = nsystem)

# subset to observed years and systems
to_keep <- !sapply(c(age_counts), is.null)
system_info <- as.numeric(system_info[to_keep])
year_info <- as.numeric(year_info[to_keep])
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

# round(cor(cbind(rrang_std, rrang_ym1_std,
#                 psprw_std, psprw_ym1_std,
#                 psumw_std, psumw_ym1_std,
#                 minwin_std, spwntmp_std),
#           use = "complete"),
#       2)

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
                       age_factor = factor(rep(seq_len(ncol(age_mat)), each = nrow(age_mat))))
data_set$system_vec <- factor(data_set$system_vec)

# hack to fill temp data for now
data_set$spwntmp_vec[is.na(data_set$spwntmp_vec)] <- mean(data_set$spwntmp_vec, na.rm = TRUE)

# fit a model
## MAKE system fixed? Then mdoel is asking: relative to age class and system, are some years higher or lower than others?
## Option: model responses at level of sites not systems? QUesiton then becomes about sites, but need site-level flow data then.
mod <- glmer(response_vec ~ age_predictor + system_vec +
               rrang_vec +
               rrang_vec : age_factor +
               rrang_ym1_vec +
               rrang_ym1_vec : age_factor +
               psprw_vec + 
               psprw_vec : age_factor + 
               psprw_ym1_vec + 
               psprw_ym1_vec : age_factor +
               psumw_vec + 
               psumw_vec : age_factor + 
               psumw_ym1_vec + 
               psumw_ym1_vec : age_factor +
               minwin_vec + spwntmp_vec +
               minwin_vec : age_factor +
               spwntmp_vec : age_factor +
               (1 | year_vec) +
               (1 | cohort_vec), 
             data = data_set,
             family = stats::poisson,
             na.action = na.fail)

mod_dredge <- dredge(mod, fixed = ~ age_predictor, beta = "sd")

