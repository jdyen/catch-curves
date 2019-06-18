# set working directory
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(greta)

# source helper functions
source("code/greta-helpers.R")
source("code/tf_functions.R")

# load compiled survey data
alldat <- readRDS("data/data-loaded-Feb19.rds")

# filter survey data to MC
to_keep <- alldat$Scientific.Name == "Maccullochella peelii"
alldat <- alldat[to_keep, ]
alldat$Common.Name <- NULL

# optional: subset to systems of interest
# systems_to_keep <- c("BROKEN", "CAMPASPE", "GOULBURN",
#                      "KING", "LODDON", "LOWERMURRAY",
#                      "OVENS")
# alldat <- alldat[alldat$SYSTEM %in% systems_to_keep, ]

# need to load flow data
flow_data <- readRDS("data/flow-data-loaded-Feb19.rds")

# filter to MC
flow_data <- flow_data[to_keep, ]

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
len_par <- normal(150, 10)
time_par <- lognormal(6, 1)
k_par <- lognormal(0.0011, 0.0001)
c_par <- -103
# len_par <- 150
# time_par <- 6
# k_par <- 0.0011
# c_par <- -103
age_vec <- inverse_growth(survey_data$length,
                          len_par, time_par, k_par, c_par)
age_vec[age_vec < 0] <- 0

# pull out indices for random effects
nsystem <- max(survey_data$system)
nsite <- max(survey_data$site)
nyear <- max(survey_data$year)
ndataset <- max(survey_data$dataset)

# settings for linear model
max_age <- ceiling(max(age_vec))

# we need binned data by site and year
age_seq <- seq(0, max_age, by = 1)
hist_fn <- function(x, breaks) {
  hist(x, breaks = breaks, plot = FALSE)$counts
}
age_counts <- tapply(age_vec, list(survey_data$system, survey_data$year), hist_fn, breaks = age_seq)
age_mat <- do.call(rbind, c(age_counts))
age_mat <- do.call(rbind, c(age_counts))
age_mat <- age_mat[, 1:8] ## truncate to four-year olds.

# include flow in year of survey only, assume cohort effects are captured in
#   survival link among years (flow affects YOY, which carries through to later years)
mannf_compiled <- tapply(flow_data[, 1], list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
cvflow_compiled <- tapply(flow_data[, 4], list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
msprf_compiled <- tapply(flow_data[, 2], list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
msumf_compiled <- tapply(flow_data[, 3], list(survey_data$system, survey_data$year), mean, na.rm = TRUE)
covsp_compiled <- tapply(flow_data[, 6], list(survey_data$system, survey_data$year), mean, na.rm = TRUE)

# pull out system and year info
system_info <- rep(rownames(age_counts), times = nyear)
year_info <- rep(colnames(age_counts), each = nsystem)

# subset to observed years and systems
to_keep <- !sapply(c(age_counts), is.null)
system_info <- as.numeric(system_info[to_keep])
year_info <- as.numeric(year_info[to_keep])
mannf_compiled <- mannf_compiled[to_keep]
cvflow_compiled <- cvflow_compiled[to_keep]
msprf_compiled <- msprf_compiled[to_keep]
msumf_compiled <- msumf_compiled[to_keep]
covsp_compiled <- covsp_compiled[to_keep]

# fill King river flow values (need to get data on these)
mannf_compiled[is.nan(mannf_compiled)] <- mean(mannf_compiled, na.rm = TRUE)
cvflow_compiled[is.nan(cvflow_compiled)] <- mean(cvflow_compiled, na.rm = TRUE)
msprf_compiled[is.nan(msprf_compiled)] <- mean(msprf_compiled, na.rm = TRUE)
msumf_compiled[is.nan(msumf_compiled)] <- mean(msumf_compiled, na.rm = TRUE)
covsp_compiled[is.nan(covsp_compiled)] <- mean(covsp_compiled, na.rm = TRUE)

# winter/spring mean or max positively associated with adult migration and condition (pre-spawning)
# previous year's flows (mean daily discharge) positively associated with adult condition
# winter low flows (median or min daily discharge) positively associated with juvenile survival and adult condition
# variability in velocity/height/discharge negatively associated with 
#    nest, egg, larval, and juvenile survival. (Oct-Nov for nests and eggs, Nov-Dec for larvae).
# summer mean/max flows associated with larval and juvenile survival (quadratic response; too much is bad, so is too little)
# years since flood/inundation matters (longer = worse for adults?); can we get bankfull thresholds for all systems?

# standardise flow values: use summer flows and CV of spawning flows (Oct-Dec)
mannf_std <- scale(msumf_compiled)
cvflow_std <- scale(covsp_compiled)

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
age_predictor <- rep(seq_len(ncol(age_mat)), each = nrow(age_mat))
system_vec <- rep(system_info, times = ncol(age_mat))
year_vec <- rep(year_info, times = ncol(age_mat))
mannf_vec <- rep(mannf_std, times = ncol(age_mat))
cvflow_vec <- rep(cvflow_std, times = ncol(age_mat))
cohort_vec <- c(cohort_mat)
response_vec <- c(age_mat)
ncohort <- length(unique(cohort_vec))

# we need variance priors before we can set random effects
sigma_mannf <- normal(0, 10, truncation = c(0, Inf))
sigma_mannf2 <- normal(0, 10, truncation = c(0, Inf))
sigma_cvflow <- normal(0, 10, truncation = c(0, Inf))
sigma_system <- normal(0, 10, truncation = c(0, Inf))
sigma_cohort <- normal(0, 10, truncation = c(0, Inf))

# now we can set priors for the intercepts and slopes
alpha <- normal(0, 5)
beta_age <- normal(0, 5)
beta_mannf <- normal(0, 5)
beta_mannf2 <- normal(0, 5)
beta_cvflow <- normal(0, 5)
gamma_mannf <- normal(0, sigma_mannf, dim = max(age_predictor))
gamma_mannf2 <- normal(0, sigma_mannf2, dim = max(age_predictor))
gamma_cvflow <- normal(0, sigma_cvflow, dim = max(age_predictor))
gamma_system <- normal(0, sigma_system, dim = max(system_vec))
gamma_cohort <- normal(0, sigma_cohort, dim = max(cohort_vec))

# set up linear predictor
mu <- alpha + beta_age * age_predictor +
  (beta_mannf + gamma_mannf[age_predictor]) * mannf_vec +
  (beta_mannf2 + gamma_mannf2[age_predictor]) * mannf_vec * mannf_vec +
  (beta_cvflow + gamma_cvflow[age_predictor]) * cvflow_vec +
  gamma_cohort[cohort_vec] + gamma_system[system_vec]

# what likelihood do we want?
distribution(response_vec) <- poisson(exp(mu))

# compile greta model
mod <- model(alpha, beta_age,
             beta_mannf, beta_cvflow, beta_mannf2,
             gamma_mannf2, gamma_mannf, gamma_cvflow,
             gamma_system, gamma_cohort,
             sigma_cohort, sigma_system, sigma_cvflow, sigma_mannf)

# sample from greta model
draws <- mcmc(mod, n_samples = 200000, warmup = 100000, thin = 50)

# summarise fitted model
mod_summary <- summary(draws)

# pull out parameters and define predictions
alpha_est <- mod_summary$quantiles[grep("alpha", rownames(mod_summary$quantiles)), "50%"]
beta_age_est <- mod_summary$quantiles[grep("beta_age", rownames(mod_summary$quantiles)), "50%"]
beta_mannf_est <- mod_summary$quantiles[grep("beta_mannf$", rownames(mod_summary$quantiles)), ]
beta_mannf2_est <- mod_summary$quantiles[grep("beta_mannf2", rownames(mod_summary$quantiles)), ]
beta_cvflow_est <- mod_summary$quantiles[grep("beta_cvflow", rownames(mod_summary$quantiles)), ]
gamma_mannf_est <- mod_summary$quantiles[grep("gamma_mannf\\[", rownames(mod_summary$quantiles)), ]
gamma_mannf2_est <- mod_summary$quantiles[grep("gamma_mannf2", rownames(mod_summary$quantiles)), ]
gamma_cvflow_est <- mod_summary$quantiles[grep("gamma_cvflow", rownames(mod_summary$quantiles)), ]
gamma_system_est <- mod_summary$quantiles[grep("gamma_system", rownames(mod_summary$quantiles)), "50%"]
gamma_cohort_est <- mod_summary$quantiles[grep("gamma_cohort", rownames(mod_summary$quantiles)), "50%"]

# plot flow effects
mannf_effects <- sweep(gamma_mannf_est, 2, beta_mannf_est, "+")
cvflow_effects <- sweep(gamma_cvflow_est, 2, beta_cvflow_est, "+")
mannf2_effects <- sweep(gamma_mannf2_est, 2, beta_mannf2_est, "+")

# calculate residuals by system and year
n_obs_tmp <- ncol(age_mat)
sys_year_obs <- sys_year_pred <- sys_year_resid <- array(NA, dim = c(nyear, n_obs_tmp, nsystem))
for (sys_set in seq_len(nsystem)) {
  age_sub <- age_mat[system_info == sys_set, ]
  year_sub <- year_info[system_info == sys_set]
  cohort_sub <- cohort_mat[system_info == sys_set, ]
  mannf_sub <- mannf_std[system_info == sys_set]
  cvflow_sub <- cvflow_std[system_info == sys_set]
  if (is.null(nrow(age_sub)))
    age_sub <- matrix(age_sub, nrow = 1)
  
  # create some vectors of predictors
  age_pred <- rep(seq_len(n_obs_tmp), each = nrow(age_sub))
  system_pred <- rep(sys_set, n_obs_tmp * nrow(age_sub))
  cohort_pred <- c(cohort_sub)
  mannf_pred <- rep(mannf_sub, times = n_obs_tmp)
  cvflow_pred <- rep(cvflow_sub, times = n_obs_tmp)
  
  # define median predicted value
  pred_mid <- exp(alpha_est + beta_age_est * age_pred +
                    (beta_mannf_est[3] + gamma_mannf_est[age_pred, 3]) * mannf_pred +
                    (beta_mannf2_est[3] + gamma_mannf2_est[age_pred, 3]) * mannf_pred * mannf_pred +
                    (beta_cvflow_est[3] + gamma_cvflow_est[age_pred, 3]) * cvflow_pred +
                    gamma_system_est[system_pred] + gamma_cohort_est[cohort_pred])
  
  sys_year_obs[year_sub, , sys_set] <- age_sub
  sys_year_pred[year_sub, , sys_set] <- pred_mid
  sys_year_resid[year_sub, , sys_set] <- age_sub - pred_mid
}
system_names <- c("Broken", "Campaspe", "Goulburn",
                  "King", "Loddon", "Lower Murray",
                  "Ovens", "Pyramid Creek")
par(mar = c(4.5, 4.5, 3.1, 1.1))
for (i in seq_len(nsystem)) {
  
  if (system_names[i] != "Pyramid Creek") {

    # pdf(file = paste0("outputs/age_class_strength_", system_names[i], ".pdf"),
    #     width = 7, height = 7)
    par(mfrow = c(3, 1))
    
    # plot observed values
    abs_range <- max(abs(range(sys_year_obs[, , i], na.rm = TRUE)))
    image_breaks <- c(seq(0, abs_range, length = 100), abs_range + 0.0001)
    col_pal_image <- colorRampPalette(c("#f7f7f7",
                                        "#2166ac", "#053061"))(99)
    col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
    sys_year_plot <- ifelse(is.na(sys_year_obs[, , i]), abs_range + 0.00005, sys_year_obs[, , i])
    fields::image.plot(sys_year_plot,
                       col = col_pal_image,
                       breaks = image_breaks,
                       xaxt = "n", yaxt = "n",
                       xlab = "Year", ylab = "Age")
    axis(1, at = seq(0.05, 1, length = 10), labels = seq(2000, 2018, by = 2))
    axis(2, at = seq(0, 1, length = n_obs_tmp), labels = seq(0, n_obs_tmp - 1, by = 1),
         las = 1)
    mtext("Observed abundances", side = 3, line = 1, adj = 1, cex = 1.1)
    mtext(system_names[i], side = 3, line = 1, adj = 0, cex = 1.35)
    
    # plot fitted values
    abs_range <- max(abs(range(sys_year_pred[, , i], na.rm = TRUE)))
    image_breaks <- c(seq(0, abs_range, length = 100), abs_range + 0.0001)
    col_pal_image <- colorRampPalette(c("#f7f7f7",
                                        "#2166ac", "#053061"))(99)
    col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
    sys_year_plot <- ifelse(is.na(sys_year_pred[, , i]), abs_range + 0.00005, sys_year_pred[, , i])
    fields::image.plot(sys_year_plot,
                       col = col_pal_image,
                       breaks = image_breaks,
                       xaxt = "n", yaxt = "n",
                       xlab = "Year", ylab = "Age")
    axis(1, at = seq(0.05, 1, length = 10), labels = seq(2000, 2018, by = 2))
    axis(2, at = seq(0, 1, length = n_obs_tmp), labels = seq(0, n_obs_tmp - 1, by = 1),
         las = 1)
    mtext("Modelled abundances", side = 3, line = 1, adj = 1, cex = 1.1)
    
    # plot residuals
    abs_range <- max(abs(range(sys_year_resid[, , i], na.rm = TRUE)))
    image_breaks <- c(seq(-abs_range, abs_range, length = 100), abs_range + 0.0001)
    col_pal_image <- colorRampPalette(c("#67001f", "#b2182b",
                                        "#f7f7f7",
                                        "#2166ac", "#053061"))(99)
    col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
    sys_year_plot <- ifelse(is.na(sys_year_resid[, , i]), abs_range + 0.00005, sys_year_resid[, , i])
    fields::image.plot(sys_year_plot,
                       col = col_pal_image,
                       breaks = image_breaks,
                       xaxt = "n", yaxt = "n",
                       xlab = "Year", ylab = "Age")
    axis(1, at = seq(0.05, 1, length = 10), labels = seq(2000, 2018, by = 2))
    axis(2, at = seq(0, 1, length = n_obs_tmp), labels = seq(0, n_obs_tmp - 1, by = 1),
         las = 1)
    mtext("Residual (age class strength)", side = 3, line = 1, adj = 1, cex = 1.1)
    
    # dev.off()
    
  }
  
}
