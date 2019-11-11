# load packages
library(rstanarm)

# load some helper functions
source("code/plot-helpers.R")

# load fitted models and data
mod <- readRDS("outputs/fitted/full-model.rds")
mod_noflow <- readRDS("outputs/fitted/noflow-model.rds")
data_set <- readRDS("outputs/fitted/data-set.rds")
additional_data <- readRDS("outputs/fitted/additional-data.rds")

# pull out fitted values
fitted_full <- posterior_predict(mod, offset = mod$offset)
fitted_noflow <- posterior_predict(mod_noflow)

cor(apply(fitted_full, 2, median), data_set$response_vec)
cor(apply(fitted_noflow, 2, median), data_set$response_vec)

# unpack additional data
flow_scales <- additional_data$flow_scales
flow_expanded <- additional_data$flow_expanded
age_mat <- additional_data$age_mat
system_info <- additional_data$system
nsystem <- length(unique(system_info))
year_info <- additional_data$year
nyear <- length(unique(year_info))

# plot some flow effects
system_names <- c("Broken", "Goulburn",
                  "King", "Murray", "Ovens")


pred_list <- c("rrang_spwn_mld",
               "prop_spr_lt_win",
               "prop_sum_lt_win",
               "maxan_mld",
               "spwntmp_c",
               "adult_cpue")

var_names <- c(
  "rrang_spwn_mld" = "Proportional variability in spawning flows",
  "prop_spr_lt_win" = "Spring flow as a proportion of long-term median flow",
  "prop_sum_lt_win" = "Summer flow as a proportion of long-term median flow",
  "maxan_mld" = "Maximum annual daily flow (ML/day)",
  "maxan_mld_ym1" = "Maximum annual daily flow in previous year (ML/day)",
  "spwntmp_c" = "Temperature during spawning period (degrees C)",
  "adult_cpue" = "Adult CPUE"
)
var_names <- var_names[pred_list]

for (i in seq_along(pred_list)) {  
  file_name <- paste0("outputs/plots/rstanarm-", pred_list[i],".jpg")
  # jpeg(file = file_name, height = 7, width = 7, units = "in", res = 300)
  plot_associations(mod,
                    variable = pred_list[i],
                    other_predictors = pred_list[!pred_list %in% pred_list[i]],
                    rescale = flow_scales,
                    xlab = var_names[i])
  # dev.off()
}

# calculate residuals by system and year
n_obs_tmp <- ncol(age_mat)
sys_year_obs <- sys_year_pred <- sys_year_resid <- array(NA, dim = c(nyear, n_obs_tmp, nsystem))
for (sys_set in seq_len(nsystem)) {
  age_sub <- age_mat[system_info == sys_set, ]
  year_sub <- year_info[system_info == sys_set]
  effort_expanded <- data_set$sampling_effort[seq_along(system_info)][system_info == sys_set]
  system_expanded <- rep(sys_set, n_obs_tmp * nrow(age_sub))
  year_expanded <- rep(year_sub, times = n_obs_tmp)
  effort_expanded <- rep(effort_expanded, times = n_obs_tmp)
  flow_sub <- flow_expanded[which(system_info == sys_set) + rep(seq(0, n_obs_tmp * (length(system_info) - 1), by = length(system_info)), each = nrow(age_sub)), ]
  flow_sub <- ifelse(flow_sub, 0, 0)
  if (is.null(nrow(age_sub)))
    age_sub <- matrix(age_sub, nrow = 1)
  pred_mid <- posterior_predict(
    mod_noflow, 
    newdata = data.frame(age_predictor = rep(seq_len(n_obs_tmp), each = nrow(age_sub)),
                         system_vec = system_expanded,
                         sampling_effort = effort_expanded)
  )
  pred_mid <- matrix(apply(pred_mid, 2, median), ncol = n_obs_tmp)
  sys_year_obs[year_sub - min(year_info) + 1, , sys_set] <- age_sub
  sys_year_pred[year_sub - min(year_info) + 1, , sys_set] <- pred_mid
  sys_year_resid[year_sub - min(year_info) + 1, , sys_set] <- age_sub - pred_mid
}
for (i in seq_len(nsystem)) {
  
  pdf(file = paste0("outputs/plots/age_class_strength_", system_names[i], ".pdf"),
      width = 8, height = 5)
  par(mfrow = c(2, 1))
  par(mar = c(4.1, 4.5, 2.5, 1.1))
  
  # plot observed values
  abs_range <- max(abs(range(sys_year_obs[, , i], na.rm = TRUE)))
  image_breaks <- c(seq(0, abs_range, length = 100), abs_range + 0.0001)
  col_pal_image <- colorRampPalette(c("#f7f7f7", "#2166ac", "#053061"))(99)
  col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
  sys_year_plot <- ifelse(is.na(sys_year_obs[, , i]), abs_range + 0.00005, sys_year_obs[, , i])
  fields::image.plot(sys_year_plot,
                     col = col_pal_image,
                     breaks = image_breaks,
                     xaxt = "n", yaxt = "n",
                     xlab = "Year", ylab = "Age")
  axis(1, at = seq(0, 1, length = 21), labels = seq(1999, 2019, by = 1))
  axis(2, at = seq(0, 1, length = n_obs_tmp), labels = seq(0, n_obs_tmp - 1, by = 1),
       las = 1)
  mtext("Observed abundances", side = 3, line = 0.5, adj = 1, cex = 1.1)
  mtext(system_names[i], side = 3, line = 0.5, adj = 0, cex = 1.35)
  
  # plot residuals
  abs_range <- max(abs(range(sys_year_resid[, , i], na.rm = TRUE)))
  image_breaks <- c(seq(-abs_range, abs_range, length = 100), abs_range + 0.0001)
  col_pal_image <- colorRampPalette(c("#67001f", "#b2182b", "#f7f7f7", "#2166ac", "#053061"))(99)
  col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
  sys_year_plot <- ifelse(is.na(sys_year_resid[, , i]), abs_range + 0.00005, sys_year_resid[, , i])
  fields::image.plot(sys_year_plot,
                     col = col_pal_image,
                     breaks = image_breaks,
                     xaxt = "n", yaxt = "n",
                     xlab = "Year", ylab = "Age")
  axis(1, at = seq(0, 1, length = 21), labels = seq(1999, 2019, by = 1))
  axis(2, at = seq(0, 1, length = n_obs_tmp), labels = seq(0, n_obs_tmp - 1, by = 1),
       las = 1)
  mtext("Residual (age class strength)", side = 3, line = 0.5, adj = 1, cex = 1.1)
  
  dev.off()
  
}

# plot variance components
PPD <- posterior_predict(mod)
vars <- apply(PPD, MARGIN = 1, FUN = var)
PPD_sys <- posterior_predict(mod, re.form = ~ (rrang_vec + rrang_ym1_vec +
                                                 psprw_vec + psprw_ym1_vec +
                                                 psumw_vec + psumw_ym1_vec + 
                                                 (psumw_vec ^ 2) + 
                                                 maxann_vec + maxann_ym1_vec +
                                                 spwntmp_vec | system_vec))
vars_sys <- apply(PPD_sys, MARGIN = 1, FUN = var)
PPD_age <- posterior_predict(mod, re.form = ~ (-1 +
                                                 rrang_vec + rrang_ym1_vec +
                                                 psprw_vec + psprw_ym1_vec +
                                                 psumw_vec + psumw_ym1_vec +
                                                 (psumw_vec ^ 2) + 
                                                 maxann_vec + maxann_ym1_vec +
                                                 spwntmp_vec | age_factor))
vars_age <- apply(PPD_age, MARGIN = 1, FUN = var)
PPD_year <- posterior_predict(mod, re.form = ~ (1 | year_vec))
vars_year <- apply(PPD_year, MARGIN = 1, FUN = var)
PPD_cohort <- posterior_predict(mod, re.form = ~ (1 | cohort_vec))
vars_cohort <- apply(PPD_cohort, MARGIN = 1, FUN = var)
PPD_0 <- posterior_predict(mod, re.form = ~ 0)
vars_0 <- apply(PPD_0, MARGIN = 1, FUN = var)
 
# pull out variance explained by each
boxplot(vars_0, vars_year, vars_cohort, vars_sys, vars_age, vars, log = "y",
        names = c("Base", "Year", "Cohort", "System", "Age", "All"))

## MAKES SENSE: LOOKS AT AMOUNT OF VARIANCE ADDED TO PREDICTED VALUES BY EACH RE.
