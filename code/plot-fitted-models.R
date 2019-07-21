# load some helper functions
source("code/plot-helpers.R")

# load fitted models and data
mod <- readRDS("outputs/fitted/full-model.rds")
data_set <- readRDS("outputs/fitted/data-set.rds")
additional_data <- readRDS("outputs/fitted/additional-data.rds")

# unpack additional data
flow_scales <- additional_data$flow_scales
age_mat <- additional_data$age_mat
nyear <- additional_data$nyear
nsystem <- additional_data$nsystem
system_info <- additional_data$system_info
year_info <- additional_data$year_info
cohort_mat <- additional_data$cohort_mat

# plot some flow effects
system_names <- c("Broken", "Goulburn",
                  "King", "Murray", "Ovens")

for (i in seq_along(system_names)) {
  
  # spring flows relative to long-term winter median
  pdf(file = paste0("outputs/plots/spring-flow-effects-", system_names[i],".pdf"), height = 8, width = 6)
  par(mfrow = c(2, 2))
  plot_associations(mod, variable = "psprw_vec", data = data_set,
                    rescale = flow_scales, xlab = "Spring flow proportional to long-term winter average",
                    system = i, cohort = data_set$cohort_vec[data_set$system_vec == i][1])
  dev.off()
  
  # previous spring flows relative to long-term winter median
  pdf(file = paste0("outputs/plots/spring-ym1-flow-effects-", system_names[i],".pdf"), height = 8, width = 6)
  par(mfrow = c(2, 2))
  plot_associations(mod, variable = "psprw_ym1_vec", data = data_set,
                    rescale = flow_scales, xlab = "Spring flow proportional to long-term winter average",
                    system = i, cohort = data_set$cohort_vec[data_set$system_vec == i][1])
  dev.off()

  # max annual flows
  pdf(file = paste0("outputs/plots/max-flow-effects-", system_names[i],".pdf"), height = 8, width = 6)
  par(mfrow = c(2, 2))
  plot_associations(mod, variable = "maxann_vec", data = data_set,
                    rescale = flow_scales, xlab = "Maximum daily flow (z score)",
                    system = i, cohort = data_set$cohort_vec[data_set$system_vec == i][1])
  dev.off()
  
  # previous max annual flows
  pdf(file = paste0("outputs/plots/max-ym1-flow-effects-", system_names[i],".pdf"), height = 8, width = 6)
  par(mfrow = c(2, 2))
  plot_associations(mod, variable = "maxann_ym1_vec", data = data_set,
                    rescale = flow_scales, xlab = "Maximum antecedent daily flow (z score)",
                    system = i, cohort = data_set$cohort_vec[data_set$system_vec == i][1])
  dev.off()
  
  # summer flows relative to long-term winter median
  pdf(file = paste0("outputs/plots/summer-flow-effects-", system_names[i],".pdf"), height = 8, width = 6)
  par(mfrow = c(2, 2))
  plot_associations(mod, variable = "psumw_vec", data = data_set,
                    rescale = flow_scales, xlab = "Summer flow proportional to long-term winter average",
                    system = i, cohort = data_set$cohort_vec[data_set$system_vec == i][1])
  dev.off()
  
  # previous summer flows relative to long-term winter median
  pdf(file = paste0("outputs/plots/summer-ym1-flow-effects-", system_names[i],".pdf"), height = 8, width = 6)
  par(mfrow = c(2, 2))
  plot_associations(mod, variable = "psumw_ym1_vec", data = data_set,
                    rescale = flow_scales, xlab = "Summer flow proportional to long-term winter average",
                    system = i, cohort = data_set$cohort_vec[data_set$system_vec == i][1])
  dev.off()
  
  # temperature effects
  pdf(file = paste0("outputs/plots/spawning-temp-effects-", system_names[i],".pdf"), height = 8, width = 6)
  par(mfrow = c(2, 2))
  plot_associations(mod, variable = "spwntmp_vec", data = data_set,
                    rescale = flow_scales, xlab = "Water temperature during spawning (C)",
                    system = i, cohort = data_set$cohort_vec[data_set$system_vec == i][1])
  dev.off()
  
  # variability in spawning flows
  pdf(file = paste0("outputs/plots/variability-spawning-flow-effects-", system_names[i],".pdf"), height = 8, width = 6)
  par(mfrow = c(2, 2))
  plot_associations(mod, variable = "rrang_vec", data = data_set,
                    rescale = flow_scales, xlab = "Maximum three-day change in flow (ML)",
                    system = i, cohort = data_set$cohort_vec[data_set$system_vec == i][1])
  dev.off()
  
  # variability in previous spawning flows
  pdf(file = paste0("outputs/plots/variability-ym1-spawning-flow-effects-", system_names[i],".pdf"), height = 8, width = 6)
  par(mfrow = c(2, 2))
  plot_associations(mod, variable = "rrang_ym1_vec", data = data_set,
                    rescale = flow_scales, xlab = "Maximum three-day change in flow (ML)",
                    system = i, cohort = data_set$cohort_vec[data_set$system_vec == i][1])
  dev.off()
  
}

# calculate residuals by system and year
n_obs_tmp <- ncol(age_mat)
sys_year_obs <- sys_year_pred <- sys_year_resid <- array(NA, dim = c(nyear, n_obs_tmp, nsystem))
for (sys_set in seq_len(nsystem)) {
  age_sub <- age_mat[system_info == sys_set, ]
  year_sub <- year_info[system_info == sys_set]
  cohort_sub <- cohort_mat[system_info == sys_set, ]
  rrang_sub <- rrang_std[system_info == sys_set]
  rrang_ym1_sub <- rrang_ym1_std[system_info == sys_set]
  psprw_sub <- psprw_std[system_info == sys_set]
  psprw_ym1_sub <- psprw_ym1_std[system_info == sys_set]
  psumw_sub <- psumw_std[system_info == sys_set]
  psumw_ym1_sub <- psumw_ym1_std[system_info == sys_set]
  maxann_sub <- maxann_std[system_info == sys_set]
  maxann_ym1_sub <- maxann_ym1_std[system_info == sys_set]
  spwntmp_sub <- spwntmp_std[system_info == sys_set]
  spwntmp_sub[is.na(spwntmp_sub)] <- 0
  if (is.null(nrow(age_sub)))
    age_sub <- matrix(age_sub, nrow = 1)
  pred_mid <- posterior_predict(mod, newdata = data.frame(age_predictor = rep(seq_len(n_obs_tmp), each = nrow(age_sub)),
                                                          age_factor = factor(rep(seq_len(n_obs_tmp), each = nrow(age_sub))),
                                                          system_vec = rep(sys_set, n_obs_tmp * nrow(age_sub)),
                                                          cohort_vec = c(cohort_sub),
                                                          year_vec = c(year_sub),
                                                          rrang_vec = rep(rrang_sub, times = n_obs_tmp),
                                                          rrang_ym1_vec = rep(rrang_ym1_sub, times = n_obs_tmp),
                                                          psprw_vec = rep(psprw_sub, times = n_obs_tmp),
                                                          psprw_ym1_vec = rep(psprw_ym1_sub, times = n_obs_tmp),
                                                          psumw_vec = rep(psumw_sub, times = n_obs_tmp),
                                                          psumw_ym1_vec = rep(psumw_ym1_sub, times = n_obs_tmp),
                                                          maxann_vec = rep(maxann_sub, times = n_obs_tmp),
                                                          maxann_ym1_vec = rep(maxann_ym1_sub, times = n_obs_tmp),
                                                          spwntmp_vec = rep(spwntmp_sub, times = n_obs_tmp)))
  pred_mid <- matrix(apply(pred_mid, 2, median), ncol = n_obs_tmp)
  sys_year_obs[year_sub, , sys_set] <- age_sub
  sys_year_pred[year_sub, , sys_set] <- pred_mid
  sys_year_resid[year_sub, , sys_set] <- age_sub - pred_mid
}
par(mar = c(4.5, 4.5, 3.1, 1.1))
for (i in seq_len(nsystem)) {
  
  pdf(file = paste0("outputs/plots/age_class_strength_", system_names[i], ".pdf"),
      width = 7, height = 7)
  par(mfrow = c(3, 1))
  
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
  axis(1, at = seq(0.05, 1, length = 10), labels = seq(2000, 2018, by = 2))
  axis(2, at = seq(0, 1, length = n_obs_tmp), labels = seq(0, n_obs_tmp - 1, by = 1),
       las = 1)
  mtext("Observed abundances", side = 3, line = 1, adj = 1, cex = 1.1)
  mtext(system_names[i], side = 3, line = 1, adj = 0, cex = 1.35)
  
  # plot fitted values
  abs_range <- max(abs(range(sys_year_pred[, , i], na.rm = TRUE)))
  image_breaks <- c(seq(0, abs_range, length = 100), abs_range + 0.0001)
  col_pal_image <- colorRampPalette(c("#f7f7f7", "#2166ac", "#053061"))(99)
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
  col_pal_image <- colorRampPalette(c("#67001f", "#b2182b", "#f7f7f7", "#2166ac", "#053061"))(99)
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
