setwd("~/Dropbox/research/catch-curves/")

## CONSIDER RE-RUNNING WITH TOTAL ABUNDANCE-BY-YEAR AS A PREDICTOR. CIRCULAR? USE CV TO TEST?
mod <- readRDS("outputs/fitted/mod-long-mcmc.rds")

# load some helper functions
source("code/methods.R")
source("code/fit_ccr.R")

# load fitted models
fitted_mods <- readRDS("outputs/fitted-models.rds")
cv_mods <- readRDS("outputs/validated-models.rds")
flow_scales <- readRDS("data/flow-standardisation.rds")

# can we run a model selection step using cross-validated r2 values?
r2_fn <- function(x, y) {
  
  log_yx <- ifelse(y > 0 & x > 0, log(y / x), 0)
  dev_full <- sum(y * log_yx - (y - x))
  log_ymean <- ifelse(y > 0, log(y / mean(y)), 0)
  dev_null <- sum(y * log_ymean)
  
  1 - (dev_full / dev_null)
  
}

r2cv <- sapply(cv_mods, function(x) r2_fn(c(exp(x$predicted)), c(x$observed)))
r2naive <- sapply(fitted_mods, function(x) cor(exp(c(apply(predict(x, year = TRUE, cohort = TRUE), c(2, 3), median))), c(x$data$response)))

# r2cv <- sapply(cv_mods, function(x) r2_fn(c(exp(x$predicted)), c(x$observed)))
r2naive <- sapply(fitted_mods, function(x) r2_fn(exp(c(apply(predict(x, year = TRUE, cohort = TRUE), c(2, 3), median))), c(x$data$response)))

# capitalise system names
systems_to_keep <- c("broken", "goulburn", "king", "murray", "ovens")
system_names <- systems_to_keep
for (i in seq_along(system_names)) {
  init <- toupper(substr(system_names[i], 1, 1))
  remainder <- substr(system_names[i], 2, nchar(system_names[i]))
  system_names[i] <- paste0(init, remainder)
}

# pull out a single model to plot
mod <- fitted_mods[[3]]

# plot some traces to diagnose poor mixing and convergence
bayesplot::mcmc_trace(mod$draws, regex_pars = "cohort\\[143,")

# which systems do we want to plot?
sys_to_plot <- c(1:5)
year_to_plot <- 10
vars_to_plot <- colnames(mod$data$predictors)
vars_to_plot <- vars_to_plot[grep("_sq", vars_to_plot, invert = TRUE)]
var_names <- c("rrang_spwn_mld" = "Rolling range of spawning flows (ML/day)",
               "prop_spr_lt_win" = "Spring flow as a proportion of long-term median flow",
               "prop_sum_lt_win" = "Summer flow as a proportion of long-term median flow",
               "prop_sum_lt_win_ym1" = "Previous year's summer flow as a proportion of long-term median flow",
               "maxan_mld" = "Maximum annual daily flow (ML/day)",
               "maxan_mld_ym1" = "Maximum annual daily flow in previous year (ML/day)",
               "spwntmp_c" = "Temperature during spawning period (degrees C)",
               "pop_abund" = "Total catch")
var_names <- var_names[vars_to_plot]

# predict to new data and plot
col_pal <- viridis::inferno(256)[seq(1, 250, length = max(sys_to_plot))]
effort_set <- 1000
for (i in seq_along(vars_to_plot)) {
  
  jpeg(file = paste0("outputs/plots/zzz_recruitment_flow_effects_", vars_to_plot[i], ".jpg"),
       width = 7.5, height = 6.8, units = "in", res = 150) 
  # pdf(file = paste0("outputs/plots/recruitment_flow_effects_", vars_to_plot[i], ".pdf"),
  #     width = 7.5, height = 6.8)
  par(mfrow = c(1, 1), mar = c(5.1, 6.1, 2.1, 1.1))
  
  test_var <- vars_to_plot[i]
  pred_data <- create_newdata(mod, var = test_var,
                              system = sys_to_plot, year = year_to_plot,
                              effort = effort_set)
  preds <- predict(mod, newdata = pred_data, cohort = FALSE, lengths = FALSE)
  pred_sum <- apply(preds, c(2, 3), quantile, p = c(0.25, 0.5, 0.75, 0.85))
  
  rm(preds)
  
  for (j in seq_along(unique(sys_to_plot))) {
    
    idx <- pred_data$system == j
    
    xplot <- pred_data$predictors[idx, test_var]
    xplot <- xplot * flow_scales["sd", test_var] + flow_scales["mean", test_var]
    
    if (j == 1) {
      plot(exp(pred_sum[2, idx, 1]) ~ xplot,
           pch = 16, bty = "l", las = 1, type = "l",
           xlab = "", ylab = "",
           col = col_pal[j], ylim = exp(range(pred_sum[, , 1])),
           log = "y", lwd = 3)
      mtext(var_names[i], side = 1, line = 2.5, adj = 0.5)
      mtext(paste0("Estimated abundance with ", effort_set, " EF seconds"),
            side = 2, line = 4.5, adj = 0.5)
      polygon(exp(c(pred_sum[1, idx, 1], rev(pred_sum[3, idx, 1]))) ~
                c(xplot, rev(xplot)),
              col = scales::alpha(col_pal[j], 0.25),
              border = NA)
    } else {
      lines(exp(pred_sum[2, idx, 1]) ~ xplot,
            col = col_pal[j], lwd = 3)
      polygon(exp(c(pred_sum[1, idx, 1], rev(pred_sum[3, idx, 1]))) ~
                c(xplot, rev(xplot)),
              col = scales::alpha(col_pal[j], 0.25),
              border = NA)
    }
    
  }
  
  legend(legend = system_names,
         x = "topright",
         lty = 1, lwd = 2, bty = "o",
         col = col_pal)
  
  dev.off()
  
}

# pull out Bayesian P-values (P(effect != 0))
var_names2 <- c("rrang_spwn_mld" = "Variability in spawning flows",
                "prop_spr_lt_win" = "Spring flow relative to long-term flow",
                "prop_sum_lt_win" = "Summer flow relative to long-term flow",
                "prop_sum_lt_win_ym1" = "Antecedent summer flow relative to long-term flow",
                "prop_sum_win_sq" = "Summer flow quadratic effect",
                "maxan_mld" = "Maximum annual daily flow",
                "maxan_mld_ym1" = "Antecedent maximum annual daily flow",
                "spwntmp_c" = "Temperature during spawning months",
                "pop_abund" = "Total catch")

pr_pos <- function(x) 
  sum(x > 0) / length(x)
boxplot_fn <- function(x, pvals = NULL, xlab = NULL, col = "black",
                       cex.lab = 0.8, cex.axis = 0.9) {
  
  vals <- apply(x, 2, quantile, p = c(0.025, 0.1, 0.5, 0.9, 0.975))
  nx <- ncol(vals)

  xplot <- seq_len(nx)  
  xrange <- c(0.5, nx + 0.5)
  yrange <- range(vals)
  yrange[1] <- yrange[1] - 0.5
  
  plot(vals["50%", ] ~ xplot,
       type = "n",
       bty = "l",
       xaxt = "n", yaxt = "n",
       las = 1,
       xlab = "", ylab = "",
       xlim = xrange, ylim = yrange)
  for (i in xplot) {
    lines(c(i, i), c(vals["2.5%", i], rev(vals["97.5%", i])),
          col = col,
          lty = 1, lwd = 1.8)
    lines(c(i, i), c(vals["10%", i], rev(vals["90%", i])),
          col = col,
          lty = 1, lwd = 3)
  } 
  points(vals["50%", ] ~ xplot,
         pch = 16, cex = 2, col = col)
  lines(c(0, nx + 1), c(0, 0), lty = 2, lwd = 1.5)
    
  mtext("System", side = 1, line = 2.7, adj = 0.5, cex = cex.lab)
  mtext("Standardised coefficient", side = 2, line = 2.7, adj = 0.5, cex = cex.lab)
  axis(2, las = 1, cex = cex.axis)
  if (!is.null(xlab))
    axis(1, at = xplot, labels = xlab, cex = cex.axis)

  if (!is.null(pvals))
    text(xplot, rep(yrange[1] + 0.1, nx), format(round(pvals, 2), digits = 2), cex = 1.2, adj = 0.5)
  
}

pdf(file = "outputs/plots/flow_effects_by_system.pdf", height = 8, width = 11)
par(mar = c(5.1, 5.1, 2.1, 1.1), mfrow = c(3, 3))

vars_to_test <- colnames(mod$data$predictors)
var_names2 <- var_names2[vars_to_test]
samples <- do.call(rbind, mod$draws)
beta_flow <- get_param(samples, "pred_effects")
age_set <- 1
system_set <- 1:5
p_vals <- matrix(NA, nrow = length(vars_to_test), ncol = length(system_set))
for (i in seq_along(vars_to_test)) {
  var_name <- vars_to_test[i]
  idx <- grep(paste0(",", match(var_name, colnames(mod$data$predictors)), "\\]"), colnames(beta_flow))
  n_age <- ncol(mod$data$length_age_matrix)
  idx <- idx[n_age * (system_set - 1) + age_set]
  beta_sub <- beta_flow[, idx]
  p_vals[i, ] <- apply(beta_sub, 2, pr_pos)
  boxplot_fn(beta_sub, pvals = p_vals[i, ], xlab = system_names)
  mtext(var_names2[i], side = 3, adj = 0, line = 0.5, cex = 0.75)
}
dev.off()

# plot proportion of recruits
range_expand <- function(x) min(x):max(x)
year_ranges <- tapply(mod$data$year, mod$data$system, range_expand)
new_cohorts <- matrix(NA, nrow = length(unlist(year_ranges)), ncol = ncol(mod$data$cohort_mat))
for (i in seq_along(year_ranges)) {
  
  idx <- mod$data$system == i
  if (i > 1) {
    idy <- (sum(sapply(year_ranges[seq_len(i - 1)], length)) + 1):sum(sapply(year_ranges[seq_len(i)], length))
  } else {
    idy <- 1:sum(sapply(year_ranges[seq_len(i)], length))
  }
    
  cohort_sub <- mod$data$cohort_mat[idx, ]
  year_sub <- mod$data$year[idx]
  
  missing <- !year_ranges[[i]] %in% year_sub
  missing_yrs <- year_ranges[[i]][missing]
  
  new_cohorts[idy[!missing], ] <- cohort_sub[order(year_sub), ]
  
  for (j in seq_along(missing_yrs)) {
    year_match <- year_ranges[[i]] == missing_yrs[j]
    if ((missing_yrs[j] - 1) %in% year_sub) {
      new_cohorts[idy[year_match], ] <- cohort_sub[year_sub == (missing_yrs[j] - 1), ] + 1
    } else {
      if ((missing_yrs[j] + 1) %in% year_sub) {
        new_cohorts[idy[year_match], ] <- cohort_sub[year_sub == (missing_yrs[j] + 1), ] - 1
      } else {
        new_cohorts[idy[year_match], ] <- rep(NA, ncol(cohort_new))
      }
    }
  }
        
}
sys_ids <- rep(1:5, times = sapply(year_ranges, length))
year_ids <- unlist(year_ranges)

# test that cohorts are all correct
# correct <- rep(NA, 10000)
# for (i in seq_len(10000)) {
#   sysx <- sample(1:5, size = 1)
#   yrx <- sample(mod$data$year[mod$data$system == sysx], size = 1)
#   correct[i] <- all(mod$data$cohort_mat[mod$data$system == sysx & mod$data$year == yrx, ] == 
#     new_cohorts[sys_ids == sysx & year_ids == yrx, ])
# }
# all(correct)

test_data <-  list(length_age_matrix = mod$data$length_age_matrix,
                   system = sys_ids,
                   year = year_ids,
                   predictors = matrix(0, nrow = length(sys_ids), ncol = ncol(mod$data$predictors)),
                   effort = rep(1000, length(sys_ids)),
                   cohort_mat = new_cohorts)
predicted_ages <- predict(mod, newdata = test_data, lengths = FALSE, cohort = TRUE, year = TRUE)
ave_abund <- exp(apply(predicted_ages, c(2, 3), mean))
recruit_prop <- apply(ave_abund, 1, function(x) x[1] / sum(x))
rm(predicted_ages)
boxplot(ave_abund[, 1] ~test_data$system, log = "", names = system_names, las = 1, ylab = "", xlab = "")
mtext("YOY CPUE per 1000 EF seconds", side = 2, line = 3.9, adj = 0.5)
mtext("System", side = 1, line = 2.7, adj = 0.5)
boxplot(ave_abund[, 1] ~test_data$system, log = "y", names = system_names, las = 1, ylab = "", xlab = "")
mtext("YOY CPUE per 1000 EF seconds", side = 2, line = 3.9, adj = 0.5)
mtext("System", side = 1, line = 2.7, adj = 0.5)

# boxplot of average numbers of recruits by system
predicted_ages <- predict(mod, lengths = FALSE, year = TRUE, cohort = TRUE, effort = 1000)
ave_abund <- exp(apply(predicted_ages, c(2, 3), mean))
pdf(file = "outputs/plots/num_recruits_boxplot.pdf", width = 7, height = 7)
par(mfrow = c(2, 1), mar = c(4.1, 5.1, 1.1, 1.1))
boxplot(ave_abund[, 1] ~ mod$data$system, log = "", names = system_names, las = 1, ylab = "", xlab = "")
mtext("YOY CPUE per 1000 EF seconds", side = 2, line = 3.9, adj = 0.5)
mtext("System", side = 1, line = 2.7, adj = 0.5)
boxplot(ave_abund[, 1] ~ mod$data$system, log = "y", names = system_names, las = 1, ylab = "", xlab = "")
mtext("YOY CPUE per 1000 EF seconds", side = 2, line = 3.9, adj = 0.5)
mtext("System", side = 1, line = 2.7, adj = 0.5)
dev.off()

pdf(file = "outputs/plots/proportion_recruits_points.pdf", height = 7, width = 6)
par(mfrow = c(3, 2), mar = c(4.1, 4.1, 2.0, 1.1))
col_pal <- viridis::inferno(256)[seq(1, 250, length = max(sys_to_plot))]
nsystem <- length(unique(mod$data$system))
for (i in seq_len(nsystem)) {
  xplot <- 1999:2019
  yplot <- rep(NA, length(xplot))
  names(yplot) <- xplot
  yplot[as.character(c(mod$data$year[mod$data$system == i] + 1998))] <-
    recruit_prop[mod$data$system == i]
  plot(yplot ~ xplot,
       type = "p",pch = 16, cex = 1.1, #lty = 1, lwd = 2,
       col = "black",
       ylim = c(0, 1),
       bty = "l",
       las = 1,
       xlab = "Year",
       ylab = "Proportion recruits")
  mtext(system_names[i], side = 3, line = 0.5, adj = 0)
}
dev.off()

# plot age-class strengths?
# "observed" is observed abundances (years in rows, age classes in columns)
# "fitted" is modelled abundances (years in rows, age classes in columns)
# "residuals" are the observed minus modelled abundances (years in rows, age classes in columns)
effort_seq <- rep(NA, length(sys_ids))
pred_set <- matrix(NA, nrow = length(sys_ids), ncol = ncol(mod$data$predictors))
for (i in seq_along(unique(sys_ids))) {
  idx <- mod$data$system == i
  idy <- sys_ids == i
  effort_seq[idy][match(mod$data$year[idx], year_ids[idy])] <- mod$data$effort[idx]
  pred_set[idy, ][match(mod$data$year[idx], year_ids[idy]), ] <- mod$data$predictors[idx, ]
}
effort_seq <- ifelse(is.na(effort_seq), 1000, effort_seq)

## MEANS DON'T SEEM TO BE NEAR ZERO. WHY?
pred_set <- ifelse(is.na(pred_set), -0.3, pred_set)


test_data <-  list(length_age_matrix = mod$data$length_age_matrix,
                   system = sys_ids,
                   year = year_ids,
                   predictors = pred_set,
                   effort = effort_seq,
                   cohort_mat = new_cohorts)
max_age <- 3
fitted_all <- predict(mod, test_data, lengths = FALSE, year = TRUE, cohort = TRUE, thin = 20)
fitted_all <- exp(apply(fitted_all, c(2, 3), median))
fitted_all <- fitted_all[, seq_len(max_age + 1)]
sys_sub <- test_data$system
year_sub <- test_data$year
observed_all <- mod$data$response
length_to_age <- get_param(mod$draws, "length_to_age")
length_to_age <- apply(length_to_age, 2, median)
length_to_age <- matrix(length_to_age, nrow = ncol(observed_all))
# length_to_age[1, ] <- c(0.96, 0.04, rep(0, 14))
observed_all <- observed_all %*% length_to_age
observed_all <- observed_all[, seq_len(max_age + 1)]
observed_tmp <- matrix(NA, nrow = length(sys_sub), ncol = ncol(observed_all))

unique_systems <- unique(mod$data$system)
for (i in seq_along(unique_systems)) {
  
  observed <- observed_tmp[sys_sub == i, ]
  observed[match(mod$data$year[mod$data$system == i], year_sub[sys_sub == i]), ] <-
    observed_all[mod$data$system == i, ]
  
  idx <- sys_sub == i
  idy <- order(year_sub[sys_sub == i])
  
  fitted <- fitted_all[idx, ]
  fitted <- fitted[idy, ]
  residuals <- observed - fitted
  
  par(mfrow = c(3, 1), mar = c(4.5, 4.5, 3.1, 1.1))
  
  # plot observed values
  plot_values <- matrix(NA, nrow = max(year_sub), ncol = ncol(fitted))
  plot_values[match(year_sub[sys_sub == i], seq_len(max(year_sub))), ] <- observed
  abs_range <- max(abs(c(range(observed, na.rm = TRUE), range(fitted, na.rm = TRUE))))
  image_breaks <- c(seq(0, abs_range, length = 100), abs_range + 0.0001)
  col_pal_image <- colorRampPalette(c("#f7f7f7",
                                      "#2166ac", "#053061"))(99)
  col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
  plot_values <- ifelse(is.na(plot_values), abs_range + 0.00005, plot_values)
  fields::image.plot(plot_values,
                     col = col_pal_image,
                     breaks = image_breaks,
                     xaxt = "n", yaxt = "n",
                     xlab = "Year", ylab = "Age")
  axis(1, at = seq(0.05, 1, length = 10), labels = seq(2000, 2018, by = 2))
  axis(2, at = seq(0, 1, length = (max_age + 1)), labels = seq(0, max_age, by = 1),
       las = 1)
  mtext("Observed abundances", side = 3, line = 1, adj = 1, cex = 1.1)
  
  # plot fitted values
  plot_values <- matrix(NA, nrow = max(year_sub), ncol = ncol(fitted))
  plot_values[match(year_sub[sys_sub == i], seq_len(max(year_sub))), ] <- fitted
  image_breaks <- c(seq(0, abs_range, length = 100), abs_range + 0.0001)
  col_pal_image <- colorRampPalette(c("#f7f7f7",
                                      "#2166ac", "#053061"))(99)
  col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
  plot_values <- ifelse(is.na(plot_values), abs_range + 0.00005, plot_values)
  fields::image.plot(plot_values,
                     col = col_pal_image,
                     breaks = image_breaks,
                     xaxt = "n", yaxt = "n",
                     xlab = "Year", ylab = "Age")
  axis(1, at = seq(0.05, 1, length = 10), labels = seq(2000, 2018, by = 2))
  axis(2, at = seq(0, 1, length = (max_age + 1)), labels = seq(0, max_age, by = 1),
       las = 1)
  mtext("Modelled abundances", side = 3, line = 1, adj = 1, cex = 1.1)
  
  # plot residuals
  plot_values <- matrix(NA, nrow = max(year_sub), ncol = ncol(fitted))
  plot_values[match(year_sub[sys_sub == i], seq_len(max(year_sub))), ] <- residuals
  abs_range <- max(abs(range(residuals, na.rm = TRUE)))
  image_breaks <- c(seq(-abs_range, abs_range, length = 100), abs_range + 0.0001)
  col_pal_image <- colorRampPalette(c("#67001f", "#b2182b",
                                      "#f7f7f7",
                                      "#2166ac", "#053061"))(99)
  col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
  plot_values <- ifelse(is.na(plot_values), abs_range + 0.00005, plot_values)
  fields::image.plot(plot_values,
                     col = col_pal_image,
                     breaks = image_breaks,
                     xaxt = "n", yaxt = "n",
                     xlab = "Year", ylab = "Age")
  axis(1, at = seq(0.05, 1, length = 10), labels = seq(2000, 2018, by = 2))
  axis(2, at = seq(0, 1, length = (max_age + 1)), labels = seq(0, max_age, by = 1),
       las = 1)
  mtext("Residual (age class strength)", side = 3, line = 1, adj = 1, cex = 1.1)
  
}
