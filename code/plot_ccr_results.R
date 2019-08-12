setwd("~/Dropbox/research/catch-curves/")

# load fitted models
fitted_mods <- readRDS("outputs/fitted-models.rds")
cv_mods <- readRDS("outputs/validated-models.rds")

#
r2cv <- sapply(cv_mods, function(x) cor(c(x$predicted), c(x$observed)) ** 2)


# capitalise system names
system_names <- systems_to_keep
for (i in seq_along(system_names)) {
  init <- toupper(substr(system_names[i], 1, 1))
  remainder <- substr(system_names[i], 2, nchar(system_names[i]))
  system_names[i] <- paste0(init, remainder)
}

# pull out a single model to plot
mod <- fitted_mods[[1]]

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
               "spwntmp_c" = "Temperature during spawning period (degrees C)")
var_names <- var_names[colnames(mod$data$predictors)]

# predict to new data and plot
col_pal <- viridis::inferno(256)[seq(1, 250, length = max(sys_to_plot))]
effort_set <- 1000
for (i in seq_along(vars_to_plot)) {
  
  jpeg(file = paste0("outputs/plots/recruitment_flow_effects_", vars_to_plot[i], ".jpg"),
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
                "spwntmp_c" = "Temperature during spawning months")
var_names2 <- var_names2[colnames(mod$data$predictors)]

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
samples <- mod$draws
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
predicted_ages <- predict(mod, lengths = FALSE)
ave_abund <- exp(apply(predicted_ages, c(2, 3), mean))
recruit_prop <- apply(ave_abund, 1, function(x) x[1] / sum(x))
rm(predicted_ages)

# boxplot of average numbers of recruits by system
boxplot(ave_abund[, 1] ~ system)

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
max_age <- 3
fitted_all <- predict(mod, lengths = FALSE)
fitted_all <- exp(apply(fitted_all, c(2, 3), median))
fitted_all <- fitted_all[, seq_len(max_age + 1)]
sys_sub <- mod$data$system
year_sub <- mod$data$year
observed_all <- mod$data$response
length_to_age <- get_param(mod$draws, "length_to_age")
length_to_age <- apply(length_to_age, 2, median)
length_to_age <- matrix(length_to_age, nrow = ncol(observed_all))
observed_all <- observed_all %*% length_to_age
observed_all <- observed_all[, seq_len(max_age + 1)]

unique_systems <- unique(mod$data$system)
for (i in seq_along(unique_systems)) {
  
  idx <- sys_sub == i
  idy <- order(year_sub[sys_sub == i])
  
  observed <- observed_all[idx, ]
  observed <- observed[idy, ]
  fitted <- fitted_all[idx, ]
  fitted <- fitted[idy, ]
  residuals <- observed - fitted
  
  par(mfrow = c(3, 1), mar = c(4.5, 4.5, 3.1, 1.1))
  
  # plot observed values
  abs_range <- max(abs(range(observed, na.rm = TRUE)))
  image_breaks <- c(seq(0, abs_range, length = 100), abs_range + 0.0001)
  col_pal_image <- colorRampPalette(c("#f7f7f7",
                                      "#2166ac", "#053061"))(99)
  col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
  plot_values <- ifelse(is.na(observed), abs_range + 0.00005, observed)
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
  abs_range <- max(abs(range(fitted, na.rm = TRUE)))
  image_breaks <- c(seq(0, abs_range, length = 100), abs_range + 0.0001)
  col_pal_image <- colorRampPalette(c("#f7f7f7",
                                      "#2166ac", "#053061"))(99)
  col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
  plot_values <- ifelse(is.na(fitted), abs_range + 0.00005, fitted)
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
  abs_range <- max(abs(range(residuals, na.rm = TRUE)))
  image_breaks <- c(seq(-abs_range, abs_range, length = 100), abs_range + 0.0001)
  col_pal_image <- colorRampPalette(c("#67001f", "#b2182b",
                                      "#f7f7f7",
                                      "#2166ac", "#053061"))(99)
  col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
  plot_values <- ifelse(is.na(residuals), abs_range + 0.00005, residuals)
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
