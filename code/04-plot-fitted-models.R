# plot associations and outputs from full fitted model

# load some helper functions
source("code/helpers.R")
source("code/plot_functions.R")
source("code/fit_ccr.R")

# load fitted models
mod <- readRDS("outputs/fitted/full-model.rds")
flow_scales <- readRDS("data/flow-standardisation.rds")

# append a model ID to plotted outputs
mod_name <- ""

# capitalise system names
systems_to_keep <- c("broken", "goulburn", "king", "murray", "ovens")
system_order <- c(4, 5, 2, 3, 1)
colour_palette_order <- c(5, 3, 4, 1, 2)
system_names <- systems_to_keep
for (i in seq_along(system_names)) {
  init <- toupper(substr(system_names[i], 1, 1))
  remainder <- substr(system_names[i], 2, nchar(system_names[i]))
  system_names[i] <- paste0(init, remainder)
}

# which systems do we want to plot?
sys_to_plot <- c(1:5)
year_to_plot <- 13
vars_to_plot <- colnames(mod$data$predictors)
var_names <- c("spawning_variability" = "Proportional variability in spawning flows",
               "prop_spring_lt" = "Spring flow relative to long-term flow",
               "prop_summer_lt" = "Summer flow relative to long-term flow",
               "prop_winter_lt" = "Winter flow relative to long-term flow",
               "prop_max_antecedent_lt" = "Antecedent maximum flow relative to long-term flow",
               "spawning_temp" = "Temperature during spawning months",
               "adult_cpue" = "Adult CPUE")
var_names <- var_names[vars_to_plot]

# predict to new data and plot
jpeg(file = paste0("outputs/plots/recruitment_flow_effects", mod_name, ".jpg"),
     width = 7.5, height = 6.8, units = "in", res = 150)
par(mfrow = c(4, 2), mar = c(3.8, 4.1, 1.1, 1.1))
col_pal <- viridis::inferno(256)[seq(1, 250, length = max(sys_to_plot))]
effort_set <- 1000
for (i in seq_along(vars_to_plot)) {
  
  test_var <- vars_to_plot[i]
  pred_data <- create_newdata(mod, var = test_var,
                              system = sys_to_plot, year = year_to_plot,
                              effort = effort_set)
  preds <- predict(mod, newdata = pred_data, thin = 100, lengths = FALSE)
  pred_sum <- apply(preds, c(2, 3), quantile, p = c(0.025, 0.5, 0.975, 0.85))
  
  rm(preds)
  
  for (j in seq_along(unique(sys_to_plot))) {
    
    idx <- pred_data$system == j
    
    xplot <- pred_data$predictors[idx, test_var]
    xplot <- xplot * flow_scales["sd", test_var] + flow_scales["mean", test_var]
    
    if (j == 1) {
      plot(exp(pred_sum[2, idx, 1]) ~ xplot,
           pch = 16, bty = "l", las = 1, type = "n",
           xlab = "", ylab = "",
           xlim = range(pred_data$predictors[, test_var] * flow_scales["sd", test_var] + flow_scales["mean", test_var]),
           ylim = exp(range(pred_sum[, , 1])))
      mtext(var_names[i], side = 1, line = 2.5, adj = 0.5, cex = 0.8)
      if (i %in% c(3)) {
        mtext(paste0("Number of recruits with ", effort_set, " EF seconds"),
              side = 2, line = 2.5, adj = 1, cex = 0.9)
      }
      polygon(exp(c(pred_sum[1, idx, 1], rev(pred_sum[3, idx, 1]))) ~
                c(xplot, rev(xplot)),
              col = scales::alpha(col_pal[colour_palette_order[j]], 0.25),
              border = NA)
      lines(exp(pred_sum[2, idx, 1]) ~ xplot,
            col = col_pal[colour_palette_order[j]], lwd = 3)
    } else {
      polygon(exp(c(pred_sum[1, idx, 1], rev(pred_sum[3, idx, 1]))) ~
                c(xplot, rev(xplot)),
              col = scales::alpha(col_pal[colour_palette_order[j]], 0.25),
              border = NA)
      lines(exp(pred_sum[2, idx, 1]) ~ xplot,
            col = col_pal[colour_palette_order[j]], lwd = 3)
    }
    
  }
  
}
plot(exp(pred_sum[2, idx, 1]) ~ xplot,
     type = "n",
     bty = "n",
     xlab = "", ylab = "",
     xlim = range(pred_data$predictors[, test_var] * flow_scales["sd", test_var] + flow_scales["mean", test_var]),
     ylim = exp(range(pred_sum[, , 1])),
     xaxt = "n",
     yaxt = "n")
legend(legend = system_names[system_order],
       x = "center",
       cex = 1.25,
       lty = 1, lwd = 2, bty = "o",
       col = col_pal,
       xpd = TRUE)

dev.off()

# pull out Bayesian P-values (P(effect != 0))
var_names2 <- c("spawning_variability" = "Proportional variability in spawning flows",
                "prop_spring_lt" = "Spring flow relative to long-term flow",
                "prop_summer_lt" = "Summer flow relative to long-term flow",
                "prop_winter_lt" = "Winter flow relative to long-term flow",
                "prop_max_antecedent_lt" = "Antecedent maximum flow relative to long-term flow",
                "spawning_temp" = "Temperature during spawning months",
                "adult_cpue" = "Adult CPUE")

pr_pos <- function(x) 
  sum(x > 0) / length(x)
boxplot_fn <- function(x, pvals = NULL, xlab = NULL, col = "black",
                       cex.lab = 0.8, cex.axis = 0.9) {
  
  vals <- apply(x, 2, quantile, p = c(0.025, 0.1, 0.5, 0.9, 0.975))
  nx <- ncol(vals)

  xplot <- seq_len(nx)  
  xrange <- c(0.5, nx + 0.5)
  yrange <- range(vals)
  yrange[1] <- yrange[1] - 0.3 * abs(yrange[1])
  yrange[1] <- ifelse(yrange[2] > 0.75 & abs(yrange[1]) < 1, -abs(yrange[2]), yrange[1])
  yrange[1] <- ifelse(abs(yrange[1]) < 1, -1, yrange[1])
  
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
    text(xplot, rep(yrange[1] + 0.09 * abs(yrange[1]), nx), format(round(pvals, 2), digits = 2), cex = 1.2, adj = 0.5)
  
}

pdf(file = paste0("outputs/plots/flow_effects_by_system", mod_name, ".pdf"), height = 8, width = 7)
par(mar = c(5.1, 5.1, 2.1, 1.1), mfrow = c(4, 2))

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
  beta_sub <- beta_flow[, idx]
  p_vals[i, ] <- apply(beta_sub, 2, pr_pos)
  boxplot_fn(beta_sub[, system_order], pvals = p_vals[i, system_order], xlab = system_names[system_order])
  mtext(var_names2[i], side = 3, adj = 0, line = 0.5, cex = 0.75)
}
dev.off()

# plot YCS
range_expand <- function(x) (min(x):max(x) - 6)
year_ranges <- tapply(mod$data$year, mod$data$system, range_expand)
sys_ids <- rep(1:5, times = sapply(year_ranges, length))
year_ids <- unlist(year_ranges)

# plot age-class strengths
n_age <- ncol(mod$data$length_age_matrix)
age_sets <- rep(seq_len(n_age), each = length(year_ids))
idx <- match(paste(sys_ids, year_ids, sep = "_"),
             paste(mod$data$system, mod$data$year - 6, sep = "_"))
pred_means <- apply(mod$data$predictors, 2, function(x) tapply(x[1:67], mod$data$system, mean))
pred_means_expanded <- pred_means[sys_ids, ]
pred_set <- do.call(rbind, lapply(seq_len(n_age), function(x) pred_means_expanded))
effort_seq <- tapply(mod$data$effort, mod$data$system, median)[sys_ids]

# create a data set over which to predict
test_data <-  list(length_age_matrix = mod$data$length_age_matrix,
                   system = sys_ids,
                   year = year_ids,
                   predictors = pred_set,
                   effort = effort_seq)
max_age <- 3
fitted_all <- predict(mod, test_data, lengths = FALSE, survey = TRUE, thin = 10)
fitted_all <- exp(apply(fitted_all, c(2, 3), median))
fitted_all <- fitted_all[, seq_len(max_age + 1)]
sys_sub <- test_data$system
year_sub <- test_data$year
observed_all <- mod$data$response
length_to_age <- get_param(samples, "length_to_age")
length_to_age <- apply(length_to_age, 2, median)
length_to_age <- matrix(length_to_age, nrow = ncol(observed_all))
length_to_age <- sweep(length_to_age, 2, colSums(length_to_age), "/")
observed_all <- observed_all %*% length_to_age
observed_all <- observed_all[, seq_len(max_age + 1)]
observed_tmp <- matrix(NA, nrow = length(sys_sub), ncol = ncol(observed_all))

scaled_abundances <- rep(1, 5)

filename <- paste0("outputs/plots/observed-all-systems", mod_name, ".jpg")
jpeg(file = filename, width = 7.5, height = 6.8, units = "in", res = 300)
par(mfrow = c(5, 1), mar = c(4, 3.8, 2.2, 1.1))

sys_abunds <- tapply(apply(mod$data$response, 1, sum),
                     mod$data$system,
                     sum)

unique_systems <- unique(mod$data$system)
for (i in unique_systems) {
  
  observed <- observed_tmp[sys_sub == i, ]
  observed[match(mod$data$year[mod$data$system == i] - 6, year_sub[sys_sub == i]), ] <-
    observed_all[mod$data$system == i, ]
  observed <- scaled_abundances[i] * observed
  
  idx <- sys_sub == i
  idy <- order(year_sub[sys_sub == i])
  
  fitted <- fitted_all[idx, ]
  fitted <- fitted[idy, ]
  residuals <- observed - fitted
  
  # plot observed
  plot_values <- matrix(NA, nrow = max(year_sub), ncol = ncol(fitted))
  plot_values[match(year_sub[sys_sub == i], seq_len(max(year_sub))), ] <- observed
  abs_range <- max(abs(range(observed, na.rm = TRUE)))
  image_breaks <- c(seq(0, abs_range, length = 100), abs_range + 0.0001)
  col_pal_image <- colorRampPalette(c("#f7f7f7",
                                      "#2166ac", "#053061"))(99)
  col_pal_image <- c(col_pal_image, ggplot2::alpha("black", 0.7))
  plot_values <- ifelse(is.na(plot_values), abs_range + 0.00005, plot_values)
  fields::image.plot(plot_values,
                     col = col_pal_image,
                     breaks = image_breaks,
                     xaxt = "n", yaxt = "n",
                     xlab = "", ylab = "")
  bin_width <- 1 / nrow(plot_values)
  axis(1, at = seq(bin_width, 1 - bin_width, length = 10), labels = seq(2000, 2018, by = 2))
  axis(2, at = seq(0, 1, length = (max_age + 1)), labels = seq(0, max_age, by = 1),
       las = 1)
  mtext(paste0(system_names[i], " (n = ", sys_abunds[i], ")"), side = 3, line = 0.5, adj = 1, cex = 1.1)
  mtext("Year", side = 1, line = 2.5, adj = 0.5, cex = 1)
  mtext("Age", side = 2, line = 2.5, adj = 0.5, cex = 1)
  
}

dev.off()

# plot year class strength
filename <- paste0("outputs/plots/ycs-all-systems", mod_name, ".jpg")
jpeg(file = filename, width = 7.5, height = 6.8, units = "in", res = 300)
par(mfrow = c(5, 1), mar = c(4, 3.8, 2.2, 1.1))

# use observed flow conditions, padded to system means in unsampled years
pred_set <- mod$data$predictors[1:67, ]
sys_year_match <- paste(mod$data$system, mod$data$year, sep = "_")
sys_year_plot <- paste(sys_ids, year_ids + 6, sep = "_")
pred_set <- pred_set[match(sys_year_plot, sys_year_match), ]
flow_sys <- apply(mod$data$predictors[1:67, ], 2, function(x) tapply(x, mod$data$system, mean))
pred_set[is.na(pred_set[, 1]), ] <- flow_sys[sys_ids[is.na(pred_set[, 1])], ]
test_data <-  list(length_age_matrix = mod$data$length_age_matrix,
                   system = sys_ids,
                   year = year_ids,
                   predictors = pred_set,
                   effort = rep(1000, nrow(pred_set)))
fitted_all <- predict(mod, test_data, lengths = FALSE, survey = FALSE, thin = 10)
fitted_upper_all <- exp(apply(fitted_all, c(2, 3), quantile, p = 0.4))
fitted_lower_all <- exp(apply(fitted_all, c(2, 3), quantile, p = 0.6))
fitted_all <- exp(apply(fitted_all, c(2, 3), median))

unique_systems <- unique(mod$data$system)
unique_systems <- c(4, 5, 2, 3, 1)
for (i in unique_systems) {

  idx <- sys_sub == i
  idy <- order(year_sub[sys_sub == i])

  fitted <- fitted_all[idx, ]
  fitted <- fitted[idy, ]
  fitted_upper <- fitted_upper_all[idx, ]
  fitted_upper <- fitted_upper[idy, ]
  fitted_lower <- fitted_lower_all[idx, ]
  fitted_lower <- fitted_lower[idy, ]
  
  # plot modelled
  plot_values <- plot_upper <- plot_lower <- matrix(NA, nrow = max(year_sub), ncol = ncol(fitted))
  plot_values[match(year_sub[sys_sub == i], seq_len(max(year_sub))), ] <- fitted
  plot_upper[match(year_sub[sys_sub == i], seq_len(max(year_sub))), ] <- fitted_upper
  plot_lower[match(year_sub[sys_sub == i], seq_len(max(year_sub))), ] <- fitted_lower
  xplot <- seq_len(nrow(plot_values))
  plot(plot_values[, 1] ~ xplot,
       type = "n",
       xlab = "", ylab = "", 
       bty = "l",
       las = 1,
       xaxt = "n",
       ylim = range(c(0, plot_values[, 1], plot_upper[, 1], plot_lower[, 1]), na.rm = TRUE))
  to_keep <- !is.na(plot_upper[, 1])
  polygon(c(xplot[to_keep], rev(xplot[to_keep])),
          c(plot_upper[to_keep, 1], rev(plot_lower[to_keep, 1])),
          col = scales::alpha("black", 0.25),
          border = NA)
  lines(plot_values[, 1] ~ xplot, col = "black", lwd = 2)
  axis(1, at = seq(2, nrow(plot_values), by = 2), labels = seq(2000, 2018, by = 2))
  mtext(system_names[i], side = 3, line = 0.5, adj = 1, cex = 1.1)
  if (i == 1)
    mtext("Year", side = 1, line = 2.5, adj = 0.5, cex = 1)
  if (i == 2)
    mtext("Number of recruits per 1000 EF seconds", side = 2, line = 2.5, adj = 0.5, cex = 1)

}

dev.off()
