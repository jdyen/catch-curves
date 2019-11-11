# set working directory
setwd("/homevol/jdyen/ccr-models/")

# load some helper functions
source("code/helpers.R")
source("code/fit_ccr.R")
source("code/methods.R")
source("code/validate_ccr.R")

# load compiled survey data
mod <- readRDS("outputs/fitted/mod_test_20190928_1636.rds")

# capitalise system names
systems_to_keep <- c("broken", "goulburn", "king", "murray", "ovens")
system_names <- systems_to_keep
for (i in seq_along(system_names)) {
  init <- toupper(substr(system_names[i], 1, 1))
  remainder <- substr(system_names[i], 2, nchar(system_names[i]))
  system_names[i] <- paste0(init, remainder)
}

# combine all chains into one matrix
samples <- do.call(rbind, lapply(mod, function(x) do.call(rbind, x$draws)))

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

# plot age-class strengths
effort_seq <- rep(NA, length(sys_ids))
pred_set <- matrix(NA, nrow = length(sys_ids), ncol = ncol(mod$data$predictors))
for (i in seq_along(unique(sys_ids))) {
  idx <- mod$data$system == i
  idy <- sys_ids == i
  effort_seq[idy][match(mod$data$year[idx], year_ids[idy])] <- mod$data$effort[idx]
  pred_set[idy, ][match(mod$data$year[idx], year_ids[idy]), ] <- mod$data$predictors[idx, ]
}

effort_seq <- ifelse(is.na(effort_seq), 1000, effort_seq)

effort_seq <- tapply(mod$data$effort, mod$data$system, median)[sys_ids]

# set some means for unobserved cases
pred_set <- ifelse(is.na(pred_set), 0, pred_set)

# create a data set over which to predict
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
length_to_age <- get_param(samples, "length_to_age")
length_to_age <- apply(length_to_age, 2, median)
length_to_age <- matrix(length_to_age, nrow = ncol(observed_all))
length_to_age <- sweep(length_to_age, 2, colSums(length_to_age), "/")
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
