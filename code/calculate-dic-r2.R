# script to calculate DIC and r2 for all fitted models

# load some helpers
source("code/helpers.R")
source("code/methods.R")
source("code/fit_ccr.R")

# load fitted models
to_load <- dir("outputs/fitted/")
to_load <- to_load[grep("mod_test_", to_load)]
to_load <- to_load[grep("mod_test_20191103", to_load, invert = TRUE)]
to_load <- to_load[grep("mod_test_20191102", to_load, invert = TRUE)]
to_load <- to_load[grep("mod_test_20191101", to_load, invert = TRUE)]
mod_list <- lapply(paste0("outputs/fitted/", to_load), readRDS)

# need to pull out fitted values from each
fitted_values <- lapply(mod_list, function(x) exp(predict(x, thin = 12, effort = x$data$effort)))

# then extract log-posterior density of every data point at each iteration (thinned)
extract_lpd <- function(x, y)
  apply(x, 1, function(x) mapply(dpois, y$data$response, x, log = TRUE))
out <- mapply(extract_lpd, fitted_values, mod_list, SIMPLIFY = FALSE)

# calculate fitted means (fitted value at mean parameters)
fitted_mean <- lapply(fitted_values,
                      function(x) apply(x, c(2, 3), mean))
lpd_bar <- apply(
  mapply(function(x, y) mapply(dpois, y$data$response, x, log = TRUE),
         fitted_mean, mod_list),
  2,
  sum
)

# calculate deviances
deviances <- sapply(out, function(x) -2 * apply(x, 2, sum))
dbar <- apply(deviances, 2, mean)
dtbar <- -2 * lpd_bar
pd <- dbar - dtbar
DIC <- pd + dbar

# can estimate r2 values
correlations <- mapply(
  cor,
  lapply(fitted_values, function(x) c(apply(x, c(2, 3), mean))),
  lapply(mod_list, function(x) c(x$data$response))
)
correlations <- correlations ** 2

# make a full table
all_vars <- c("rrang_spwn_mld", "prop_spr_lt", "prop_maxan_lt_ym1", "prop_sum_lt", "prop_win_lt", "spwntmp_c")
var_sets <- c(
  list(all_vars),
  combn(all_vars, m = 5, simplify = FALSE),
  combn(all_vars, m = 4, simplify = FALSE),
  combn(all_vars, m = 3, simplify = FALSE)
)
var_sets <- lapply(var_sets, function(x) c(x, "adult_cpue"))

full_var_names <- var_sets[[1]]
output_table <- do.call(
  rbind,
  lapply(var_sets, function(x, var_names) var_names %in% x, var_names = all_vars)
)
output_table[output_table] <- 1
colnames(output_table) <- all_vars
output_table <- cbind(output_table, "DIC" = DIC, "r2" = correlations)
output_table <- output_table[order(DIC), ]
write.csv(output_table, file = "outputs/tables/model-selection-details.csv", row.names = FALSE)
