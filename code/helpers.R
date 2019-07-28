# histogram function to return counts only
hist_fn <- function(x, breaks) {
  hist(x, breaks = breaks, plot = FALSE)$counts
}

# create a pairwise matrix matching classes from two variables
classify <- function(x, y) {
  out <- matrix(0, nrow = max(x, na.rm = TRUE), ncol = max(y, na.rm = TRUE))
  for (i in seq_along(x)) {
    if (!is.na(x[i]) & !is.na(y[i]))
      out[x[i], y[i]] <- out[x[i], y[i]] + 1
  }
  out
}

# define a function to convert size to age based on model in Todd & Koehn 2008
inverse_growth <- function(x, length_inf, time_zero, k_param, c_param) {
  
  par1 <- 1 / k_param
  ratio_param <- 1 - (time_zero / length_inf)
  par2 <- 1 / (1 - c_param * ratio_param)
  par3 <- (length_inf - time_zero) / (length_inf - x)
  par4 <- c_param * ratio_param
  
  par1 * log(par2 * (par3 - par4))
  
}

# rebase a vector intended as an index to start at 1 and count up with no gaps
rebase_index <- function(x)
  as.integer(as.factor(x))

# rescale data and store the standardisations
extract_standards <- function(x) {
  c("mean" = mean(x, na.rm = TRUE), "sd" = sd(x, na.rm = TRUE))
}

# subset a data object from a fitted CCR model
subset_data <- function(x, idx) {
  
  list(response = x$response[idx, ],
       length_age_matrix = x$length_age_matrix,
       system = x$system[idx],
       year = x$year[idx],
       predictors = x$predictors[idx, ],
       cohort_mat = x$cohort_mat[idx, ])
  
}
