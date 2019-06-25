# histogram function to return counts only
hist_fn <- function(x, breaks) {
  hist(x, breaks = breaks, plot = FALSE)$counts
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

