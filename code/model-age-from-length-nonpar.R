# set working directory
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(greta)

# source helper functions
source("code/greta-helpers.R")

# load compiled survey data
oti_data <- readRDS("data/otolith-data-loaded.rds")

# define a function to convert size to age based on model in Todd & Koehn 2008
inverse_growth <- function(x, length_inf, time_zero, k_param, c_param) {
  
  par1 <- 1 / k_param
  ratio_param <- 1 - (time_zero / length_inf)
  par2 <- 1 / (1 - c_param * ratio_param)
  par3 <- (length_inf - time_zero) / (length_inf - x)
  par4 <- c_param * ratio_param
  
  par1 * log(par2 * (par3 - par4))
  
}

# define model
oti_analysis_data <- data.frame(length = oti_data$T_Length..mm.,
                                age = oti_data$AGE,
                                site = as.integer(as.factor(oti_data$Site)),
                                year = as.integer(as.factor(oti_data$Year)))
oti_analysis_data <- oti_analysis_data[apply(oti_analysis_data, 1, function(x) !any(is.na(x))), ]
oti_analysis_data$site <- as.integer(as.factor(oti_analysis_data$site))
oti_analysis_data$year <- as.integer(as.factor(oti_analysis_data$year))

# convert observed lengths to ages
len_par <- normal(150, 10)
time_par <- normal(6, 1, truncation = c(0, Inf))
k_par <- normal(0.0011, 0.001, truncation = c(0, Inf))
c_par <- normal(-103, 10)
age_vec <- inverse_growth(oti_analysis_data$length / 10,
                          len_par, time_par, k_par, c_par)

# define likelihood
sigma_main <- lognormal(0, 0.5)
distribution(oti_analysis_data$age) <- normal(age_vec, sigma_main, truncation = c(0, Inf))

# compile model
mod <- model(len_par, time_par, k_par, c_par)

# sample
init_set <- initials(len_par = 150, time_par = 6, k_par = 0.0011, c_par = -103)
draws <- mcmc(mod, n_samples = 5000, warmup = 5000, initial_values = init_set)

# make up a prediction function from draws
age_pred <- function(c_samples, t_samples, k_samples, len_samples, lengths, n = 1000) {

  idx <- sample(seq_len(length(beta_samples)), size = n, replace = TRUE)

  out <- matrix(NA, nrow = n, ncol = length(lengths))
  for (i in seq_along(lengths))
    out[, i] <- inverse_growth(lengths[i], len_samples[idx], t_samples[idx], k_samples[idx], c_samples[idx])
  
  out

}

c_samples <- c(sapply(draws, function(x) x[, "c_par"]))
t_samples <- c(sapply(draws, function(x) x[, "time_par"]))
k_samples <- c(sapply(draws, function(x) x[, "k_par"]))
len_samples <- c(sapply(draws, function(x) x[, "len_par"]))

lengths <- seq(0, 1000, by = 50)

to_plot <- age_pred(c_samples, t_samples, k_samples, len_samples, lengths / 10)
to_plot[to_plot < 0] <- 0
to_plot <- round(to_plot)

# proportions function
prop_equal <- function(x, range) {

  out <- rep(NA, length(range))
  for (i in seq_along(range))
    out[i] <- sum(x == range[i]) / length(x)

  out

}

alk <- t(apply(to_plot, 2, prop_equal, c(0:max(to_plot))))
rownames(alk) <- lengths
colnames(alk) <- c(0:max(to_plot))

FSA::alkPlot(alk, type = "area", pal = "grey")

plot(oti_analysis_data$age,
     inverse_growth(oti_analysis_data$length / 10,
                    mean(len_samples), mean(t_samples),
                    mean(k_samples), mean(c_samples)),
     bty = "l", pch = 16,
     xlab = "True age", ylab = "Modelled age", col = "gray50", las = 1)

for(i in c(0, seq_len(30)))
  lines(c(0, 2000), c(i, i), lty = 1, lwd = 1, col = ggplot2::alpha("grey50", 0.5))
