# set working directory
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(greta)

# source helper functions
source("code/greta-helpers.R")

# need to load all otolith data
oti_data <- readRDS("data/otolith-data-loaded.rds")

# filter to MC
oti_data <- oti_data[oti_data$SPECIES == "Maccullochella peelii", ]

# define model
oti_analysis_data <- data.frame(length = oti_data$T_Length..mm.,
                                age = oti_data$AGE,
                                site = as.integer(as.factor(oti_data$Site)),
                                year = as.integer(as.factor(oti_data$Year)))
oti_analysis_data <- oti_analysis_data[apply(oti_analysis_data, 1, function(x) !any(is.na(x))), ]
oti_analysis_data$site <- as.integer(as.factor(oti_analysis_data$site))
oti_analysis_data$year <- as.integer(as.factor(oti_analysis_data$year))

# length to age (vals for CT MC mod)
age_alpha <- normal(0, 10)
age_beta <- normal(0, 10)

# define a function to convert size to age based on model in Todd & Koehn 2008
inverse_growth <- function(x, int, slope) {
  
  int + slope * x
  
}

# random intercepts
# sigma_site_oti <- normal(0, sd(oti_analysis_data$age), truncation = c(0, Inf))
# sigma_year_oti <- normal(0, sd(oti_analysis_data$age), truncation = c(0, Inf))
# gamma_site_oti <- normal(0, sigma_site_oti, dim = length(unique(oti_analysis_data$site)))
# gamma_year_oti <- normal(0, sigma_year_oti, dim = length(unique(oti_analysis_data$year)))

# linear predictor
length_scaled <- scale(oti_analysis_data$length)
age_est_log <- inverse_growth(length_scaled,
                          age_alpha, age_beta)

# add likelihood for age model
sigma_oti <- normal(0, 10, truncation = c(0, Inf))
log_age_p1 <- log(oti_analysis_data$age + 1)
distribution(log_age_p1) <- normal(age_est_log, sigma_oti, truncation = c(0, Inf))

# compile model
mod <- model(age_alpha, age_beta)

# sample
draws <- mcmc(mod, n_samples = 5000, warmup = 5000)

# make up a prediction function from draws
age_pred <- function(alpha_samples, beta_samples, lengths, n = 1000) {
  
  idx <- sample(seq_len(length(beta_samples)), size = n, replace = TRUE)
  
  out <- sweep(beta_samples[idx] %o% lengths, 1, alpha_samples[idx], "+")
  
  exp(out) - 1
  
}

alpha_samples <- c(sapply(draws, function(x) x[, "age_alpha"]))
beta_samples <- c(sapply(draws, function(x) x[, "age_beta"]))

lengths <- seq(0, 1000, by = 50)
lengths <- (lengths - attributes(length_scaled)$`scaled:center`) / 
  attributes(length_scaled)$`scaled:scale`
  
to_plot <- age_pred(alpha_samples, beta_samples, lengths)
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
rownames(alk) <- round(lengths * attributes(length_scaled)$`scaled:scale` + 
  attributes(length_scaled)$`scaled:center`)
colnames(alk) <- c(0:max(to_plot))

FSA::alkPlot(alk, type = "area", pal = "grey")


# define a function to convert size to age based on model in Todd & Koehn 2008
inverse_growth2 <- function(x, length_inf, time_zero, k_param, c_param) {
  
  par1 <- 1 / k_param
  ratio_param <- 1 - (time_zero / length_inf)
  par2 <- 1 / (1 - c_param * ratio_param)
  par3 <- (length_inf - time_zero) / (length_inf - x)
  par4 <- c_param * ratio_param
  
  par1 * log(par2 * (par3 - par4))
  
}
age_vec2 <- inverse_growth2(oti_analysis_data$length / 10,
                           150, 6, 0.0011, -103)

plot(age_vec2 ~ oti_analysis_data$length, type = "n", bty = "l", xlab = "Length", ylab = "Age")
for(i in c(0, seq_len(30)))
  lines(c(0, 2000), c(i, i), lty = 1, lwd = 1, col = ggplot2::alpha("grey50", 0.5))
lines(age_vec2[order(oti_analysis_data$length)] ~ sort(oti_analysis_data$length), lwd = 2, col = "black")
