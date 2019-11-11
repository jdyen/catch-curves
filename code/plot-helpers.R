plot_associations <- function(
  mod,
  variable,
  other_predictors,
  rescale = NULL,
  nplot = 100,
  year = 18,
  effort = 1000,
  systems = c(4, 5, 2, 3, 1),
  xlab = "Predictor",
  ylab = "Abundance",
  col_pal = viridis::inferno(256),
  ...
) {
  
  # pull out some indices
  n_age <- max(mod$data$age_predictor)
  
  # create a dummy matrix to fill with appropriate values to plot
  mat_tmp <- matrix(0, nrow = nplot, ncol = ncol(mod$data))
  colnames(mat_tmp) <- colnames(mod$data)
  data_tmp <- as.data.frame(mat_tmp)
  data_tmp$year_vec <- rep(year, nplot)
  data_tmp$sampling_effort <- rep(log(effort), nplot)
  
  # subset colour palette
  col_pal <- col_pal[floor(seq(1, 250, length = length(systems)))]
  
  # add the age class we want (YOY in this model)
  data_tmp$age_predictor <- rep(1, nplot)
  data_tmp$age_factor <- factor(data_tmp$age_predictor)
  
  # expand over all systems
  n_system <- length(unique(mod$data$system_vec))
  data_tmp <- data_tmp[rep(seq_len(nplot), times = n_system), ]
  data_tmp$system_vec <- rep(systems, each = nplot)
    
  # set chosen variable as a sequence rather than holding it at its mean value
  for (i in seq_along(systems)) {
    idx <- data_tmp$system_vec == systems[i]
    data_tmp[variable][idx, ] <- seq(min(mod$data[variable][mod$data$system_vec == systems[i], 1], na.rm = TRUE),
                                     max(mod$data[variable][mod$data$system_vec == systems[i], 1], na.rm = TRUE),
                                     length = nplot)
  }
  
  # set other variables to their system means to avoid weird intercepts
  for (i in seq_along(other_predictors)) {
    var_tmp <- other_predictors[i]
    sys_means <- tapply(mod$data[var_tmp][, 1], mod$data$system_vec, mean)
    data_tmp[var_tmp] <- sys_means[data_tmp$system_vec]
  }
  
  # calculate posterior predictions
  # flow_pred <- predict(mod_freq, newdata = data_tmp, offset = log(data_tmp$sampling_effort), type = "response")
  flow_pred <- posterior_predict(mod, newdata = data_tmp,
                                 offset = log(data_tmp$sampling_effort))

  # summarise the predicted catch curves as mean and 95% credible intervals
  # flow_mean <- flow_pred
  # flow_upper <- flow_lower <- flow_pred
  flow_mean <- apply(flow_pred, 2, mean, na.rm = TRUE)
  flow_sd <- apply(flow_pred, 2, sd, na.rm = TRUE)
  flow_lower <- flow_mean - flow_sd
  flow_lower <- ifelse(flow_lower < 0, 0, flow_lower)
  flow_upper <- flow_mean + flow_sd
  # flow_lower <- apply(flow_pred, 2, quantile, p = 0.025, na.rm = TRUE)
  # flow_upper <- apply(flow_pred, 2, quantile, p = 0.975, na.rm = TRUE)

  # plot it
  idx <- data_tmp$system_vec == systems[1]
  
  # set up a plotting x variable
  x_set <- seq(min(data_tmp[variable][idx, ]), max(data_tmp[variable][idx, ]), length = nplot)
  if (!is.null(rescale))
    x_set <- x_set * flow_scales["sd", variable] + flow_scales["mean", variable]
  
  plot(flow_mean[idx] ~ x_set,
       type = "n", las = 1, bty = "l",
       xlab = xlab, ylab = ylab,
       lwd = 2,
       xlim = c(range(data_tmp[variable]) * flow_scales["sd", variable] + flow_scales["mean", variable]),
       ylim = range(c(flow_mean, flow_lower, flow_upper)),
       ...)
  polygon(c(x_set, rev(x_set)),
          c(flow_lower[idx], rev(flow_upper[idx])),
          col = scales::alpha(col_pal[1], 0.25), border = NA)
  lines(flow_mean[idx] ~ x_set, lwd = 2, col = col_pal[1])
  
  for (i in seq_along(systems)[-1]) {
    
    idx <- data_tmp$system_vec == systems[i]

    # set up a plotting x variable
    x_set <- seq(min(data_tmp[variable][idx, ]), max(data_tmp[variable][idx, ]), length = nplot)
    if (!is.null(rescale))
      x_set <- x_set * flow_scales["sd", variable] + flow_scales["mean", variable]
    
    polygon(c(x_set, rev(x_set)),
            c(flow_lower[idx], rev(flow_upper[idx])),
            col = scales::alpha(col_pal[i], 0.25), border = NA)
    lines(flow_mean[idx] ~ x_set, lwd = 2, col = col_pal[i])
  }
  
  sys_names <- c("Broken", "Goulburn", "King", "Murray", "Ovens")
  legend(legend = sys_names[systems],
         x = "topright",
         lty = 1, lwd = 2, bty = "o",
         col = col_pal)
  
}
