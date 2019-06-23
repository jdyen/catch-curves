# something
plot_associations <- function(mod,
                              variable, data, rescale = NULL,
                              nplot = 100,
                              system = 5, cohort = 1,  year = 15,
                              labels = NULL,
                              xlab = "Flow (ML / day)", ylab = "Abundance",
                              ...) {

  data_tmp <- data.frame(rrang_vec = rep(0, nplot),
                         rrang_ym1_vec = rep(0, nplot),
                         psprw_vec = rep(0, nplot),
                         psprw_ym1_vec = rep(0, nplot),
                         psumw_vec = rep(0, nplot),
                         psumw_ym1_vec = rep(0, nplot),
                         minwin_vec = rep(0, nplot),
                         spwntmp_vec = rep(0, nplot),
                         system_vec = factor(rep(system, nplot)),
                         year_vec = rep(year, nplot),
                         cohort_vec = rep(cohort, nplot))
  
  # set chosen variable as a sequence rather than holding it at its mean value
  data_tmp[variable] <- seq(min(data[variable][data$system_vec == system, 1], na.rm = TRUE),
                            max(data[variable][data$system_vec == system, 1], na.rm = TRUE),
                            length = nplot)
  
  # set some plot labels
  if (is.null(labels))
    label_set <- c("Young-of-year", paste0(seq_len(max(data$age_predictor)), " year olds"))

  # calculate fitted catch curves for each age class
  flow_pred <- vector("list", length = max(data$age_predictor))
  for (i in seq_len(max(data$age_predictor))) {
    age_set <- rep(i, nplot)
    data_tmp$age_predictor <- age_set
    data_tmp$age_factor <- factor(age_set)
    flow_pred[[i]] <- posterior_predict(mod, newdata = data_tmp)
  } 
  
  # summarise the predicted catch curves as mean and 95% credible intervals
  flow_mean <- lapply(flow_pred, function(x) apply(x, 2, median, na.rm = TRUE))
  flow_lower <- lapply(flow_pred, function(x) apply(x, 2, quantile, p = 0.025, na.rm = TRUE))
  flow_upper <- lapply(flow_pred, function(x) apply(x, 2, quantile, p = 0.975, na.rm = TRUE))

  # set up a plotting x variable
  x_set <- seq(min(data_tmp[variable]), max(data_tmp[variable]), length = nplot)
  if (!is.null(rescale))
    x_set <- x_set * rescale[[variable]]$sd + rescale[[variable]]$mean
  
  # plot it
  for (i in seq_along(flow_mean)) {
    plot(flow_mean[[i]] ~ x_set,
         type = "n", las = 1, bty = "l",
         xlab = xlab, ylab = ylab,
         lwd = 2,
         ylim = range(c(flow_mean[[i]], flow_lower[[i]], flow_upper[[i]])),
         ...)
    polygon(c(x_set, rev(x_set)),
            c(flow_lower[[i]], rev(flow_upper[[i]])),
            col = "gray50", border = NA)
    lines(flow_mean[[i]] ~ x_set, lwd = 2, col = "black")
    mtext(label_set[i], side = 3, adj = 1, line = 1, cex = 1)
  } 
  
}
