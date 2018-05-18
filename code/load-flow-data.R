# R code to load flow data for each river reach

# load predictor (flow) files into a list from all files with "to_use" suffix
pred_list <- dir("./data/flow-data/")[grep("to_use", dir("./data/flow-data/"))]
predictors <- vector("list", length = length(pred_list))
for (i in seq_along(pred_list)) 
  predictors[[i]] <- read.csv(paste0("./data/flow-data/", pred_list[i]))
names(predictors) <- sapply(strsplit(pred_list, "_"), function(x) paste(x[1], x[2], sep = "_"))

# define some helper functions
max_fun <- function(x) {
  
  out <- NA
  
  if (any(!is.na(x)))
    out <- max(x, na.rm = TRUE)
  
  out
  
}

# format dates and pull out month/year
mean_flow_stats <- vector("list", length = length(predictors))
for (i in seq_along(predictors)) {
  
  predictors[[i]]$date_formatted <- parse_date_time(predictors[[i]]$date,
                                                    orders = c("dmy_HM", "dmy_HMS", "dmy"))
  predictors[[i]]$date_formatted <- format(predictors[[i]]$date_formatted, format = "%d-%m-%Y")
  predictors[[i]]$day <- as.numeric(sapply(strsplit(predictors[[i]]$date_formatted, "-"),
                                           function(x) x[1]))
  predictors[[i]]$month <- as.numeric(sapply(strsplit(predictors[[i]]$date_formatted, "-"),
                                             function(x) x[2]))
  predictors[[i]]$year <- as.numeric(sapply(strsplit(predictors[[i]]$date_formatted, "-"),
                                            function(x) x[3]))
  
  if (any(predictors[[i]]$year > 2017)) {
    predictors[[i]]$year[which(predictors[[i]]$year > 2017)] <- predictors[[i]]$year[which(predictors[[i]]$year > 2017)] - 100
  }
  
  if (length(predictors[[i]]$discharge_ml_d.1)) {
    if (any(is.na(predictors[[i]]$discharge_ml_d))) {
      subset <- which(is.na(predictors[[i]]$discharge_ml_d))
      predictors[[i]]$discharge_ml_d[subset] <- predictors[[i]]$discharge_ml_d.1[subset]
    }
  }

  predictors[[i]]$season <- ifelse(predictors[[i]]$month < 3, "summer",
                                   ifelse(predictors[[i]]$month < 6, "autumn",
                                          ifelse(predictors[[i]]$month < 9, "winter",
                                                 ifelse(predictors[[i]]$month < 12, "spring",
                                                        "summer"))))
  
  # calculate flow characteristics
  mean_flow_stats[[i]] <- data.frame(year = unique(predictors[[i]]$year))
  mean_flow_stats[[i]]$mean_ann_flow <- tapply(predictors[[i]]$discharge_ml_d, predictors[[i]]$year, mean, na.rm = TRUE)
  mean_seasonal_flow <- tapply(predictors[[i]]$discharge_ml_d,
                               list(predictors[[i]]$season,
                                    predictors[[i]]$year),
                               mean, na.rm = TRUE)
  spawning_flow_data <- predictors[[i]][predictors[[i]]$month %in% c(10, 11, 12), ]
  mean_flow_stats[[i]]$mean_spr_flow <- mean_seasonal_flow["spring", ]
  mean_flow_stats[[i]]$mean_sum_flow <- mean_seasonal_flow["summer", ]
  mean_flow_stats[[i]]$max_ann_flow <- tapply(predictors[[i]]$discharge_ml_d, predictors[[i]]$year, max_fun)
  mean_flow_stats[[i]]$cov_ann_flow <- tapply(predictors[[i]]$discharge_ml_d, predictors[[i]]$year,
                                              function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)))
  mean_spwn_flow_tmp <- tapply(spawning_flow_data$discharge_ml_d,
                               spawning_flow_data$year,
                               mean, na.rm = TRUE)
  cov_spwn_flow_tmp <- tapply(spawning_flow_data$discharge_ml_d,
                              spawning_flow_data$year,
                              function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)))
  mean_flow_stats[[i]]$mean_spwn_flow <- rep(NA, nrow(mean_flow_stats[[i]]))
  mean_flow_stats[[i]]$cov_spwn_flow <- rep(NA, nrow(mean_flow_stats[[i]]))
  mean_flow_stats[[i]]$mean_spwn_flow[match(names(mean_spwn_flow_tmp),
                                            mean_flow_stats[[i]]$year)] <- mean_spwn_flow_tmp
  mean_flow_stats[[i]]$cov_spwn_flow[match(names(cov_spwn_flow_tmp),
                                           mean_flow_stats[[i]]$year)] <- cov_spwn_flow_tmp
  
} 

sys_name <- sapply(strsplit(names(predictors), "_"), function(x) toupper(x[1]))
reach_no <- c(3, 5, 5, 2, 3, 1, 2, 3, 4, 5, 3, 4, 1, 2, 1, 3, 5)

predictors <- lapply(predictors, function(x) x[which(x$year > 1999), ])
mean_flow_stats <- lapply(mean_flow_stats, function(x) x[which(x$year > 1999), ])
sys_tmp <- rep(sys_name, times = sapply(mean_flow_stats, nrow))
reach_tmp <- rep(reach_no, times = sapply(mean_flow_stats, nrow))
mean_flow_stats <- do.call("rbind", mean_flow_stats)
mean_flow_stats$system <- sys_tmp
mean_flow_stats$reach <- reach_tmp

mean_flow_stats$system <- gsub("MURRAY", "LOWERMURRAY", mean_flow_stats$system)

# tidy workspace
rm(i, pred_list, mean_seasonal_flow, sys_tmp, reach_tmp)
