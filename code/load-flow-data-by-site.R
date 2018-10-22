# R code to load flow data for each river reach

# load predictor (flow) files into a list from all files with 'to_use' suffix
pred_list <- dir('./data/flow-data/')[grep('to_use', dir('./data/flow-data/'))]
predictors <- vector('list', length = length(pred_list))
date_parse <- c(rep('dmy_HM', 3), 'dmy', rep('dmy_HM', 4), 'dmY', rep('dmy_HM', 8))
for (i in seq_along(pred_list)) {
  predictors[[i]] <- read.csv(paste0('./data/flow-data/', pred_list[i]),
                              stringsAsFactors = FALSE)
  if (i == 9) {
    predictors[[i]]$day <- day(parse_date_time(predictors[[i]]$date,
                                               orders = c('dmy_HM')))
    predictors[[i]]$day[c(21389:nrow(predictors[[i]]))] <-
      as.numeric(predictors[[i]]$date[c(21389:nrow(predictors[[i]]))])
    predictors[[i]]$date <- paste(predictors[[i]]$day,
                                  predictors[[i]]$month,
                                  predictors[[i]]$year,
                                  sep = '/')
  }
  predictors[[i]]$Event_Date <- parse_date_time2(predictors[[i]]$date,
                                                orders = date_parse[i],
                                                cutoff_2000 = 18)
  predictors[[i]] <- predictors[[i]][year(predictors[[i]]$Event_Date) > 1990, ]
  predictors[[i]] <- predictors[[i]][!is.na(predictors[[i]]$date), ]
}
names(predictors) <- sapply(strsplit(pred_list, '_'), function(x) paste(x[1], x[2], sep = '_'))

# calculate flow metrics
calc_flow_metrics <- function(x, na_thresh = 0.2) {
  
  if (is.null(x)) {
    out <- rep(NA, 12)
  } else {
    
    depth_data <- x[, grep('water_level', colnames(x))]
    if (is.matrix(depth_data) | is.data.frame(depth_data)) {
      if (ncol(depth_data) > 1)
        depth_data <- depth_data[, -grep('qc', colnames(depth_data))]
    }
    if (is.null(depth_data))
      stop('something wrong with depth colnames', call. = FALSE)
    
    flow_data <- x[, grep('discharge_ml', colnames(x))]
    if (is.matrix(flow_data) | is.data.frame(flow_data)) {
      if (ncol(flow_data) > 1)
        flow_data <- flow_data[, -grep('qc', colnames(flow_data))]
    }
    if (is.null(flow_data))
      stop('something wrong with discharge colnames', call. = FALSE)

    if (any(is.na(flow_data))) {
      if ((sum(is.na(flow_data)) / length(flow_data)) > na_thresh) {
        maf <- NA
        covaf <- NA
        maxan <- NA
        minan <- NA
        mspr <- NA
        msum <- NA
        mspwn <- NA
        covsp <- NA
      } else {
        maf <- mean(flow_data, na.rm = TRUE)
        covaf <- sd(flow_data, na.rm = TRUE) / maf
        maxan <- max(flow_data, na.rm = TRUE)
        minan <- min(flow_data, na.rm = TRUE)
        mspr <- mean(flow_data[month(x$Event_Date) %in% c(9:11)], na.rm = TRUE)
        msum <- mean(flow_data[month(x$Event_Date) %in% c(12, 1:2)], na.rm = TRUE)
        mspwn <- mean(flow_data[month(x$Event_Date) %in% c(10:12)], na.rm = TRUE)
        covsp <- sd(flow_data[month(x$Event_Date) %in% c(10:12)], na.rm = TRUE) / mspwn
      }
    } else {
      maf <- mean(flow_data, na.rm = TRUE)
      covaf <- sd(flow_data, na.rm = TRUE) / maf
      maxan <- max(flow_data, na.rm = TRUE)
      minan <- min(flow_data, na.rm = TRUE)
      mspr <- mean(flow_data[month(x$Event_Date) %in% c(9:11)], na.rm = TRUE)
      msum <- mean(flow_data[month(x$Event_Date) %in% c(12, 1:2)], na.rm = TRUE)
      mspwn <- mean(flow_data[month(x$Event_Date) %in% c(10:12)], na.rm = TRUE)
      covsp <- sd(flow_data[month(x$Event_Date) %in% c(10:12)], na.rm = TRUE) / mspwn
    }
    
    if (any(is.na(depth_data))) {
      if ((sum(is.na(depth_data)) / length(depth_data)) > na_thresh) {
        madp <- NA
        cvdp <- NA
        maxdp <- NA
        mindp <- NA
      } else {
        madp <- mean(depth_data, na.rm = TRUE)
        cvdp <- sd(depth_data, na.rm = TRUE) / madp
        maxdp <- max(depth_data, na.rm = TRUE)
        mindp <- min(depth_data, na.rm = TRUE)
      }
    } else {
      madp <- mean(depth_data, na.rm = TRUE)
      cvdp <- sd(depth_data, na.rm = TRUE) / madp
      maxdp <- max(depth_data, na.rm = TRUE)
      mindp <- min(depth_data, na.rm = TRUE)
    }
    
    out <- c(maf, mspr, msum, covaf,
             mspwn, covsp, maxan, minan,
             madp, cvdp, maxdp, mindp)

  }
  
  names(out) <- c('mannf_mld', 'msprf_mld', 'msumf_mld', 'covaf_mld',
                  'mspwn_mld', 'covsp_mld', 'maxan_mld', 'minan_mld',
                  'madpth_m', 'cvdpth_m', 'maxdpth_m', 'mindpth_m')
  
  out
  
}

# pull out preceding 12 months of data
flow_tm1 <- function(x, flow, na_thresh = 0.2) {
  
  systmp <- paste0(tolower(x$SYSTEM), '_r', x$Reach)
  
  if (systmp %in% c('campaspe_r3', 'campaspe_r4'))
    systmp <- 'campaspe_r34'
  if (systmp %in% c('ovens_r1'))
    systmp <- 'ovens_wangaratta'
  if (systmp %in% c('lowermurray_r1'))
    systmp <- 'murray_yarrawonga'
  if (systmp %in% c('ovens_r1'))
    systmp <- 'ovens_wangaratta'
  if (systmp %in% c('loddon_r5'))
    systmp <- 'loddon_r4'
  if (systmp %in% c('broken_r4'))
    systmp <- 'broken_r3'
  if (systmp %in% c('goulburn_r1'))
    systmp <- 'goulburn_r4'
  
  if (any(names(flow) == systmp)) {
    flow_sub <- flow[names(flow) == systmp][[1]]
    tmp <- flow_sub[flow_sub$Event_Date < x$Event_Date, ]
    tmp <- tail(tmp, 365)
  } else {
    tmp <- NULL
  }
  
  calc_flow_metrics(tmp, na_thresh = na_thresh)

}

calc_flow_fn <- function(i, data, predictors, na_thresh = 0.2) {
  flow_tm1(data[i, ], predictors, na_thresh)
}

future::plan(future::multicore)
out_test <- future.apply::future_sapply(seq_len(nrow(alldat)),
                                        calc_flow_fn,
                                        data = alldat,
                                        predictors = predictors,
                                        na_thresh = 0.2)
flow_data <- t(out_test)

# clean up workspace
rm(i, out_test, date_parse, predictors, pred_list)
