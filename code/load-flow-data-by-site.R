# R code to load flow data for each river reach

# load predictor (flow) files into a list from all files with "to_use" suffix
pred_list <- dir("data/flow-data/")[grep("to_use", dir("data/flow-data/"))]
predictors <- vector("list", length = length(pred_list))

# make sure dates parse correctly
date_parse <- c(rep("dmy_HM", 2), "dmy", rep("dmy_HM", 4), "dmY", rep("dmy_HM", 9))

for (i in seq_along(pred_list)) {
  predictors[[i]] <- read.csv(paste0("./data/flow-data/", pred_list[i]),
                              stringsAsFactors = FALSE)
  if (i == 8) {
    predictors[[i]]$day <- day(parse_date_time(predictors[[i]]$date,
                                               orders = c("dmy_HM")))
    predictors[[i]]$day[c(21389:nrow(predictors[[i]]))] <-
      as.numeric(predictors[[i]]$date[c(21389:nrow(predictors[[i]]))])
    predictors[[i]]$date <- paste(predictors[[i]]$day,
                                  predictors[[i]]$month,
                                  predictors[[i]]$year,
                                  sep = "/")
  }
  predictors[[i]]$Event_Date <- parse_date_time2(predictors[[i]]$date,
                                                orders = date_parse[i],
                                                cutoff_2000 = 18)
  predictors[[i]] <- predictors[[i]][year(predictors[[i]]$Event_Date) > 1990, ]
  predictors[[i]] <- predictors[[i]][!is.na(predictors[[i]]$date), ]
}
names(predictors) <- sapply(strsplit(pred_list, "_"), function(x) paste(x[1], x[2], sep = "_"))

# backfill Ovens@Wang with temperature data from logger at Peechelba.
predictors$ovens_wangaratta$water_temp_c <- 
  predictors$ovens_peechelba$water_temp_c[match(predictors$ovens_wangaratta$Event_Date,
                                                predictors$ovens_peechelba$Event_Date)]

# calculate flow metrics
calc_flow_metrics <- function(x, na_thresh = 0.2) {
  
  if (is.null(x)) {
    out <- rep(NA, 20)
  } else {
    
    depth_data <- x[, grep("water_level", colnames(x))]
    if (is.matrix(depth_data) | is.data.frame(depth_data)) {
      if (ncol(depth_data) > 1)
        depth_data <- depth_data[, -grep("qc", colnames(depth_data))]
    }
    if (is.null(depth_data))
      stop("something wrong with depth colnames", call. = FALSE)
    
    flow_data <- x[, grep("discharge_ml", colnames(x))]
    if (is.matrix(flow_data) | is.data.frame(flow_data)) {
      if (ncol(flow_data) > 1)
        flow_data <- flow_data[, -grep("qc", colnames(flow_data))]
    }
    if (is.null(flow_data))
      stop("something wrong with discharge colnames", call. = FALSE)

    temp_data <- x[, grep("temp", colnames(x))]
    if (is.matrix(temp_data) | is.data.frame(temp_data)) {
      if (ncol(temp_data) > 1)
        temp_data <- temp_data[, -grep("qc", colnames(temp_data))]
      if (!is.null(dim(temp_data))) {
        if (ncol(temp_data) == 0)
          temp_data <- NULL
      }
    }
    
    if (any(is.na(flow_data))) {
      if ((sum(is.na(flow_data)) / length(flow_data)) > na_thresh) {
        maf <- NA
        covaf <- NA
        maxan <- NA
        minan <- NA
        mspr <- NA
        msum <- NA
        mspwn <- NA
        mwin <- NA
        medwin <- NA
        medspr <- NA
        medsum <- NA
        covsp <- NA
        ransp <- NA
        ransm <- NA
        minwin <- NA
      } else {
        maf <- mean(flow_data, na.rm = TRUE)
        covaf <- sd(flow_data, na.rm = TRUE) / maf
        maxan <- max(flow_data, na.rm = TRUE)
        minan <- min(flow_data, na.rm = TRUE)
        mspr <- mean(flow_data[month(x$Event_Date) %in% c(9:11)], na.rm = TRUE)
        msum <- mean(flow_data[month(x$Event_Date) %in% c(12, 1:2)], na.rm = TRUE)
        mspwn <- mean(flow_data[month(x$Event_Date) %in% c(10:12)], na.rm = TRUE)
        mwin <- mean(flow_data[month(x$Event_Date) %in% c(5:7)], na.rm = TRUE)
        medwin <- median(flow_data[month(x$Event_Date) %in% c(6:8)], na.rm = TRUE)
        medspr <- median(flow_data[month(x$Event_Date) %in% c(9:11)], na.rm = TRUE)
        medsum <- median(flow_data[month(x$Event_Date) %in% c(12, 1:2)], na.rm = TRUE)
        covsp <- sd(flow_data[month(x$Event_Date) %in% c(10:12)], na.rm = TRUE) / mspwn
        ransp <- rolling_range(flow_data[month(x$Event_Date) %in% c(10:12)], 3)
        ransm <- rolling_range(flow_data[month(x$Event_Date) %in% c(12, 1:2)], 3)
        minwin <- min(flow_data[month(x$Event_Date) %in% c(6:8)], na.rm = TRUE)
      }
    } else {
      maf <- mean(flow_data, na.rm = TRUE)
      covaf <- sd(flow_data, na.rm = TRUE) / maf
      maxan <- max(flow_data, na.rm = TRUE)
      minan <- min(flow_data, na.rm = TRUE)
      mspr <- mean(flow_data[month(x$Event_Date) %in% c(9:11)], na.rm = TRUE)
      msum <- mean(flow_data[month(x$Event_Date) %in% c(12, 1:2)], na.rm = TRUE)
      mspwn <- mean(flow_data[month(x$Event_Date) %in% c(10:12)], na.rm = TRUE)
      mwin <- mean(flow_data[month(x$Event_Date) %in% c(5:7)], na.rm = TRUE)
      medwin <- median(flow_data[month(x$Event_Date) %in% c(6:8)], na.rm = TRUE)
      medspr <- median(flow_data[month(x$Event_Date) %in% c(9:11)], na.rm = TRUE)
      medsum <- median(flow_data[month(x$Event_Date) %in% c(12, 1:2)], na.rm = TRUE)
      covsp <- sd(flow_data[month(x$Event_Date) %in% c(10:12)], na.rm = TRUE) / mspwn
      ransp <- rolling_range(flow_data[month(x$Event_Date) %in% c(10:12)], 3)
      ransm <- rolling_range(flow_data[month(x$Event_Date) %in% c(12, 1:2)], 3)
      minwin <- min(flow_data[month(x$Event_Date) %in% c(6:8)], na.rm = TRUE)
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
    
    if (!is.null(temp_data)) {
      if (any(is.na(temp_data))) {
        if ((sum(is.na(temp_data)) / length(temp_data)) > na_thresh) {
          spwntmp <- NA
        } else {
          spwntmp <- mean(temp_data[month(x$Event_Date) %in% c(10:12)], na.rm = TRUE)
        }
      } else { 
        spwntmp <- mean(temp_data[month(x$Event_Date) %in% c(10:12)], na.rm = TRUE)
      }
    } else {
      spwntmp <- NA
    }
    
    out <- c(maf, mspr, msum, covaf,
             mspwn, mwin,
             covsp, maxan, minan,
             madp, cvdp, maxdp, mindp,
             ransp, ransm,
             medwin, medspr, medsum,
             minwin, spwntmp)

  }
  
  names(out) <- c("mannf_mld", "msprf_mld", "msumf_mld", "covaf_mld",
                  "mspwn_mld", "mwin_mld",
                  "covsp_mld", "maxan_mld", "minan_mld",
                  "madpth_m", "cvdpth_m", "maxdpth_m", "mindpth_m",
                  "rrang_spwn_mld", "rrang_sum_mld",
                  "median_win_mld", "median_spr_mld", "median_sum_mld",
                  "minwin_mld", "spwntmp_c")
  
  out
  
}

# pull out preceding 12 months of data
flow_tm1 <- function(x, flow, na_thresh = 0.2, year_lag = 0) {
  
  systmp <- paste0(tolower(x$SYSTEM), "_r", x$Reach)
  
  if (systmp %in% c("campaspe_r3", "campaspe_r4"))
    systmp <- "campaspe_r34"
  if (systmp %in% c("ovens_r1"))
    systmp <- "ovens_wangaratta"
  if (systmp %in% c("lowermurray_r1"))
    systmp <- "murray_yarrawonga"
  if (systmp %in% c("ovens_r1"))
    systmp <- "ovens_wangaratta"
  if (systmp %in% c("loddon_r5"))
    systmp <- "loddon_r4"
  if (systmp %in% c("broken_r4"))
    systmp <- "broken_r3"
  if (systmp %in% c("goulburn_r1"))
    systmp <- "goulburn_r4"
  
  if (any(names(flow) == systmp)) {
    flow_sub <- flow[names(flow) == systmp][[1]]
    tmp <- flow_sub[flow_sub$Event_Date < (x$Event_Date - years(year_lag)), ]
    tmp <- tail(tmp, 365)
    
    flow_yr <- tmp[, grep("discharge_ml", colnames(tmp))]
    if (is.matrix(flow_yr) | is.data.frame(flow_yr)) {
      if (ncol(flow_yr) > 1)
        flow_yr <- flow_yr[, -grep("qc", colnames(flow_yr))]
    }
    
  } else {
    tmp <- NULL
  }
  
  out <- calc_flow_metrics(tmp, na_thresh = na_thresh)
  
  if (any(names(flow) == systmp)) {
    x <- flow_sub
    flow_tmp <- x[, grep("discharge_ml", colnames(x))]
    if (is.matrix(flow_tmp) | is.data.frame(flow_tmp)) {
      if (ncol(flow_tmp) > 1)
        flow_tmp <- flow_tmp[, -grep("qc", colnames(flow_tmp))]
    }
    lt_medwin <- median(flow_tmp[month(x$Event_Date) %in% c(6:8)], na.rm = TRUE)
    lt_qwin <- quantile(flow_tmp[month(x$Event_Date) %in% c(5:7)], p = 0.1, na.rm = TRUE)
    st_medwin <- out["median_win_mld"]
    new_vars <- c(out["median_spr_mld"] / st_medwin,
                  out["median_sum_mld"] / st_medwin,
                  out["median_spr_mld"] / lt_medwin,
                  out["median_sum_mld"] / lt_medwin,
                  sum(flow_yr[month(tmp$Event_Date) %in% c(5:7)] < lt_qwin))
    names(new_vars) <- c("prop_spr_win", "prop_sum_win",
                         "prop_spr_lt_win", "prop_sum_lt_win",
                         "numlow_days") 
  } else {
    new_vars <- rep(NA, 5)
  }
  
  out <- c(out, new_vars)
  
  out

}

calc_flow_fn <- function(i, data, predictors, na_thresh = 0.2, year_lag = 0) {
  flow_tm1(data[i, ], predictors, na_thresh, year_lag = year_lag)
}

future::plan(future::multicore)
out_test <- future.apply::future_sapply(seq_len(nrow(alldat)),
                                        calc_flow_fn,
                                        data = alldat,
                                        predictors = predictors,
                                        na_thresh = 0.2)
flow_data <- t(out_test)

# repeat to get data lagged by 1 and 2 years
out_ym1 <- future.apply::future_sapply(seq_len(nrow(alldat)),
                                       calc_flow_fn,
                                       data = alldat,
                                       predictors = predictors,
                                       na_thresh = 0.2,
                                       year_lag = 1)
out_ym2 <- future.apply::future_sapply(seq_len(nrow(alldat)),
                                       calc_flow_fn,
                                       data = alldat,
                                       predictors = predictors,
                                       na_thresh = 0.2,
                                       year_lag = 2)
out_ym1 <- t(out_ym1)
out_ym2 <- t(out_ym2)
colnames(out_ym1) <- paste0(colnames(flow_data), "_ym1")
colnames(out_ym2) <- paste0(colnames(flow_data), "_ym2")
flow_data <- cbind(flow_data, out_ym1, out_ym2)

# easier to work with data.frame
flow_data <- as.data.frame(flow_data)

# clean up workspace
rm(i, out_test, date_parse, predictors, pred_list)
