
# calculate flow metrics
calc_flow_metrics_ccr <- function(x, sample_year, na_thresh = 0.2) {
  
  month_ids <- define_months()
  
  # initialise empty output
  out <- list(spawning_variability = NA,
              median_spring = NA,
              median_summer = NA,
              median_winter = NA,
              max_antecedent = NA,
              spawning_temp = NA)
  
  if (!is.null(x)) {
    
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
    
    # discharge calculations filtered to <= na_thresh missing values
    if (!is.null(flow_data)) {
      flow_contemporary <- flow_data[year(x$Event_Date) == sample_year]
      flow_prior <- flow_data[year(x$Event_Date) == (sample_year - 1L)]
      flow_antecedent <- flow_data[year(x$Event_Date) == (sample_year - 2L)]
      flow_spanning <- c(flow_data[month(x$Event_Date) == 12 & year(x$Event_Date) == (sample_year - 1L)],
                         flow_data[month(x$Event_Date) < 5 & year(x$Event_Date) == sample_year])
      if ((sum(is.na(flow_prior)) / length(flow_prior)) <= na_thresh) {
        out$spawning_variability <- rolling_range(flow_prior[month(x$Event_Date) %in% month_ids$spawning_months], 3)
        out$median_spring <- median(flow_prior[month(x$Event_Date) %in% month_ids$spring_months], na.rm = TRUE)
      }
      if (length(flow_spanning) > 0) {
        if ((sum(is.na(flow_spanning)) / length(flow_spanning)) <= na_thresh)
          out$median_summer <- median(flow_spanning[month(x$Event_Date) %in% month_ids$summer_months], na.rm = TRUE)
      }
      if (length(flow_contemporary) > 0) {
        if ((sum(is.na(flow_contemporary)) / length(flow_contemporary)) <= na_thresh)
          out$median_winter <- median(flow_contemporary[month(x$Event_Date) %in% month_ids$cooler_months], na.rm = TRUE)
      }
      if ((sum(is.na(flow_antecedent)) / length(flow_antecedent)) <= na_thresh)
        out$max_antecedent <- max(flow_antecedent, na.rm = TRUE)
    }
    
    if (!is.null(temp_data)) {
      temp_prior <- temp_data[year(x$Event_Date) == (sample_year - 1L) &
                                month(x$Event_Date) %in% month_ids$spawning_months]
      if (length(temp_prior) > 0) {
        if ((sum(is.na(temp_prior)) / length(temp_prior)) <= na_thresh)
          out$spawning_temp <- mean(temp_prior, na.rm = TRUE)
      }
    }
    
  }
  
  out
  
}

# pull out preceding 12 months of data
flow_ccr <- function(x, flow, na_thresh = 0.2, year_lag = 0) {
  
  month_ids <- define_months()
  
  systmp <- paste0(tolower(x$system), "_r", x$reach)
  
  if (systmp %in% c("campaspe_r3", "campaspe_r4"))
    systmp <- "campaspe_r34"
  if (systmp %in% c("ovens_r1"))
    systmp <- "ovens_wangaratta"
  if (systmp %in% c("murray_r1"))
    systmp <- "murray_yarrawonga"
  if (systmp %in% c("loddon_r5"))
    systmp <- "loddon_r4"
  if (systmp %in% c("broken_r4"))
    systmp <- "broken_r3"
  if (systmp %in% c("goulburn_r1"))
    systmp <- "goulburn_r5"
  if (systmp %in% c("goulburn_r4"))
    systmp <- "goulburn_r5"
  
  if (any(names(flow) == systmp)) {
    flow_sub <- flow[names(flow) == systmp][[1]]
    sample_year <- year(x$date_formatted)
    relevant_years <- sample_year - c(2:0)
    tmp <- flow_sub[year(flow_sub$Event_Date) %in% relevant_years, ]
  } else {
    tmp <- NULL
  }
  
  out <- calc_flow_metrics_ccr(tmp, sample_year = sample_year, na_thresh = na_thresh)
  
  if (any(names(flow) == systmp)) {
    x <- flow_sub
    flow_tmp <- x[, grep("discharge_ml", colnames(x))]
    if (is.matrix(flow_tmp) | is.data.frame(flow_tmp)) {
      if (ncol(flow_tmp) > 1)
        flow_tmp <- flow_tmp[, -grep("qc", colnames(flow_tmp))]
    }
    lt_med <- median(flow_tmp[month(x$Event_Date) %in% month_ids$all_months], na.rm = TRUE)
    out$prop_spring_lt <- out$median_spring / lt_med
    out$prop_summer_lt <- out$median_summer / lt_med
    out$prop_max_antecedent_lt <- out$max_antecedent / lt_med
    out$prop_winter_lt <- out$median_winter / lt_med
  } else {
    out$prop_spring_lt <- NA
    out$prop_summer_lt <- NA
    out$prop_max_antecedent_lt <- NA
    out$prop_winter_lt <- NA
  }
  
  out
  
}

calc_flow <- function(i, data, predictors, na_thresh = 0.2, year_lag = 0) {
  flow_ccr(data[i, ], predictors, na_thresh, year_lag = year_lag)
}

# calculate flow metrics
calc_flow_metrics <- function(x, na_thresh = 0.2) {
  
  month_ids <- define_months()
  
  if (is.null(x)) {
    out <- rep(NA, 22)
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
        mcool <- NA
        mwin <- NA
        medcool <- NA
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
        mspr <- mean(flow_data[month(x$Event_Date) %in% month_ids$spring_months], na.rm = TRUE)
        msum <- mean(flow_data[month(x$Event_Date) %in% month_ids$summer_months], na.rm = TRUE)
        mspwn <- mean(flow_data[month(x$Event_Date) %in% month_ids$spawning_months], na.rm = TRUE)
        mcool <- mean(flow_data[month(x$Event_Date) %in% month_ids$cooler_months], na.rm = TRUE)
        mwin <- mean(flow_data[month(x$Event_Date) %in% month_ids$winter_months], na.rm = TRUE)
        medcool <- median(flow_data[month(x$Event_Date) %in% month_ids$cooler_months], na.rm = TRUE)
        medwin <- median(flow_data[month(x$Event_Date) %in% month_ids$winter_months], na.rm = TRUE)
        medspr <- median(flow_data[month(x$Event_Date) %in% month_ids$spring_months], na.rm = TRUE)
        medsum <- median(flow_data[month(x$Event_Date) %in% month_ids$summer_months], na.rm = TRUE)
        covsp <- sd(flow_data[month(x$Event_Date) %in% month_ids$spawning_months], na.rm = TRUE) / mspwn
        ransp <- rolling_range(flow_data[month(x$Event_Date) %in% month_ids$spawning_months], 3)
        ransm <- rolling_range(flow_data[month(x$Event_Date) %in% month_ids$summer_months], 3)
        minwin <- min(flow_data[month(x$Event_Date) %in% month_ids$winter_months], na.rm = TRUE)
      }
    } else {
      maf <- mean(flow_data, na.rm = TRUE)
      covaf <- sd(flow_data, na.rm = TRUE) / maf
      maxan <- max(flow_data, na.rm = TRUE)
      minan <- min(flow_data, na.rm = TRUE)
      mspr <- mean(flow_data[month(x$Event_Date) %in% month_ids$spring_months], na.rm = TRUE)
      msum <- mean(flow_data[month(x$Event_Date) %in% month_ids$summer_months], na.rm = TRUE)
      mspwn <- mean(flow_data[month(x$Event_Date) %in% month_ids$spawning_months], na.rm = TRUE)
      mcool <- mean(flow_data[month(x$Event_Date) %in% month_ids$cooler_months], na.rm = TRUE)
      mwin <- mean(flow_data[month(x$Event_Date) %in% month_ids$winter_months], na.rm = TRUE)
      medcool <- median(flow_data[month(x$Event_Date) %in% month_ids$cooler_months], na.rm = TRUE)
      medwin <- median(flow_data[month(x$Event_Date) %in% month_ids$winter_months], na.rm = TRUE)
      medspr <- median(flow_data[month(x$Event_Date) %in% month_ids$spring_months], na.rm = TRUE)
      medsum <- median(flow_data[month(x$Event_Date) %in% month_ids$summer_months], na.rm = TRUE)
      covsp <- sd(flow_data[month(x$Event_Date) %in% month_ids$spawning_months], na.rm = TRUE) / mspwn
      ransp <- rolling_range(flow_data[month(x$Event_Date) %in% month_ids$spawning_months], 3)
      ransm <- rolling_range(flow_data[month(x$Event_Date) %in% month_ids$summer_months], 3)
      minwin <- min(flow_data[month(x$Event_Date) %in% month_ids$winter_months], na.rm = TRUE)
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
          spwntmp <- mean(temp_data[month(x$Event_Date) %in% month_ids$spawning_months], na.rm = TRUE)
        }
      } else { 
        spwntmp <- mean(temp_data[month(x$Event_Date) %in% month_ids$spawning_months], na.rm = TRUE)
      }
    } else {
      spwntmp <- NA
    }
    
    out <- c(maf, mspr, msum, covaf,
             mspwn, mcool, mwin,
             covsp, maxan, minan,
             madp, cvdp, maxdp, mindp,
             ransp, ransm,
             medcool, medwin, medspr, medsum,
             minwin, spwntmp)
    
  }
  
  names(out) <- c("mannf_mld", "msprf_mld", "msumf_mld", "covaf_mld",
                  "mspwn_mld", "mcool_mld", "mwin_mld",
                  "covsp_mld", "maxan_mld", "minan_mld",
                  "madpth_m", "cvdpth_m", "maxdpth_m", "mindpth_m",
                  "rrang_spwn_mld", "rrang_sum_mld",
                  "median_cool_mld", "median_win_mld",
                  "median_spr_mld", "median_sum_mld",
                  "minwin_mld", "spwntmp_c")
  
  out
  
}

# pull out preceding 12 months of data
flow_tm1 <- function(x, flow, na_thresh = 0.2, year_lag = 0) {
  
  month_ids <- define_months()
  
  systmp <- paste0(tolower(x$system), "_r", x$reach)
  
  if (systmp %in% c("campaspe_r3", "campaspe_r4"))
    systmp <- "campaspe_r34"
  if (systmp %in% c("ovens_r1"))
    systmp <- "ovens_wangaratta"
  if (systmp %in% c("murray_r1"))
    systmp <- "murray_yarrawonga"
  if (systmp %in% c("loddon_r5"))
    systmp <- "loddon_r4"
  if (systmp %in% c("broken_r4"))
    systmp <- "broken_r3"
  if (systmp %in% c("goulburn_r1"))
    systmp <- "goulburn_r4"
  
  if (any(names(flow) == systmp)) {
    flow_sub <- flow[names(flow) == systmp][[1]]
    tmp <- flow_sub[flow_sub$Event_Date < (x$date_formatted - years(year_lag)), ]
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
    lt_med <- median(flow_tmp[month(x$Event_Date) %in% month_ids$all_months], na.rm = TRUE)
    lt_qcool <- quantile(flow_tmp[month(x$Event_Date) %in% month_ids$cooler_months], p = 0.1, na.rm = TRUE)
    st_medwin <- out["median_win_mld"]
    new_vars <- c(out["median_spr_mld"] / st_medwin,
                  out["median_sum_mld"] / st_medwin,
                  out["median_spr_mld"] / lt_med,
                  out["median_sum_mld"] / lt_med,
                  out["maxan_mld"] / lt_med,
                  out["median_cool_mld"] / lt_med,
                  out["median_win_mld"] / lt_med,
                  sum(flow_yr[month(tmp$Event_Date) %in% month_ids$cooler_months] < lt_qcool))
    names(new_vars) <- c("prop_spr_win", "prop_sum_win",
                         "prop_spr_lt", "prop_sum_lt",
                         "prop_maxan_lt",
                         "prop_cool_lt", "prop_win_lt",
                         "numlow_days") 
  } else {
    new_vars <- rep(NA, 8)
  }
  
  out <- c(out, new_vars)
  
  out
  
}

calc_flow_fn <- function(i, data, predictors, na_thresh = 0.2, year_lag = 0) {
  flow_tm1(data[i, ], predictors, na_thresh, year_lag = year_lag)
}
