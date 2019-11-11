# helpers for catch curve analysis

# helper function to match cols and fill empties
col_match <- function(x, names, classes = NULL) {
  
  idx <- match(names, colnames(x))
  if (any(is.na(idx)))
    idx[is.na(idx)] <- ncol(x) + 1
  x <- cbind(x, rep(NA, nrow(x)))
  
  out <- x[, idx]
  colnames(out) <- names
  
  out  
  
}

add_scientific_names <- function(x) {
  
  sci_names <- c(maccullochellapeelii = "Maccullochella peelii",
                 maccullochellamacquariensis = "Maccullochella macquariensis",
                 cyprinuscarpio = "Cyprinus carpio",
                 macquariaambigua = "Macquaria ambigua",
                 macquariaaustralasica = "Macquaria australasica",
                 melanotaeniafluviatilis = "Melanotaenia fluviatilis",
                 bidyanusbidyanus = "Bidyanus bidyanus")
  
  sci_names[x]
  
}

add_common_names <- function(x) {
  
  common_names <- c(maccullochellapeelii = "Murray cod",
                    maccullochellamacquariensis = "Trout cod",
                    cyprinuscarpio = "Common carp",
                    macquariaambigua = "Golden perch",
                    macquariaaustralasica = "Macquarie perch",
                    melanotaeniafluviatilis = "Murray river rainbowfish",
                    bidyanusbidyanus = "Silver perch") 
  
  common_names[x]
  
}

# filter a matrix down by column to remove correlated columns
remove_correlated <- function(flow) {
  
  # calculate correlations > 0.7 and remove one at a time, taking left-most of the columns
  #   correlated with the most other variables
  cor_vals <- cor(flow, use = "complete")
  sum_correlated <- apply(cor_vals, 2, function(x) sum(abs(x) > 0.7) - 1)
  while(any(sum_correlated > 0)) {
    cor_vals <- cor_vals[-max(which(sum_correlated == max(sum_correlated))), -max(which(sum_correlated == max(sum_correlated)))]
    sum_correlated <- apply(cor_vals, 2, function(x) sum(abs(x) > 0.7) - 1)
  }
  flow[, match(names(sum_correlated), colnames(flow))]
  
}

# define some helper functions
max_fun <- function(x) {
  
  out <- NA
  
  if (any(!is.na(x)))
    out <- max(x, na.rm = TRUE)
  
  out
  
}

na_replace_fun <- function(x) {
  
  if (any(is.na(x))) 
    x[is.na(x)] <- mean(x, na.rm = TRUE)

  x
  
}

calculate_length_conversions <- function(data) {
  
  # pull out species names
  sp_names <- unique(data$Common.Name)
  
  # set up output list
  length_weight_conversion <- list()
  
  # loop through each species
  for (i in seq_along(sp_names)) {
    
    # subset for species i
    dat <- data[which(data$Common.Name == sp_names[i]), ]
    
    # log transform length and weight
    rows_to_rm <- NULL
    if (any(is.na(dat$totallength)))
      rows_to_rm <- c(rows_to_rm, which(is.na(dat$totallength)))
    if (any(is.na(dat$WEIGHT)))
      rows_to_rm <- c(rows_to_rm, which(is.na(dat$WEIGHT)))
    if (any(dat$totallength <= 0, na.rm = TRUE))
      rows_to_rm <- c(rows_to_rm, which(dat$totallength <= 0))
    if (any(dat$WEIGHT <= 0, na.rm = TRUE))
      rows_to_rm <- c(rows_to_rm, which(dat$WEIGHT <= 0))
    if (length(rows_to_rm))
      dat <- dat[-rows_to_rm, ]
    log_length <- log(dat$totallength)
    log_weight <- log(dat$WEIGHT)
    
    # fit linear model to log-log transformed data
    if (length(log_weight) > 10)
      mod_tmp <- lm(log_weight ~ log_length)
    
    # pull out coefs and other stats
    if (length(log_weight) > 10) {
      coefs_tmp <- coef(mod_tmp)
      n_obs <- length(log_length)
      r2_vals <- summary(mod_tmp)$r.squared
      resid_tmp <- residuals(mod_tmp)
    } else {
      coefs_tmp <- NA
      n_obs <- 0
      r2_vals <- NA
      resid_tmp <- NA
    }
    
    # save outputs
    length_weight_conversion[[sp_names[i]]] <-
      list(n = n_obs,
           coef = coefs_tmp,
           length = dat$totallength,
           weight = dat$WEIGHT,
           r2 = r2_vals,
           resid = resid_tmp)
    
  }
  
  coefs_all <- lapply(length_weight_conversion, function(x) x$coef)
  coefs_all <- coefs_all[-which(is.na(coefs_all))]
  
  length_weight_conversion$generic <- apply(matrix(unlist(coefs_all), ncol = 2, byrow = TRUE),
                                            2, mean)
  
  length_weight_conversion
  
}

impute_weights <- function(data, length_conversions) {
  
  # fill missing weights
  sp_tmp <- unique(data$Common.Name)
  for (i in seq_along(sp_tmp)) {
    
    # subset to single species
    dat_tmp <- data[which(data$Common.Name == sp_tmp[i]), ]
    
    # check to make sure some weights are missing      
    if (any(is.na(dat_tmp$WEIGHT))) {
      
      # use species-specific equation if it exists, generic otherwise
      if (length_conversions[[sp_tmp[i]]]$n) {
        coefs <- length_conversions[[sp_tmp[i]]]$coef
      } else {
        coefs <- length_conversions$generic
      }
      
      # subset NA observations
      na_sub <- which(is.na(dat_tmp$WEIGHT))
      
      # estimate weight from length
      dat_tmp$WEIGHT[na_sub] <- exp(coefs[1] + coefs[2] * log(dat_tmp$totallength[na_sub]))
      
      # return estimated data to full data set
      data[which(data$Common.Name == sp_tmp[i]), ] <- dat_tmp
      
    }
  }
  
  data
  
}

clean_reaches <- function(data) {
  
  # arrange alldat to have more consistent set of variables
  data <- data.frame(date = data$Date,
                     year = data$YEAR,
                     site = data$SITE_CODE,
                     system = data$SYSTEM,
                     reach = data$Reach,
                     species = data$Common.Name,
                     length = data$totallength,
                     weight = data$WEIGHT,
                     abundance = data$Total.Sampled,
                     intensity = (data$total_no_passes * data$seconds))
  
  data$reach_alt <- data$reach
  data$reach_alt <- ifelse(data$system == "BROKEN",
                           ifelse(data$reach_alt == 4, 5, data$reach_alt), # downstream 
                           data$reach_alt)
  data$reach_alt <- ifelse(data$system == "THOMSON",
                           ifelse(data$reach_alt == 2, 3,             # downstream 
                                  ifelse(data$reach_alt == 6, 5,      # upstream
                                         data$reach_alt)),
                           data$reach_alt)
  data$reach_alt <- ifelse(data$system == "LODDON",
                           ifelse(data$reach_alt == 2, 3,             # downstream
                                  ifelse(data$reach_alt == 5, 4,      # upstream 
                                         data$reach_alt)),
                           data$reach_alt)
  
  data
  
}

# internal functon for rolling_range (below)
get_range <- function(x) {
  
  out <- NA
  
  if (any(!is.na(x)))
    out <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    
  out
  
}

# alternative internal functon for rolling_range (below)
get_range2 <- function(x) {
  
  out <- NA
  
  if (any(!is.na(x))) {
    out <- 0
    if (min(x, na.rm = TRUE) > 0)
      out <- max(x, na.rm = TRUE) / min(x, na.rm = TRUE)
  }
  
  out
  
}

# calculate maximum rolling range in a variable over a set lag or lead time
rolling_range <- function(data, lag, variable = NULL){ 
  
  if (length(dim(data)) < 2) {
    data <- as.matrix(data, ncol = 1)
    variable <- 1
  }
  
  if (ncol(data) == 1)
    variable <- 1
  
  if (is.null(variable))
    stop("`variable` must be provided if `data` has more than one column", call. = FALSE)
  
  nrows <- nrow(data)
  
  idx <- sapply(rev(seq_len(lag)) - 1, function(x) rep(x, nrows))
  idx <- sweep(idx, 1, seq_len(nrows), "+")
  idx <- ifelse(idx > nrows, NA, idx)
  
  df <- matrix(data[c(idx), variable], nrow = nrows)
  
  diff <- apply(df, 1, get_range2)
  
  max(diff, na.rm = TRUE)
  
}   

# function to rename river systems from sites
system_switch_fun <- function(x) {
  
  x <- as.character(x)
  
  out <- rep(NA, length(x))
  for (i in seq_along(x)) {
    
    out[i] <- switch(substr(x[i], 1, 2),
                     "LO" = "LODDON",
                     "GO" = "GOULBURN",
                     "BR" = "BROKEN",
                     "PC" = "PYRAMIDCK",
                     "CA" = "CAMPASPE")
    
  }
  
  out
  
}

# set month definitions
define_months <- function() {
  list(
    all_months = c(1:12),
    summer_months = c(12, 1:3),
    cooler_months = c(5:7),
    winter_months = c(6:8),
    spawning_months = c(10:12),
    spring_months = c(9:11)
  )
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
