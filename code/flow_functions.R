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

# calculate metrics based on flow data from survey year and two years prior
calculate_flow_metrics <- function(i, survey_data, flow_data, na_thresh = 0.2, year_lag = 0) {

  # which months correspond to which seasons?
  month_ids <- define_months()
  
  # which system and reach do we need data for?
  systmp <- paste0(tolower(survey_data$system[i]), "_r", survey_data$reach[i])

  # do these data exist? If not, replace with nearby reach  
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
  
  # pull out flow data for target year and two years prior
  if (any(names(flow_data) == systmp)) {
    flow_sub <- flow_data[names(flow_data) == systmp][[1]]
    sample_year <- year(survey_data$date_formatted[i])
    relevant_years <- sample_year - c(2:0)
    tmp <- flow_sub[year(flow_sub$Event_Date) %in% relevant_years, ]
  } else {
    tmp <- NULL
  }
  
  # initialise empty output
  out <- list(
    spawning_variability = NA,
    median_spring = NA,
    median_summer = NA,
    median_winter = NA,
    max_antecedent = NA,
    spawning_temp = NA,
    prop_spring_lt = NA,
    prop_summer_lt = NA,
    prop_max_antecedent_lt = NA,
    prop_winter_lt = NA,
    number_low_days = NA
  )
  
  # do the data exist for this year/system?
  if (!is.null(tmp)) {
    
    # pull out flow data if it exists, error otherwise
    flow_vector <- tmp[, grep("discharge_ml", colnames(tmp))]
    flow_all_years <- flow_sub[, grep("discharge_ml", colnames(flow_sub))]
    if (is.matrix(flow_vector) | is.data.frame(flow_vector)) {
      if (ncol(flow_vector) > 1)
        flow_vector <- flow_vector[, -grep("qc", colnames(flow_vector))]
      if (ncol(flow_all_years) > 1)
        flow_all_years <- flow_all_years[, -grep("qc", colnames(flow_all_years))]
    }
    if (is.null(flow_vector))
      stop("something wrong with discharge colnames", call. = FALSE)

    # pull out temperature data if it exists, set to NULL otherwise
    temp_vector <- tmp[, grep("temp", colnames(tmp))]
    if (is.matrix(temp_vector) | is.data.frame(temp_vector)) {
      if (ncol(temp_vector) > 1)
        temp_vector <- temp_vector[, -grep("qc", colnames(temp_vector))]
      if (!is.null(dim(temp_vector))) {
        if (ncol(temp_vector) == 0)
          temp_vector <- NULL
      }
    }
    
    # set up some subsets based on years
    contemporary <- year(tmp$Event_Date) == sample_year
    prior <- year(tmp$Event_Date) == (sample_year - 1L)
    antecedent <- year(tmp$Event_Date) == (sample_year - 2L)
    spanning <- (month(tmp$Event_Date) == 12 & year(tmp$Event_Date) == (sample_year - 1L)) |
      (month(tmp$Event_Date) %in% c(1:3) & year(tmp$Event_Date) == sample_year)
    
    # also need some based on seasons
    summer <- month(tmp$Event_Date) %in% month_ids$summer_months
    spawning <- month(tmp$Event_Date) %in% month_ids$spawning_months
    spring <- month(tmp$Event_Date) %in% month_ids$spring_months
    cool_season <- month(tmp$Event_Date) %in% month_ids$cooler_months
    
    # calculate flow metrics for different combinations of years and seasons
    if ((sum(is.na(flow_vector[prior])) / sum(prior)) <= na_thresh) {
      if (sum(prior & spawning) > 0)
        out$spawning_variability <- rolling_range(flow_vector[prior & spawning], 3)
      if (sum(prior & spring) > 0)
        out$median_spring <- median(flow_vector[prior & spring], na.rm = TRUE)
    }
    if (sum(spanning) > 0) {
      if ((sum(is.na(flow_vector[spanning])) / sum(spanning)) <= na_thresh)
        out$median_summer <- median(flow_vector[spanning & summer], na.rm = TRUE)
    }
    if (sum(contemporary) > 0) {
      if ((sum(is.na(flow_vector[contemporary])) / sum(contemporary)) <= na_thresh)
        out$median_winter <- median(flow_vector[contemporary & cool_season], na.rm = TRUE)
    }
    if ((sum(is.na(flow_vector[antecedent])) / sum(antecedent)) <= na_thresh)
      out$max_antecedent <- max(flow_vector[antecedent], na.rm = TRUE)
    
    # standardise some of the flow metrics by long-term averages
    full_year_all_years <- month(flow_sub$Event_Date) %in% month_ids$all_months
    cool_season_all_years <- month(flow_sub$Event_Date) %in% month_ids$cooler_months
    lt_med <- median(flow_all_years[full_year_all_years], na.rm = TRUE)
    lt_qcool <- quantile(flow_all_years[cool_season_all_years], p = 0.1, na.rm = TRUE)
    out$prop_spring_lt <- out$median_spring / lt_med
    out$prop_summer_lt <- out$median_summer / lt_med
    out$prop_max_antecedent_lt <- out$max_antecedent / lt_med
    out$prop_winter_lt <- out$median_winter / lt_med
    out$number_low_days <- sum(flow_vector[prior & cool_season] < lt_qcool)
    
    # temperature calculations filtered to <= na_thresh missing vlaues
    if (!is.null(temp_vector)) {
      temp_prior <- temp_vector[prior & spawning]
      if (length(temp_prior) > 0) {
        if ((sum(is.na(temp_prior)) / length(temp_prior)) <= na_thresh)
          out$spawning_temp <- mean(temp_prior, na.rm = TRUE)
      }
    }

  } 
  
  # return
  out
  
}
