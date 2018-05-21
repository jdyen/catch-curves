# helpers for catch curve analysis
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

calc_flow_pc <- function(flow, scale = FALSE, nkeep = 3) {
  
  flow <- apply(flow, 2, na_replace_fun)
  if (scale)
    flow <- scale(flow)
  pc_out <- princomp(flow)
  
  pc_out
  
}

# do this for all ages, then just focus on 1-3 year olds for some analyses
catch_curve_fun <- function(x, sp) {
  
  x_sub <- x[which(x$species == sp), ]
  
  flow_tmp <- x_sub[, grep("mannf$", colnames(x_sub)):grep("cspwn_tm2", colnames(x_sub))]
  
  bins <- c(0, seq_len(max(x_sub$age, na.rm = TRUE) + 1)) + 0.5
  obs_unique <- paste(x_sub$system, paste0("reach", x_sub$reach), x_sub$year, sep = "_")
  
  obs_vec <- unique(obs_unique)
  
  out <- matrix(NA, nrow = length(obs_vec), ncol = (length(bins) - 1))
  flow_data <- matrix(NA, nrow = length(obs_vec), ncol = ncol(flow_tmp))
  for (i in seq_along(obs_vec)) {
    dat_sub <- x_sub[which(obs_unique == obs_vec[i]), ]
    out[i, ] <- hist(dat_sub$age, breaks = bins, plot = FALSE)$counts
    flow_data[i, ] <- apply(flow_tmp[which(obs_unique == obs_vec[i]), ], 2, mean, na.rm = TRUE)
  }
  colnames(flow_data) <- colnames(flow_tmp)
  colnames(out) <- seq_len(ncol(out))
  rownames(out) <- obs_vec
  site_split <- strsplit(obs_vec, "_")
  site_info <- data.frame(system = sapply(site_split, function(x) x[1]),
                          reach = sapply(site_split, function(x) as.numeric(substr(x[2],
                                                                                   start = nchar(x[2]),
                                                                                   stop = nchar(x[2])))),
                          year = sapply(site_split, function(x) as.numeric(x[3])))
  
  out <- list(age_dist = out,
              info = site_info,
              flow = flow_data,
              bins = bins)
  
}

# do this for all ages, then just focus on 1-3 year olds for some analyses
size_dist_fun <- function(x, sp, nbins = 5) {
  
  x_sub <- x[which(x$species == sp), ]
  
  bins <- c(0, 100, 250, 400, 600, 1000, 1500, 15000, max(x_sub$weight, na.rm = TRUE))
  obs_unique <- paste(x_sub$system, paste0("reach", x_sub$reach), x_sub$year, sep = "_")
  
  obs_vec <- unique(obs_unique)
  
  out <- matrix(NA, nrow = length(obs_vec), ncol = (length(bins) - 1))
  for (i in seq_along(obs_vec)) {
    dat_sub <- x_sub[which(obs_unique == obs_vec[i]), ]
    out[i, ] <- hist(dat_sub$weight, breaks = bins, plot = FALSE)$counts
  }
  colnames(out) <- seq_len(ncol(out))
  rownames(out) <- obs_vec
  site_split <- strsplit(obs_vec, "_")
  site_info <- data.frame(system = sapply(site_split, function(x) x[1]),
                          reach = sapply(site_split, function(x) as.numeric(substr(x[2],
                                                                                   start = nchar(x[2]),
                                                                                   stop = nchar(x[2])))),
                          year = sapply(site_split, function(x) as.numeric(x[3])))
  
  out <- list(size_dist = out,
              info = site_info,
              bins = bins)
  
}

list_full_spp_names <- function() {
  
  # set up species names for plots
  sp_names <- data.frame(full = c("Australian bass",
                                  "Australian grayling",
                                  "Australian smelt",
                                  "Black bream",
                                  "Bony bream",
                                  "Brown trout",
                                  "Carp",
                                  "Carp gudgeon",
                                  "Common galaxias",
                                  "Dwarf flathead gudgeon",
                                  "Eastern gambusia",
                                  "Eel",
                                  "Estuary perch",
                                  "Flathead galaxias",
                                  "Flathead gudgeon",
                                  "Golden perch",
                                  "Long-finned eel",
                                  "Luderick",
                                  "Mountain galaxias",
                                  "Murray cod",
                                  "Murray river rainbowfish",
                                  "Obscure galaxias",
                                  "Oriental weather loach",
                                  "Pouched lamprey",
                                  "Pygmy perch",
                                  "Rainbow trout",
                                  "Redfin",
                                  "River blackfish",
                                  "River garfish",
                                  "Roach",
                                  "Sea mullet",
                                  "Short finned eel",
                                  "Short headed lamprey",
                                  "Silver perch",
                                  "Southern pygmy perch",
                                  "Tench",
                                  "Trout cod",
                                  "Tupong",
                                  "Two spined blackfish",
                                  "Unspecked hardyhead",
                                  "Variegated pygmy perch",
                                  "Western carp gudgeon",
                                  "Yarra pygmy perch",
                                  "Yellow eyed mullet"),
                         code = as.character(levels(alldat$species)))
  
  sp_names
  
}

filter_fun <- function(x) {
  
  x_id <- unique(x$ID)
  x_len <- rep(NA, length(x_id))
  x_age <- rep(NA, length(x_id))
  
  for (i in seq_along(x_id)) {
    
    x_sub <- x[which(x$ID == x_id[i]), ]
    x_len[i] <- max(x_sub$Length)
    x_age[i] <- max(x_sub$Age)
    
  }
  
  x_filtered <- data.frame(id = x_id,
                           length = x_len,
                           age = x_age)
  
  x_filtered
  
}

coef_extract <- function(x_filtered) {
  
  params <- coef(lm(x_filtered$age ~ x_filtered$length))
  names(params) <- c("intercept", "slope")
  
  params
  
}

length_to_age <- function(x, params) {
  
  out <- round(params[1] + params[2] * x)
  out <- ifelse(out <= 0, NA, out)
  
  out
  
}

# functions to read in and load data

get_data <- function(file) {
  
  read.csv(file)
  
}

compile_data <- function(vefmap_data, snags_data, ovens_data) {
  
  # clean snags data set
  snags_data$date_new <- format(dmy_hms(snags_data$surveydate), format = "%d/%m/%Y")
  snags_data$YEAR <- sapply(strsplit(snags_data$date_new, "/"),
                            function(x) x[3])
  snags_data$taxonname <- as.character(snags_data$taxonname)
  snags_data$taxonname <- ifelse(snags_data$taxonname == "Yellowbelly",
                                 "Golden perch",
                                 snags_data$taxonname)
  snags_data$taxonname <- factor(snags_data$taxonname)
  snags_data <- data.frame(SYSTEM = rep("LOWERMURRAY", nrow(snags_data)),
                           SITE_CODE = paste0("Lm", snags_data$idsite),
                           Reach = rep(1, nrow(snags_data)),
                           geartype = factor(rep("EF/Boat"), nrow(snags_data)),
                           Event_Date = snags_data$date_new,
                           Pass.No = rep(1, nrow(snags_data)),
                           total_no_passes = rep(1, nrow(snags_data)),
                           seconds = snags_data$seconds,
                           Common.Name = snags_data$taxonname,
                           Scientific.Name = snags_data$Scientific.Name,
                           totallength = snags_data$totallength,
                           WEIGHT = snags_data$weight,
                           Total.Sampled = rep(1, nrow(snags_data)),
                           VEFMAP.Stage = rep(NA, nrow(snags_data)),
                           YEAR = as.integer(snags_data$YEAR))
  
  # clean ovens data
  ovens_data$date_new <- format(dmy(ovens_data$date), format = "%d/%m/%Y")
  ovens_data$YEAR <- sapply(strsplit(ovens_data$date_new, "/"),
                            function(x) x[3])
  ovens_data$species <- as.character(ovens_data$species)
  ovens_data$species <- ifelse(ovens_data$species == "Maccullochella peelii ",
                               "Maccullochella peelii peelii",
                               ovens_data$species)
  ovens_data$species <- ifelse(ovens_data$species == "Maccullochella peelii",
                               "Maccullochella peelii peelii",
                               ovens_data$species)
  ovens_data$species <- factor(ovens_data$species)
  ovens_data$common_name <- vefmap_data$Common.Name[match(ovens_data$species, vefmap_data$Scientific.Name)]
  ovens_data <- data.frame(SYSTEM = rep("OVENS", nrow(ovens_data)),
                           SITE_CODE = paste0("Ov", ovens_data$site),
                           Reach = rep(1, nrow(ovens_data)),
                           geartype = ovens_data$gear_type,
                           Event_Date = ovens_data$date_new,
                           Pass.No = rep(1, nrow(ovens_data)),
                           total_no_passes = rep(1, nrow(ovens_data)),
                           seconds = ovens_data$electro_seconds,
                           Common.Name = ovens_data$common_name,
                           Scientific.Name = ovens_data$species,
                           totallength = ovens_data$total_length_mm,
                           WEIGHT = ovens_data$weight_g,
                           Total.Sampled = ovens_data$no_collected,
                           VEFMAP.Stage = rep(NA, nrow(ovens_data)),
                           YEAR = as.integer(ovens_data$YEAR))
  
  # compile all data
  data <- rbind(vefmap_data, ovens_data, snags_data)
  
  # reformat dates
  data$Date <- dmy(data$Event_Date)
  
  data
  
}

subset_systems <- function(data, systems = NULL) {
  
  if (is.null(systems)) {  
    # set systems of interest
    systems <- c("BROKEN",
                 "LOWERMURRAY",
                 "CAMPASPE",
                 "GLENELG",
                 "GOULBURN",
                 "LODDON",
                 "THOMSON",
                 "OVENS")
  } 
  
  data[-which(is.na(match(data$SYSTEM, systems))), ]
  
} 

clean_spp <- function(data) {
  
  # clean up common names
  data$Common.Name <- tolower(data$Common.Name)
  data$Common.Name <- gsub(" ", "", data$Common.Name)
  data$Common.Name <- gsub("-[[:digit:]]*", "", data$Common.Name)
  data$Common.Name <- gsub("/", "", data$Common.Name)
  data$Common.Name <- gsub("sp\\.", "", data$Common.Name)
  data$Common.Name <- gsub("flatheadedgudgeon", "flatheadgudgeon", data$Common.Name)
  data$Common.Name <- gsub("europeancarp", "carp", data$Common.Name)
  data$Common.Name <- gsub("redfinperch", "redfin", data$Common.Name)
  data$Common.Name <- ifelse(data$Common.Name == "weatherloach",
                             "orientalweatherloach",
                             data$Common.Name)
  data$Common.Name <- ifelse(data$Common.Name == "rainbowfish",
                             "murrayriverrainbowfish",
                             data$Common.Name)
  data$Common.Name <- ifelse(data$Common.Name == "hardyhead",
                             "unspeckedhardyhead",
                             data$Common.Name)
  data$Common.Name <- ifelse(data$Common.Name == "blackfish",
                             "riverblackfish", 
                             data$Common.Name)
  
  # remove spp not of interest
  sp_to_rm <- c("carpgoldfishhybrid",
                "carpgudgeoncomplex",
                "commonyabby",
                "easternsnakeneckedturtle",
                "freshwatershrimp",
                "glenelgspinyfreshwatercrayfish",
                "goldfish",
                "hypseleostris",
                "murrayspinycrayfish",
                "nofish",
                "unidentifiedcod",
                "unknown",
                "",
                "codspp.",
                "yabbie",
                "longneckturtle",
                "gudgeons",
                "murraycray",
                "murraycodtroutcodhybrid")
  data <- data[-which(!is.na(match(data$Common.Name, sp_to_rm))), ]
  if (any(is.na(data$Common.Name)))
    data <- data[-which(is.na(data$Common.Name)), ]
  
  data
  
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
    length_weight_conversion[[sp_names[i]]] <- list(n = n_obs,
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

calculate_flow_stats <- function(files) {
  
  # load predictor (flow) files into a list from all files with "to_use" suffix
  pred_list <- files
  predictors <- vector("list", length = length(pred_list))
  for (i in seq_along(pred_list)) 
    predictors[[i]] <- read.csv(paste0("./data/flow-data/", pred_list[i]))
  names(predictors) <- sapply(strsplit(pred_list, "_"), function(x) paste(x[1], x[2], sep = "_"))
  
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
  
  mean_flow_stats
  
}

add_flow_data <- function(data, mean_flow_stats) {
  
  mean_ann_flow_vec <- mean_spr_flow_vec <- mean_sum_flow_vec <- rep(NA, nrow(data))
  maxm_ann_flow_vec <- covn_ann_flow_vec <- rep(NA, nrow(data))
  mean_spwn_flow_vec <- cov_spwn_flow_vec <- rep(NA, nrow(data))
  mean_spwn_flow_tm1 <- covn_ann_flow_tm1 <- cov_spwn_flow_tm1 <- rep(NA, nrow(data))
  mean_spwn_flow_tm2 <- covn_ann_flow_tm2 <- cov_spwn_flow_tm2 <- rep(NA, nrow(data))
  mean_ann_flow_tm1 <- mean_spr_flow_tm1 <- mean_sum_flow_tm1 <- maxm_ann_flow_tm1 <- rep(NA, nrow(data))
  mean_ann_flow_tm2 <- mean_spr_flow_tm2 <- mean_sum_flow_tm2 <- maxm_ann_flow_tm2 <- rep(NA, nrow(data))
  for (i in seq_len(nrow(data))) {
    row_tmp <- which((mean_flow_stats$system == data$system[i]) &
                       (mean_flow_stats$reach == data$reach[i]) &
                       (mean_flow_stats$year == data$year[i]))
    row_tmp2 <- which((mean_flow_stats$system == data$system[i]) &
                        (mean_flow_stats$reach == data$reach[i]) &
                        (mean_flow_stats$year == (data$year[i] - 1)))
    row_tmp3 <- which((mean_flow_stats$system == data$system[i]) &
                        (mean_flow_stats$reach == data$reach[i]) & 
                        (mean_flow_stats$year == (data$year[i] - 2)))
    if (!length(row_tmp)) {
      row_tmp <- which((mean_flow_stats$system == data$system[i]) &
                         (mean_flow_stats$reach == data$reach_alt[i]) &
                         (mean_flow_stats$year == data$year[i]))
    }
    if (!length(row_tmp2)) {
      row_tmp2 <- which((mean_flow_stats$system == data$system[i]) &
                          (mean_flow_stats$reach == data$reach_alt[i]) &
                          (mean_flow_stats$year == (data$year[i] - 1)))
    }
    if (!length(row_tmp3)) {
      row_tmp3 <- which((mean_flow_stats$system == data$system[i]) &
                          (mean_flow_stats$reach == data$reach_alt[i]) & 
                          (mean_flow_stats$year == (data$year[i] - 2)))
    }
    
    if (length(row_tmp)) {
      mean_ann_flow_vec[i] <- mean_flow_stats$mean_ann_flow[row_tmp]
      mean_spr_flow_vec[i] <- mean_flow_stats$mean_spr_flow[row_tmp]
      mean_sum_flow_vec[i] <- mean_flow_stats$mean_sum_flow[row_tmp]
      maxm_ann_flow_vec[i] <- mean_flow_stats$max_ann_flow[row_tmp]
      covn_ann_flow_vec[i] <- mean_flow_stats$cov_ann_flow[row_tmp]
      mean_spwn_flow_vec[i] <- mean_flow_stats$mean_spwn_flow[row_tmp]
      cov_spwn_flow_vec[i] <- mean_flow_stats$cov_spwn_flow[row_tmp]
    }
    if (length(row_tmp2)) {
      mean_ann_flow_tm1[i] <- mean_flow_stats$mean_ann_flow[row_tmp2]
      mean_spr_flow_tm1[i] <- mean_flow_stats$mean_spr_flow[row_tmp2]
      mean_sum_flow_tm1[i] <- mean_flow_stats$mean_sum_flow[row_tmp2]
      maxm_ann_flow_tm1[i] <- mean_flow_stats$max_ann_flow[row_tmp2]
      covn_ann_flow_tm1[i] <- mean_flow_stats$cov_ann_flow[row_tmp2]
      mean_spwn_flow_tm1[i] <- mean_flow_stats$mean_spwn_flow[row_tmp2]
      cov_spwn_flow_tm1[i] <- mean_flow_stats$cov_spwn_flow[row_tmp2]
    }
    if (length(row_tmp3)) {
      mean_ann_flow_tm2[i] <- mean_flow_stats$mean_ann_flow[row_tmp3]
      mean_spr_flow_tm2[i] <- mean_flow_stats$mean_spr_flow[row_tmp3]
      mean_sum_flow_tm2[i] <- mean_flow_stats$mean_sum_flow[row_tmp3]
      maxm_ann_flow_tm2[i] <- mean_flow_stats$max_ann_flow[row_tmp3]
      covn_ann_flow_tm2[i] <- mean_flow_stats$cov_ann_flow[row_tmp3]
      mean_spwn_flow_tm2[i] <- mean_flow_stats$mean_spwn_flow[row_tmp3]
      cov_spwn_flow_tm2[i] <- mean_flow_stats$cov_spwn_flow[row_tmp3]
    }
  }
  data$mannf <- mean_ann_flow_vec
  data$msprf <- mean_spr_flow_vec
  data$msumf <- mean_sum_flow_vec
  data$maxaf <- maxm_ann_flow_vec
  data$covaf <- covn_ann_flow_vec
  data$mspwn <- mean_spwn_flow_vec
  data$cspwn <- cov_spwn_flow_vec
  
  data$mannf_tm1 <- mean_ann_flow_tm1
  data$msprf_tm1 <- mean_spr_flow_tm1
  data$msumf_tm1 <- mean_sum_flow_tm1
  data$maxaf_tm1 <- maxm_ann_flow_tm1
  data$covaf_tm1 <- covn_ann_flow_tm1
  data$mspwn_tm1 <- mean_spwn_flow_tm1
  data$cspwn_tm1 <- cov_spwn_flow_tm1
  
  data$mannf_tm2 <- mean_ann_flow_tm2
  data$msprf_tm2 <- mean_spr_flow_tm2
  data$msumf_tm2 <- mean_sum_flow_tm2
  data$maxaf_tm2 <- maxm_ann_flow_tm2
  data$covaf_tm2 <- covn_ann_flow_tm2
  data$mspwn_tm2 <- mean_spwn_flow_tm2
  data$cspwn_tm2 <- cov_spwn_flow_tm2
  
  data
  
} 

subset_spp <- function(data, spp_subset) {
  
  # clean up NAs in reach column
  if (any(is.na(data$reach))) {
    data <- data[-which(is.na(data$reach)), ]
  }
  
  # subset species
  data <- data[which(!is.na(match(data$species, spp_subset))), ]
  
  data
  
}

load_otolith_data <- function(file_list) {
  
  otolith_data <- vector("list", length = length(file_list))
  for (i in seq_along(file_list))
    otolith_data[[i]] <- read.csv(paste0(file_list[i])) 
  
  names(otolith_data) <- file_list
  name_match <- match(c("./data/MC.csv", "./data/TC.csv",
                        "./data/YB.csv", "./data/SP.csv"),
                      file_list)
  names(otolith_data) <- c("murraycod", "troutcod", "goldenperch", "silverperch")[name_match]
  
  otolith_data <- t(sapply(otolith_data, function(x) coef_extract(filter_fun(x))))
  
  otolith_data
  
} 

add_age_data <- function(data, otolith_data) {
  
  # add length/weight to age calculations
  data$age <- rep(NA, nrow(data))
  for (i in seq_along(unique(data$species))) {
    
    row_sub <- which(data$species == unique(data$species)[i])
    if (any(rownames(otolith_data) == unique(data$species)[i])) {
      data$age[row_sub] <- length_to_age(data$length[row_sub],
                                         params = otolith_data[which(rownames(otolith_data) == unique(data$species)[i]), ])
    } else { 
      data$age[row_sub] <- rep(NA, length(row_sub))
    } 
    
  } 
  
  data
  
} 

create_catch_curves <- function(data) {
  
  sp_list <- unique(data$species)
  catch_curves <- vector("list", length = length(sp_list))
  for (i in seq_along(sp_list))
    catch_curves[[i]] <- catch_curve_fun(data, sp = sp_list[i])
  names(catch_curves) <- sp_list
  
  catch_curves
  
}

# prepare species x age matrices
prepare_analysis_data <- function(catch_curves, n_class = 6) {
  
  site_list <- sapply(catch_curves, function(x) apply(x$info, 1, paste, collapse = "_"))
  site_flat <- unique(unlist(site_list))
  
  n_sp <- length(catch_curves)
  
  age_out <- matrix(NA, nrow = length(site_flat), ncol = (n_class * n_sp))
  flow_out <- matrix(NA, nrow = length(site_flat), ncol = ncol(catch_curves[[1]]$flow))
  system <- reach <- year <- rep(NA, length(site_flat))
  for (i in seq_len(length(site_flat))) {
    
    row_ids <- lapply(site_list, function(x) which(x == site_flat[i]))
    if (any(sapply(row_ids, length) > 1)) {
      spp <- names(catch_curves)[sapply(row_ids) > 1]
      row_tmp <- row_ids[sapply(row_ids) > 1]
      stop(paste0("Observation ", row_tmp, " is being double counted for ", spp),
           call. = FALSE)
    }
    age_tmp <- rep(NA, ncol(age_out))
    flow_tmp <- rep(NA, ncol(flow_out))
    for (j in seq_along(catch_curves)) {
      if (length(row_ids[[j]])) {
        age_tmp[((j - 1) * n_class + 1):(j * n_class)] <- catch_curves[[j]]$age_dist[row_ids[[j]], seq_len(n_class)]
        system[i] <- as.character(catch_curves[[j]]$info$system[row_ids[[j]]])
        reach[i] <- catch_curves[[j]]$info$reach[row_ids[[j]]]
        year[i] <- catch_curves[[j]]$info$year[row_ids[[j]]]
        flow_tmp <- catch_curves[[j]]$flow[row_ids[[j]], ]
      }
    }
    age_out[i, ] <- age_tmp
    flow_out[i, ] <- flow_tmp
    
  }
  colnames(flow_out) <- colnames(catch_curves[[1]]$flow)
  bins_all <- catch_curves[[1]]$bins[seq_len(n_class)]
  info_out <- data.frame(system = system, reach = reach, year = year)
  
  # convert NAs to zeros
  age_out <- ifelse(is.na(age_out), 0, age_out)

  # collate output
  out <- list(age = age_out,
              flow = flow_out,
              info = info_out,
              bins = bins_all,
              n_class = n_class,
              n_sp = n_sp,
              n_obs = nrow(age_out))
  
  out
  
}
