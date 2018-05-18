
# fill gaps in rainbows with averages for fish "identified" but not "caught" (double check this for all species)
### (use the gap filling approach from the ISD studies)
## ONLY have these data for the Ovens.

# load data
alldat <- read.csv("./data/VEFMAP_FISH_20171024.csv")

# load snags data set
snags_data <- read.csv("./data/SNAGS_FISH_20171205.csv")
snags_data$date_new <- format(dmy_hms(snags_data$surveydate), format = "%d/%m/%Y")
snags_data$YEAR <- sapply(strsplit(snags_data$date_new, "/"),
                          function(x) x[3])
snags_data$taxonname <- as.character(snags_data$taxonname)
snags_data$taxonname <- ifelse(snags_data$taxonname == "Yellowbelly",
                               "Golden perch",
                               snags_data$taxonname)
snags_data$taxonname <- factor(snags_data$taxonname)
snags_data2 <- data.frame(SYSTEM = rep("LOWERMURRAY", nrow(snags_data)),
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

# load ovens data and combine with alldat
ovens_data <- read.table("./data/vba_ovens_2008_2017.csv", sep = "\t", header = TRUE)
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
ovens_data$common_name <- alldat$Common.Name[match(ovens_data$species, alldat$Scientific.Name)]
ovens_data2 <- data.frame(SYSTEM = rep("OVENS", nrow(ovens_data)),
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
alldat <- rbind(alldat, ovens_data2, snags_data2)

# load flow data
# source("./code/load-flow-data.R")

# set systems of interest
system_sub <- c("BROKEN",
                "LOWERMURRAY",
                "CAMPASPE",
                "GLENELG",
                "GOULBURN",
                "LODDON",
                "THOMSON",
                "OVENS")
alldat <- alldat[-which(is.na(match(alldat$SYSTEM, system_sub))), ]

# clean up common names
alldat$Common.Name <- tolower(alldat$Common.Name)
alldat$Common.Name <- gsub(" ", "", alldat$Common.Name)
alldat$Common.Name <- gsub("-[[:digit:]]*", "", alldat$Common.Name)
alldat$Common.Name <- gsub("/", "", alldat$Common.Name)
alldat$Common.Name <- gsub("sp\\.", "", alldat$Common.Name)
alldat$Common.Name <- gsub("flatheadedgudgeon", "flatheadgudgeon", alldat$Common.Name)
alldat$Common.Name <- gsub("europeancarp", "carp", alldat$Common.Name)
alldat$Common.Name <- gsub("redfinperch", "redfin", alldat$Common.Name)
alldat$Common.Name <- ifelse(alldat$Common.Name == "weatherloach",
                             "orientalweatherloach",
                             alldat$Common.Name)
alldat$Common.Name <- ifelse(alldat$Common.Name == "rainbowfish",
                             "murrayriverrainbowfish",
                             alldat$Common.Name)
alldat$Common.Name <- ifelse(alldat$Common.Name == "hardyhead",
                             "unspeckedhardyhead",
                             alldat$Common.Name)
alldat$Common.Name <- ifelse(alldat$Common.Name == "blackfish",
                             "riverblackfish",
                             alldat$Common.Name)

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
alldat <- alldat[-which(!is.na(match(alldat$Common.Name, sp_to_rm))), ]
if (any(is.na(alldat$Common.Name)))
  alldat <- alldat[-which(is.na(alldat$Common.Name)), ]

# fill missing weights
sp_tmp <- unique(alldat$Common.Name)
for (i in seq_along(sp_tmp)) {
  
  # subset to single species
  dat_tmp <- alldat[which(alldat$Common.Name == sp_tmp[i]), ]
  
  # check to make sure some weights are missing      
  if (any(is.na(dat_tmp$WEIGHT))) {
    
    # use species-specific equation if it exists, generic otherwise
    if (length_weight_conversion[[sp_tmp[i]]]$n) {
      coefs <- length_weight_conversion[[sp_tmp[i]]]$coef
    } else {
      coefs <- length_weight_conversion$generic
    }
    
    # subset NA observations
    na_sub <- which(is.na(dat_tmp$WEIGHT))
    
    # estimate weight from length
    dat_tmp$WEIGHT[na_sub] <- exp(coefs[1] + coefs[2] * log(dat_tmp$totallength[na_sub]))
    
    # return estimated data to full data set
    alldat[which(alldat$Common.Name == sp_tmp[i]), ] <- dat_tmp
    
  }
}

# reformat dates
alldat$Date <- dmy(alldat$Event_Date)

# arrange alldat to have more consistent set of variables
alldat <- data.frame(date = alldat$Date,
                     year = alldat$YEAR,
                     site = alldat$SITE_CODE,
                     system = alldat$SYSTEM,
                     reach = alldat$Reach,
                     species = alldat$Common.Name,
                     length = alldat$totallength,
                     weight = alldat$WEIGHT,
                     abundance = alldat$Total.Sampled,
                     intensity = (alldat$total_no_passes * alldat$seconds))

alldat$reach_alt <- alldat$reach
alldat$reach_alt <- ifelse(alldat$system == "BROKEN",
                           ifelse(alldat$reach_alt == 4, 5, alldat$reach_alt), # downstream
                           alldat$reach_alt)
alldat$reach_alt <- ifelse(alldat$system == "THOMSON",
                           ifelse(alldat$reach_alt == 2, 3,             # downstream
                                  ifelse(alldat$reach_alt == 6, 5,      # upstream
                                         alldat$reach_alt)),
                           alldat$reach_alt)
alldat$reach_alt <- ifelse(alldat$system == "LODDON",
                           ifelse(alldat$reach_alt == 2, 3,             # downstream
                                  ifelse(alldat$reach_alt == 5, 4,      # upstream
                                         alldat$reach_alt)),
                           alldat$reach_alt)

mean_ann_flow_vec <- mean_spr_flow_vec <- mean_sum_flow_vec <- rep(NA, nrow(alldat))
maxm_ann_flow_vec <- covn_ann_flow_vec <- rep(NA, nrow(alldat))
mean_spwn_flow_vec <- cov_spwn_flow_vec <- rep(NA, nrow(alldat))
mean_spwn_flow_tm1 <- covn_ann_flow_tm1 <- cov_spwn_flow_tm1 <- rep(NA, nrow(alldat))
mean_spwn_flow_tm2 <- covn_ann_flow_tm2 <- cov_spwn_flow_tm2 <- rep(NA, nrow(alldat))
mean_ann_flow_tm1 <- mean_spr_flow_tm1 <- mean_sum_flow_tm1 <- maxm_ann_flow_tm1 <- rep(NA, nrow(alldat))
mean_ann_flow_tm2 <- mean_spr_flow_tm2 <- mean_sum_flow_tm2 <- maxm_ann_flow_tm2 <- rep(NA, nrow(alldat))
for (i in seq_len(nrow(alldat))) {
  row_tmp <- which((mean_flow_stats$system == alldat$system[i]) &
                     (mean_flow_stats$reach == alldat$reach[i]) &
                     (mean_flow_stats$year == alldat$year[i]))
  row_tmp2 <- which((mean_flow_stats$system == alldat$system[i]) &
                      (mean_flow_stats$reach == alldat$reach[i]) &
                      (mean_flow_stats$year == (alldat$year[i] - 1)))
  row_tmp3 <- which((mean_flow_stats$system == alldat$system[i]) &
                      (mean_flow_stats$reach == alldat$reach[i]) & 
                      (mean_flow_stats$year == (alldat$year[i] - 2)))
  if (!length(row_tmp)) {
    row_tmp <- which((mean_flow_stats$system == alldat$system[i]) &
                       (mean_flow_stats$reach == alldat$reach_alt[i]) &
                       (mean_flow_stats$year == alldat$year[i]))
  }
  if (!length(row_tmp2)) {
    row_tmp2 <- which((mean_flow_stats$system == alldat$system[i]) &
                        (mean_flow_stats$reach == alldat$reach_alt[i]) &
                        (mean_flow_stats$year == (alldat$year[i] - 1)))
  }
  if (!length(row_tmp3)) {
    row_tmp3 <- which((mean_flow_stats$system == alldat$system[i]) &
                        (mean_flow_stats$reach == alldat$reach_alt[i]) & 
                        (mean_flow_stats$year == (alldat$year[i] - 2)))
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
alldat$mannf <- mean_ann_flow_vec
alldat$msprf <- mean_spr_flow_vec
alldat$msumf <- mean_sum_flow_vec
alldat$maxaf <- maxm_ann_flow_vec
alldat$covaf <- covn_ann_flow_vec
alldat$mspwn <- mean_spwn_flow_vec
alldat$cspwn <- cov_spwn_flow_vec

alldat$mannf_tm1 <- mean_ann_flow_tm1
alldat$msprf_tm1 <- mean_spr_flow_tm1
alldat$msumf_tm1 <- mean_sum_flow_tm1
alldat$maxaf_tm1 <- maxm_ann_flow_tm1
alldat$covaf_tm1 <- covn_ann_flow_tm1
alldat$mspwn_tm1 <- mean_spwn_flow_tm1
alldat$cspwn_tm1 <- cov_spwn_flow_tm1

alldat$mannf_tm2 <- mean_ann_flow_tm2
alldat$msprf_tm2 <- mean_spr_flow_tm2
alldat$msumf_tm2 <- mean_sum_flow_tm2
alldat$maxaf_tm2 <- maxm_ann_flow_tm2
alldat$covaf_tm2 <- covn_ann_flow_tm2
alldat$mspwn_tm2 <- mean_spwn_flow_tm2
alldat$cspwn_tm2 <- cov_spwn_flow_tm2

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

# clean up NAs in reach column
if (any(is.na(alldat$reach))) {
  alldat <- alldat[-which(is.na(alldat$reach)), ]
}

# subset species
sp_sub <- c("goldenperch", "murraycod",
            "troutcod", "silverperch")
alldat <- alldat[which(!is.na(match(alldat$species, sp_sub))), ]

# add length/weight to age calculations
alldat$age <- rep(NA, nrow(alldat))
for (i in seq_along(unique(alldat$species))) {
  
  row_sub <- which(alldat$species == unique(alldat$species)[i])
  if (any(rownames(params_all) == unique(alldat$species)[i])) {
    alldat$age[row_sub] <- length_to_age(alldat$length[row_sub],
                                         params = params_all[which(rownames(params_all) == unique(alldat$species)[i]), ])
  } else {
    alldat$age[row_sub] <- rep(NA, length(row_sub))
  }
  
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

mc_catch_curve <- catch_curve_fun(alldat, sp = "murraycod")
gp_catch_curve <- catch_curve_fun(alldat, sp = "goldenperch")
sp_catch_curve <- catch_curve_fun(alldat, sp = "silverperch")
tc_catch_curve <- catch_curve_fun(alldat, sp = "troutcod")

rm(catch_curve_fun, coef_extract, coefs, cov_spwn_flow_tm1, cov_spwn_flow_tm2,
   cov_spwn_flow_tmp, cov_spwn_flow_vec, covn_ann_flow_tm1, covn_ann_flow_tm2,
   covn_ann_flow_vec, dat_tmp, filter_fun, i, length_to_age,
   length_weight_conversion, length_weight_conversion_sra, max_fun,
   maxm_ann_flow_tm1, maxm_ann_flow_tm2, maxm_ann_flow_vec,
   mean_ann_flow_tm1, mean_ann_flow_tm2, mean_ann_flow_vec,
   mean_flow_stats, mean_spr_flow_tm1, mean_spr_flow_tm2,
   mean_spr_flow_vec, mean_spwn_flow_tm1, mean_spwn_flow_tm2,
   mean_spwn_flow_tmp, mean_spwn_flow_vec, mean_sum_flow_tm1,
   mean_sum_flow_tm2, mean_sum_flow_vec, na_sub, oti_data,
   ovens_data, ovens_data2, params_all, predictors,
   reach_no, row_sub, row_tmp, row_tmp2, row_tmp3,
   size_dist_fun, snags_data, snags_data2, sp_sub, sp_tmp, sp_to_rm,
   spawning_flow_data, subset, sys_name, system_sub)
