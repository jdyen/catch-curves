# R code to relate weight to length in MDB fishes

files_cur <- ls()

# load data
alldat <- read.csv("../catch-curves/data/VEFMAP_ALL_FISH_190731.csv")

# load snags data set
snags_data <- read.csv("../catch-curves/data/SNAGS_FISH_20190827.csv")
snags_data$date_new <- parse_date_time(snags_data$event_date, orders = c("dmy_HM"))
snags_data$YEAR <- year(snags_data$date_new)
snags_data$taxonname <- as.character(snags_data$common_name)
snags_data$taxonname <- ifelse(snags_data$taxonname == "Yellowbelly",
                               "Golden perch",
                               snags_data$taxonname)
snags_data$taxonname <- factor(snags_data$taxonname)
snags_data2 <- data.frame(sample_id = rep(NA, nrow(snags_data)),
                          system = rep("LOWERMURRAY", nrow(snags_data)),
                          site = paste0("Lm", snags_data$site),
                          reach = rep(1, nrow(snags_data)),
                          gear_type = factor(rep("EF/Boat"), nrow(snags_data)),
                          event_date = snags_data$date_new,
                          pass_no = rep(1, nrow(snags_data)),
                          total_no_passes = rep(1, nrow(snags_data)),
                          ef_seconds_total = snags_data$ef_seconds_total,
                          taxon_id = rep(NA, nrow(snags_data)),
                          common_name = snags_data$taxonname,
                          scientific_name = snags_data$scientific_name,
                          total_length_mm = snags_data$total_length_mm,
                          weight_g = snags_data$weight_g,
                          no_collected = rep(1, nrow(snags_data)),
                          no_observed = rep(NA, nrow(snags_data)),
                          vefmap_stage = rep(NA, nrow(snags_data)),
                          year = as.integer(snags_data$YEAR),
                          notes = rep(NA, nrow(snags_data)))

# load ovens data and combine with alldat
ovens_data <- read.csv("../catch-curves/data/OVENS_VBA_2008_2017.csv", header = TRUE)
ovens_data$date_new <- parse_date_time(ovens_data$event_date, orders = c("dmy"))
ovens_data$YEAR <- year(ovens_data$date_new)
ovens_data$species <- as.character(ovens_data$scientific_name)
ovens_data$species <- ifelse(ovens_data$species == "Maccullochella peelii ",
                             "Maccullochella peelii peelii",
                             ovens_data$species)
ovens_data$species <- ifelse(ovens_data$species == "Maccullochella peelii",
                             "Maccullochella peelii peelii",
                             ovens_data$species)
ovens_data$species <- factor(ovens_data$species)
ovens_data$common_name <- alldat$common_name[match(ovens_data$species, alldat$scientific_name)]
ovens_data2 <- data.frame(sample_id = rep(NA, nrow(ovens_data)),
                          system = rep("OVENS", nrow(ovens_data)),
                          site = paste0("Ov", ovens_data$site),
                          reach = rep(1, nrow(ovens_data)),
                          gear_type = ovens_data$gear_type,
                          event_date = ovens_data$date_new,
                          pass_no = rep(1, nrow(ovens_data)),
                          total_no_passes = rep(1, nrow(ovens_data)),
                          ef_seconds_total = ovens_data$ef_seconds_total,
                          taxon_id = rep(NA, nrow(ovens_data)),
                          common_name = ovens_data$common_name,
                          scientific_name = ovens_data$species,
                          total_length_mm = ovens_data$total_length_mm,
                          weight_g = ovens_data$weight_g,
                          no_collected = ovens_data$no_collected,
                          no_observed = rep(NA, nrow(ovens_data)),
                          vefmap_stage = rep(NA, nrow(ovens_data)),
                          year = as.integer(ovens_data$YEAR),
                          notes = rep(NA, nrow(ovens_data)))

alldat <- rbind(alldat, ovens_data2, snags_data2)

# set systems of interest
system_sub <- c("BROKEN",
                "LOWERMURRAY",
                "CAMPASPE",
                "GLENELG",
                "GOULBURN",
                "LODDON",
                "THOMSON",
                "OVENS")
alldat <- alldat[-which(is.na(match(alldat$system, system_sub))), ]

# clean up common names
alldat$Common.Name <- tolower(alldat$common_name)
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
alldat$Common.Name <- gsub("carp", "commoncarp", alldat$Common.Name)

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
alldat <- alldat[-which(is.na(alldat$Common.Name)), ]

# reformat dates
alldat$Date <- alldat$event_date

# pull out species names
sp_names <- unique(alldat$Common.Name)

# set up output list
length_weight_conversion <- list()

# loop through each species
for (i in seq_along(sp_names)) {
  
  # subset for species i
  dat <- alldat[which(alldat$Common.Name == sp_names[i]), ]
  
  # log transform length and weight
  rows_to_rm <- NULL
  if (any(is.na(dat$total_length_mm)))
    rows_to_rm <- c(rows_to_rm, which(is.na(dat$total_length_mm)))
  if (any(is.na(dat$weight_g)))
    rows_to_rm <- c(rows_to_rm, which(is.na(dat$weight_g)))
  if (any(dat$total_length_mm <= 0, na.rm = TRUE))
    rows_to_rm <- c(rows_to_rm, which(dat$total_length_mm <= 0))
  if (any(dat$weight_g <= 0, na.rm = TRUE))
    rows_to_rm <- c(rows_to_rm, which(dat$weight_g <= 0))
  if (length(rows_to_rm))
    dat <- dat[-rows_to_rm, ]
  log_length <- log(dat$total_length_mm)
  log_weight <- log(dat$weight_g)
  
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

rm(list = ls()[-match(c(files_cur, "length_weight_conversion"), ls())])
