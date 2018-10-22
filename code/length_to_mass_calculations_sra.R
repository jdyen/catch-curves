coefs_sra <- read.csv("./data/sra-length-weight-relation.csv")

# load data
alldat <- read.csv("./data/VEFMAP_FISH_20171024.csv")

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
alldat$Date <- dmy(alldat$Event_Date)

# pull out species names
sp_names <- unique(alldat$Common.Name)

# set up output list
length_weight_conversion_sra <- list()
for (i in seq_along(sp_names)) {
  
  # subset for species i
  coefs_row_id <- which(coefs_sra$best_guess_common == sp_names[i])
  if (length(coefs_row_id) > 1) {
    print(sp_names[i])
    coefs_row_id <- coefs_row_id[1]
  }

  # pull out coefs
  coefs_tmp <- NULL
  if (length(coefs_row_id)) {
    coefs_tmp <- c(intercept = coefs_sra$LW_constant[coefs_row_id],
                   slope = coefs_sra$LW_slope[coefs_row_id],
                   length_scale_fac = coefs_sra$length_factor[coefs_row_id],
                   fork_scale_fac = coefs_sra$fork_factor[coefs_row_id])
  }

  # save outputs
  length_weight_conversion_sra[[sp_names[i]]] <- coefs_tmp
  
}

length_weight_conversion_sra$generic <- apply(matrix(unlist(length_weight_conversion_sra),
                                                     ncol = 4, byrow = TRUE),
                                              2, mean)
names(length_weight_conversion_sra$generic) <- names(length_weight_conversion_sra[[1]])

rm(list = ls()[-grep("length_weight_conversion", ls())])