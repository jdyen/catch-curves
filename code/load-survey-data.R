# load data
vefmap17 <- read.csv("./data/VEFMAP_FISH_20171024.csv")
vefmap17new <- read.csv("./data/VEFMAP_FISH_180726.csv")

# convert date to string
vefmap17new$Event_Date <- dmy_hms(vefmap17new$Event_Date)
vefmap17new$SYSTEM <- gsub("GOULBURN ", "GOULBURN", vefmap17new$SYSTEM)
vefmap17$Event_Date <- dmy(vefmap17$Event_Date)

# remove unneeded columns from vefmap17new
colnames(vefmap17new)[c(4, 9)] <- c("Reach", "seconds")
vefmap17new <- vefmap17new[, match(colnames(vefmap17), colnames(vefmap17new))]

# add in a recorded column to both data sets
vefmap17$recorded <- rep(1, nrow(vefmap17))
vefmap17new$recorded <- rep(1, nrow(vefmap17new))

# load MC conversion from length to mass
mc_conv <- read.csv("./data/murray_cod_length_to_mass.csv")
colnames(mc_conv) <- c("length_cm", "weight_kg")
mc_conv$length_mm <- mc_conv$length_cm * 10
mc_conv$weight_g <- mc_conv$weight_kg * 1000

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
                          Event_Date = dmy(snags_data$date_new),
                          Pass.No = rep(1, nrow(snags_data)),
                          total_no_passes = rep(1, nrow(snags_data)),
                          seconds = snags_data$seconds,
                          Common.Name = snags_data$taxonname,
                          Scientific.Name = snags_data$Scientific.Name,
                          totallength = snags_data$totallength,
                          WEIGHT = snags_data$weight,
                          Total.Sampled = rep(1, nrow(snags_data)),
                          VEFMAP.Stage = rep(NA, nrow(snags_data)),
                          YEAR = as.integer(snags_data$YEAR),
                          recorded = rep(1, nrow(snags_data)))

# load ovens data and combine with vefmap17
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
ovens_data$common_name <- vefmap17$Common.Name[match(ovens_data$species, vefmap17$Scientific.Name)]
ovens_data2 <- data.frame(SYSTEM = rep("OVENS", nrow(ovens_data)),
                          SITE_CODE = paste0("Ov", ovens_data$site),
                          Reach = rep(1, nrow(ovens_data)),
                          geartype = ovens_data$gear_type,
                          Event_Date = dmy(ovens_data$date_new),
                          Pass.No = rep(1, nrow(ovens_data)),
                          total_no_passes = rep(1, nrow(ovens_data)),
                          seconds = ovens_data$electro_seconds,
                          Common.Name = ovens_data$common_name,
                          Scientific.Name = ovens_data$species,
                          totallength = ovens_data$total_length_mm,
                          WEIGHT = ovens_data$weight_g,
                          Total.Sampled = ovens_data$no_collected,
                          VEFMAP.Stage = rep(NA, nrow(ovens_data)),
                          YEAR = as.integer(ovens_data$YEAR),
                          recorded = ovens_data$no_obs)

# load the 2018 vefmap data and combine with alldat
vefmap18 <- read.csv("./data/vefmap-2018.csv")
vefmap18_clean <- data.frame(SYSTEM = system_switch_fun(vefmap18$SITE.NUMBER),
                             SITE_CODE = vefmap18$SITE.NUMBER,
                             Reach = vefmap17$Reach[match(vefmap18$SITE.NUMBER, vefmap17$SITE_CODE)],
                             geartype = vefmap18$Gear.Type,
                             Event_Date = dmy(vefmap18$Date.Start),
                             Pass.No = vefmap18$Electro.fishing.Pass.No,
                             total_no_passes = vefmap18$Total.No.Passes,
                             seconds = vefmap18$Duration,
                             Common.Name = vefmap17$Common.Name[match(vefmap18$Scientific.Name,
                                                                      vefmap17$Scientific.Name)],
                             Scientific.Name = vefmap18$Scientific.Name,
                             totallength = vefmap18$TOTAL.LENGTH,
                             WEIGHT = vefmap18$WEIGHT,
                             Total.Sampled = vefmap18$Total.collected.for.each.particular.gear,
                             VEFMAP.Stage = rep(6, nrow(vefmap18)),
                             YEAR = year(dmy(vefmap18$Date.Start)),
                             recorded = vefmap18$No.Observed)
vefmap18_clean$Total.Sampled <- ifelse(as.character(vefmap18_clean$Scientific.Name) == "Melanotaenia fluviatilis",
                                       ifelse(is.na(vefmap18_clean$Total.Sampled),
                                              vefmap18$No.Observed,
                                              vefmap18_clean$Total.Sampled),
                                       vefmap18_clean$Total.Sampled)
vefmap18_clean$Common.Name <- ifelse(vefmap18_clean$Scientific.Name == "Maccullochella peelii",
                                     "Murray cod", vefmap18_clean$Common.Name)
vefmap18_clean <- vefmap18_clean[vefmap18_clean$Scientific.Name != "Misc No Fish", ]
vefmap18_clean <- vefmap18_clean[vefmap18_clean$Scientific.Name != "Misc No fish", ]

# add newest Goulburn data
vefmap18_additional <- read.csv("./data/Stage5DatawithWKGoulburndataAdded_20180710.csv")
colnames(vefmap18_additional)[1] <- "SYSTEM"
vefmap18_additional_clean <- data.frame(SYSTEM = vefmap18_additional$SYSTEM,
                             SITE_CODE = vefmap18_additional$SITE_CODE,
                             Reach = vefmap18_additional$REACH,
                             geartype = vefmap18_additional$SAMP_TECH,
                             Event_Date = dmy(vefmap18_additional$Date),
                             Pass.No = vefmap18_additional$SAMP_NUMBER,
                             total_no_passes = vefmap18_additional$SAMP_UNITS,
                             seconds = vefmap18_additional$SAMP_TIME,
                             Common.Name = vefmap18_additional$Common.Name,
                             Scientific.Name = vefmap18_additional$Scientific.Name,
                             totallength = vefmap18_additional$LENGTH,
                             WEIGHT = vefmap18_additional$WEIGHT,
                             Total.Sampled = rep(1, nrow(vefmap18_additional)),
                             VEFMAP.Stage = rep(6, nrow(vefmap18_additional)),
                             YEAR = year(dmy(vefmap18_additional$Date)),
                             recorded = vefmap18_additional$ABUND)
vefmap18_additional_clean$Scientific.Name <- gsub("Melanotaenia fluviatilis ",
                                                  "Melanotaenia fluviatilis",
                                                  vefmap18_additional_clean$Scientific.Name)
vefmap18_additional_clean$Scientific.Name <- gsub("Bidyanus bidyanus ",
                                                  "Bidyanus bidyanus",
                                                  vefmap18_additional_clean$Scientific.Name)
vefmap18_additional_clean$Scientific.Name <- gsub("Macquaria ambigua ",
                                                  "Macquaria ambigua",
                                                  vefmap18_additional_clean$Scientific.Name)
vefmap18_additional_clean$Scientific.Name <- gsub("philypnodon grandiceps",
                                                  "Philypnodon grandiceps",
                                                  vefmap18_additional_clean$Scientific.Name)
vefmap18_additional_clean$SYSTEM <- gsub("BROKEN ",
                                         "BROKEN",
                                         vefmap18_additional_clean$SYSTEM)
vefmap18_additional_clean$SYSTEM <- gsub("GOULBURN ",
                                         "GOULBURN",
                                         vefmap18_additional_clean$SYSTEM)

# load king river data
king_data <- read.csv("data/king-2009-2016.csv", stringsAsFactors = FALSE)
king_data$SPECIES <- ifelse(king_data$SPECIES == "Maccullochella peelii ",
                            "Maccullochella peelii peelii",
                            king_data$SPECIES)
king_data$SPECIES <- ifelse(king_data$SPECIES == "Maccullochella peelii",
                            "Maccullochella peelii peelii",
                            king_data$SPECIES)
king_data$SPECIES <- ifelse(king_data$SPECIES == "Maccullochella peelii peelii ",
                            "Maccullochella peelii peelii",
                            king_data$SPECIES)
king_data$SPECIES <- ifelse(king_data$SPECIES == "maccullochella macquariensis",
                            "Maccullochella macquariensis",
                            king_data$SPECIES)
king_data$Sampled <- apply(cbind(king_data$No..Coll., king_data$No..obs.), 1, sum, na.rm = TRUE)
king_clean <- data.frame(SYSTEM = rep("KING", nrow(king_data)),
                         SITE_CODE = king_data$Site.No.,
                         Reach = rep(NA, nrow(king_data)),
                         geartype = king_data$Gear.Type..EF.BP..EF.LB..EF.MB..EF.SB..BT.,
                         Event_Date = dmy(king_data$Date),
                         Pass.No = king_data$Operation,
                         total_no_passes = rep(NA, nrow(king_data)),
                         seconds = king_data$Electro.Seconds,
                         Common.Name = rep(NA, nrow(king_data)),
                         Scientific.Name = king_data$SPECIES,
                         totallength = king_data$Total.Length..mm.,
                         WEIGHT = king_data$Weight..gms.,
                         Total.Sampled = king_data$No..Coll.,
                         VEFMAP.Stage = rep(NA, nrow(king_data)),
                         YEAR = year(dmy(king_data$Date)),
                         recorded = king_data$Sampled)

# load mid-ovens data
mid_ovens_data <- read.csv("data/mid-ovens-2012-2018-data.csv", stringsAsFactors = FALSE)
mid_ovens_data$SPECIES <- ifelse(mid_ovens_data$SPECIES == "Maccullochella macquariensis ",
                                 "Maccullochella macquariensis",
                                 mid_ovens_data$SPECIES)
mid_ovens_data$SPECIES <- ifelse(mid_ovens_data$SPECIES == "Maccullochella peelii ",
                                 "Maccullochella peelii",
                                 mid_ovens_data$SPECIES)
mid_ovens_data$SPECIES <- ifelse(mid_ovens_data$SPECIES == "cyprinus carpio",
                                 "Cyprinus carpio",
                                 mid_ovens_data$SPECIES)
mid_ovens_data$SPECIES <- ifelse(mid_ovens_data$SPECIES == "Cyprinus carpio ",
                                 "Cyprinus carpio",
                                 mid_ovens_data$SPECIES)
mid_ovens_data$Sampled <- apply(cbind(mid_ovens_data$Collected, mid_ovens_data$Observed), 1, sum, na.rm = TRUE)
mid_ovens_clean <- data.frame(SYSTEM = rep("MIDOVENS", nrow(mid_ovens_data)),
                              SITE_CODE = mid_ovens_data$Location,
                              Reach = rep(NA, nrow(mid_ovens_data)),
                              geartype = mid_ovens_data$Gear.Type,
                              Event_Date = dmy(mid_ovens_data$Date..DD.MM.YYYY.),
                              Pass.No = mid_ovens_data$Shot,
                              total_no_passes = rep(NA, nrow(mid_ovens_data)),
                              seconds = mid_ovens_data$Electro.Seconds,
                              Common.Name = rep(NA, nrow(mid_ovens_data)),
                              Scientific.Name = mid_ovens_data$SPECIES,
                              totallength = mid_ovens_data$length_mm,
                              WEIGHT = mid_ovens_data$weight_g,
                              Total.Sampled = mid_ovens_data$Collected,
                              VEFMAP.Stage = rep(NA, nrow(mid_ovens_data)),
                              YEAR = year(dmy(mid_ovens_data$Date)),
                              recorded = mid_ovens_data$Sampled)

# load ovens 2018 data
ovens2018_data <- read.csv("data/ovens-2018-data.csv", stringsAsFactors = FALSE)
ovens2018_data$Sampled <- apply(cbind(ovens2018_data$Total.number.of.each.species.collected, ovens2018_data$No..Observed), 1, sum, na.rm = TRUE)
ovens2018_clean <- data.frame(SYSTEM = rep("OVENS", nrow(ovens2018_data)),
                             SITE_CODE = ovens2018_data$Site..,
                             Reach = rep(NA, nrow(ovens2018_data)),
                             geartype = ovens2018_data$Gear.Type,
                             Event_Date = dmy(ovens2018_data$Date.Start..DD.MM.YYYY.),
                             Pass.No = ovens2018_data$Shot.Number,
                             total_no_passes = rep(NA, nrow(ovens2018_data)),
                             seconds = ovens2018_data$Electro.Time.On..secs.,
                             Common.Name = rep(NA, nrow(ovens2018_data)),
                             Scientific.Name = ovens2018_data$SPECIES..scientific.names.,
                             totallength = ovens2018_data$LENGTH..mm.....for.each.individual.animal..,
                             WEIGHT = ovens2018_data$WEIGHT..gms.....for.each.individual.animal.,
                             Total.Sampled = ovens2018_data$Total.number.of.each.species.collected,
                             VEFMAP.Stage = rep(NA, nrow(ovens2018_data)),
                             YEAR = year(dmy(ovens2018_data$Date.Start..DD.MM.YYYY.)),
                             recorded = ovens2018_data$Sampled)

# combine all data sets
alldat <- rbind(vefmap17new, ovens_data2, snags_data2, vefmap18_clean, vefmap18_additional_clean,
                ovens2018_clean, mid_ovens_clean)
alldat$dataset <- c(rep("vefmap", nrow(vefmap17new)),
                    rep("ovens", nrow(ovens_data2)),
                    rep("snags", nrow(snags_data2)),
                    rep("vefmap", nrow(vefmap18_clean)),
                    rep("vefmap", nrow(vefmap18_additional_clean)),
                    rep("ovens", nrow(ovens2018_clean)),
                    rep("ovens", nrow(mid_ovens_clean)))

# add reach IDs for odd sites
alldat$Reach[grep("Th4", alldat$SITE_CODE)] <- 4
alldat$Reach[grep("TH4", alldat$SITE_CODE)] <- 4
alldat$Reach[grep("LO3", alldat$SITE_CODE)] <- 3
alldat$Reach[grep("Lo3", alldat$SITE_CODE)] <- 3
alldat$Reach[grep("LMR", alldat$SITE_CODE)] <- 1
alldat$Reach[grep("PCK0", alldat$SITE_CODE)] <- 1
alldat$Reach[grep("TC01", alldat$SITE_CODE)] <- 1
alldat$Reach[grep("Cable Hole", alldat$SITE_CODE)] <- 1

alldat$SITE_CODE <- gsub("Cable Hole", "GOCH1", alldat$SITE_CODE)
alldat$SITE_CODE <- gsub("Lo", "LO", alldat$SITE_CODE)
alldat$SITE_CODE <- gsub("Ca", "CA", alldat$SITE_CODE)
alldat$SITE_CODE <- gsub("Gl", "GL", alldat$SITE_CODE)
alldat$SITE_CODE <- gsub("Lo", "LO", alldat$SITE_CODE)
alldat$SITE_CODE <- gsub("Lm", "LM", alldat$SITE_CODE)
alldat$SITE_CODE <- gsub("Th", "TH", alldat$SITE_CODE)
alldat$SITE_CODE <- gsub("Ov", "OV", alldat$SITE_CODE)
alldat$SITE_CODE <- gsub(" ", "", alldat$SITE_CODE)

alldat$Scientific.Name <- as.character(alldat$Scientific.Name)
alldat$Scientific.Name <- ifelse(alldat$Scientific.Name == "Maccullochella peelii ",
                                 "Maccullochella peelii",
                                 alldat$Scientific.Name)
alldat$Scientific.Name <- ifelse(alldat$Scientific.Name == "Maccullochella peelii  ",
                                 "Maccullochella peelii",
                                 alldat$Scientific.Name) 
alldat$Scientific.Name <- ifelse(alldat$Scientific.Name == "Maccullochella peelii peelii",
                                 "Maccullochella peelii", 
                                 alldat$Scientific.Name)
alldat$Scientific.Name <- ifelse(alldat$Scientific.Name == "Maccullochella macquariensis ",
                                 "Maccullochella macquariensis",
                                 alldat$Scientific.Name)
alldat$Scientific.Name <- ifelse(alldat$Scientific.Name == "Melanotaenia fluviatilis ",
                                 "Melanotaenia fluviatilis",
                                 alldat$Scientific.Name)
alldat$Scientific.Name <- ifelse(alldat$Scientific.Name == "Macquaria ambigua ",
                                 "Macquaria ambigua", 
                                 alldat$Scientific.Name)
alldat$Scientific.Name <- ifelse(alldat$Scientific.Name == "bidyanus bidyanus",
                                 "Bidyanus bidyanus", 
                                 alldat$Scientific.Name)
alldat$Scientific.Name <- ifelse(alldat$Scientific.Name == "Bidyanus bidyanus ",
                                 "Bidyanus bidyanus",
                                 alldat$Scientific.Name) 
alldat$Scientific.Name <- ifelse(alldat$Scientific.Name == "Cyprinus carpio  ",
                                 "Cyprinus carpio", 
                                 alldat$Scientific.Name)
alldat$Scientific.Name <- ifelse(alldat$Scientific.Name == "Macquaria australasica ",
                                 "Macquaria australasica", 
                                 alldat$Scientific.Name)

saveRDS(alldat, file = "data/data-loaded.rds")
