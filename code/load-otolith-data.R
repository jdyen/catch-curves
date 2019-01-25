# load otolith data for all species
all_oti <- read.csv("data/MurrayNtribs_otolith_data.csv", stringsAsFactors = FALSE)
all_oti$X.SPECIES[all_oti$X.SPECIES == "MC "] <- "MC"
all_oti$X.SPECIES[all_oti$X.SPECIES == "mc"] <- "MC"
all_oti$X.SPECIES[all_oti$X.SPECIES == "silver perch"] <- "Silver perch"
all_oti$X.SPECIES[all_oti$X.SPECIES == "Silver Perch"] <- "Silver perch"
all_oti$X.SPECIES[all_oti$X.SPECIES == "Sp"] <- "SP"
all_oti$X.SPECIES[all_oti$X.SPECIES == "Tc"] <- "TC"
all_oti <- all_oti[all_oti$X.SPECIES != "BT", ]

switch_names <- function(x) {
  
  out <- rep(NA, length(x))
  for (i in seq_along(x)) {
    out[i] <- 
      switch(x[i],
             "MC" = "Maccullochella peelii",
             "Murray Cod" = "Maccullochella peelii",
             "SP" = "Bidyanus bidyanus",
             "Silver perch" = "Bidyanus bidyanus",
             "TC" = "Maccullochella macquariensis",
             "YB" = "Macquaria ambigua",
             "GP" = "Macquaria ambigua",
             "Golden Perch" = "Macquaria ambigua",
             "RF" = "Melanotaenia fluviatilis")
  }
  
  out
  
}

all_oti$SPECIES <- switch_names(all_oti$X.SPECIES)

# ages and lengths are character valued for some reason
all_oti$T_Length..mm. <- as.numeric(all_oti$T_Length..mm.)
all_oti$AGE <- as.numeric(all_oti$AGE)

saveRDS(all_oti, file = "data/otolith-data-loaded.rds")
