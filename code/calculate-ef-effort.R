
## send Zeb effort data for all sites, especially non-VEFMAP sites. Format in SITE x year table?
##  OR just match EF audit file format.

# set working directory
setwd("~/Dropbox/research/catch-curves/")

# load pacakges
library(lubridate)

# load compiled survey data
alldat <- readRDS("data/data-loaded-Jul19.rds")
to_keep <- alldat$Scientific.Name == "Maccullochella peelii"
alldat <- alldat[to_keep, ]
alldat$Common.Name <- NULL

# pull out the effort data
seconds_tmp <- tapply(alldat$seconds, list(alldat$SITE_CODE, alldat$YEAR), unique)
passes_tmp <- tapply(alldat$total_no_passes, list(alldat$SITE_CODE, alldat$YEAR), unique)
n_vals <- sapply(seconds_tmp, length)
sites_expanded <- rep(rownames(seconds_tmp), times = ncol(seconds_tmp))
years_expanded <- rep(colnames(seconds_tmp), each = nrow(seconds_tmp))
seconds_data <- data.frame(sites = rep(sites_expanded, times = n_vals),
                           years = rep(years_expanded, times = n_vals),
                           seconds = unlist(seconds_tmp))
seconds_data$system <- alldat$SYSTEM[match(seconds_data$sites, alldat$SITE_CODE)]
seconds_data$reach <- alldat$Reach[match(seconds_data$sites, alldat$SITE_CODE)]
seconds_data$gear_type <- alldat$geartype[match(seconds_data$sites, alldat$SITE_CODE)]

write.csv(seconds_data, file = "data/ef_effort_data.csv", row.names = FALSE)
