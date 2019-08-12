# set working directory
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(lubridate)
library(greta)

# load some helper functions
source("code/helpers.R")
source("code/fit_ccr.R")
source("code/methods.R")
source("code/validate_ccr.R")

## RELOAD DATA -- IS OVENS VBA DATE FORMATTED CORRECTLY?
## SOME YEARS MISSING IN MID-OVENS SURVEYS?? CHECK LATEST FILES

# load compiled survey data
alldat <- readRDS("data/data-loaded-Aug19.rds")

# filter survey data to MC
to_keep <- alldat$scientific_name == "Maccullochella peelii"
alldat <- alldat[to_keep, ]

# load otolith data
oti_data <- readRDS("data/otolith-data-loaded.rds")

# filter to MC
oti_data <- oti_data[oti_data$SPECIES == "Maccullochella peelii", ]

# need to load flow data
flow_data <- readRDS("data/flow-data-loaded-Aug19.rds")

# filter to MC
flow_data <- flow_data[to_keep, ]

# optional: subset to systems of interest
systems_to_keep <- c("broken", "goulburn", "king", "murray", "ovens")
flow_data <- flow_data[alldat$system %in% systems_to_keep, ]
alldat <- alldat[alldat$system %in% systems_to_keep, ]

# hack for now because ovens site data incomplete
alldat$site[is.na(alldat$site)] <- 1

# prepare survey data
survey_data <- data.frame(length_mm = alldat$total_length_mm,
                          system = alldat$system,
                          system_id = rebase_index(alldat$system),
                          site = alldat$site,
                          site_id = rebase_index(alldat$site),
                          date = alldat$date_formatted,
                          dataset = alldat$dataset,
                          dataset_id = rebase_index(alldat$dataset),
                          effort = alldat$ef_seconds_total,
                          stringsAsFactors = FALSE)

# add years
survey_data$year <- year(survey_data$date)
survey_data$year_id <- rebase_index(survey_data$year)

flow_data <- flow_data[apply(survey_data, 1, function(x) !any(is.na(x))), ]
survey_data <- survey_data[apply(survey_data, 1, function(x) !any(is.na(x))), ]

# bin otolith data by age and size: offset by -0.4 (0-0.6 = YOY, 0.6-1.6 = 1YO, etc.)
oti_data$age_class <- cut(oti_data$AGE, breaks = c(-0.4:ceiling(max(oti_data$AGE, na.rm = TRUE))),
                          labels = FALSE)
size_breaks <- c(0, 150, 276, 469, 544, 610, 667, 718, 760, max(survey_data$length_mm, na.rm = TRUE))
oti_data$len_class <- cut(oti_data$T_Length..mm., breaks = size_breaks, labels = FALSE)
length_age_matrix <- classify(oti_data$len_class, oti_data$age_class)

# how many length classes do we need to keep to give all 0-4 year olds?
n_len <- max(which(length_age_matrix[, 5] > 0))

# how many age classes do we need to keep to give all 0:n_len length individuals?
n_age <- max(which(length_age_matrix[n_len, ] > 0))

# let's just keep those from the age_length_matrix
length_age_matrix <- length_age_matrix[seq_len(n_len), seq_len(n_age)]

# need to bin the survey data by lengths
response_matrix <- do.call(
  rbind, tapply(survey_data$length_mm,
                list(survey_data$system, survey_data$year),
                hist_fn, breaks = size_breaks))

# filter length class survey data to the length classes we care about (1:n_len)
response_matrix <- response_matrix[, seq_len(n_len)]

# define system/year predictors
system <- c(tapply(survey_data$system_id,
                   list(survey_data$system, survey_data$year),
                   unique))
system <- system[!is.na(system)]
year <- c(tapply(survey_data$year_id,
                 list(survey_data$system, survey_data$year),
                 unique))
year <- year[!is.na(year)]

# how much effort went into each survey? This should be summed over all sits in a given system
#   (but not double-count effort for each individual fish)
sys_site <- paste(survey_data$system, survey_data$site, sep = "_")
effort_site <- tapply(survey_data$effort,
                      list(sys_site, survey_data$year),
                      unique)
year_ids <- rep(colnames(effort_site), each = nrow(effort_site))
system_ids <- sapply(strsplit(rownames(effort_site), "_"), function(x) x[1])
system_ids <- rep(system_ids, times = ncol(effort_site))
effort_site <- sapply(c(effort_site), sum)
effort <- tapply(effort_site, list(system_ids, year_ids), sum)
effort <- effort[effort > 0]

# add a quadratic summer effect to the flow predictors
flow_data$prop_sum_win_sq <- flow_data$prop_sum_lt_win ^ 2
flow_data$prop_sum_win_sq_ym1 <- flow_data$prop_sum_lt_win_ym1 ^ 2

# compile flow predictors
var_sets <- list(c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "prop_sum_lt_win_ym1", "prop_sum_win_sq", "maxan_mld", "maxan_mld_ym1", "spwntmp_c"),
                 c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "prop_sum_win_sq", "maxan_mld", "maxan_mld_ym1", "spwntmp_c"),
                 c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "maxan_mld", "maxan_mld_ym1", "spwntmp_c"),
                 c("rrang_spwn_mld", "prop_spr_lt_win", "prop_sum_lt_win", "maxan_mld", "spwntmp_c"),
                 c("rrang_spwn_mld", "spwntmp_c"),
                 c("prop_spr_lt_win", "spwntmp_c"),
                 c("prop_sum_lt_win", "spwntmp_c"),
                 c("maxan_mld", "spwntmp_c"))
mod <- list()
mod_cv <- list()
for (i in seq_along(var_sets)) {
  
  vars_to_include <- var_sets[[i]]
  flow_compiled <- sapply(vars_to_include,
                          function(x) tapply(get(x, flow_data),
                                             list(survey_data$system_id, survey_data$year_id),
                                             mean, na.rm = TRUE))
  flow_compiled <- flow_compiled[!is.na(flow_compiled[, 1]), ]
  
  # replace missing temperature data:
  #   - King = Ovens (2011, 2012, 2017, 2018, 2019 filled with mean of adjacent years)
  #   - Murray 1999-2002 filled with average of Murray 2003-2008
  if ("spwntmp_c" %in% vars_to_include) {
    flow_compiled[system == 4 & year %in% c(1:4), "spwntmp_c"] <- mean(flow_compiled[system == 4 & year %in% c(5:10), "spwntmp_c"])
    flow_compiled[system == 3, "spwntmp_c"] <- flow_compiled[match(paste0("5", year[system == 3]), paste0(system, year)), "spwntmp_c"]
    flow_compiled[system == 3 & year %in% c(12:13), "spwntmp_c"] <- mean(flow_compiled[system == 3 & year %in% c(11, 14), "spwntmp_c"])
    flow_compiled[system == 3 & year %in% c(19), "spwntmp_c"] <- mean(flow_compiled[system == 5 & year %in% c(18, 20), "spwntmp_c"])
  }
  
  # standardised flow data
  flow_std <- apply(flow_compiled, 2, scale)
  flow_scales <- apply(flow_compiled, 2, extract_standards)
  
  # fit model
  mod[[i]] <- fit_ccr(response = response_matrix,
                 length_age_matrix = length_age_matrix,
                 predictors = flow_std,
                 effort = effort,
                 system = system, year = year,
                 mcmc_settings = list(n_samples = 10000, warmup = 5000))
  
  # validate model
  mod_cv[[i]] <- validate(mod[[i]], folds = 10, year = TRUE, cohort = TRUE)
  
}

saveRDS(mod, file = "outputs/fitted-models.rds")
saveRDS(mod_cv, file = "outputs/validated-models.rds")
