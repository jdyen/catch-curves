# script to run model selection on all combinations of three or more variables 

# load packages
library(lubridate)
library(greta)
library(future.apply)

# load some helper functions
source("code/helpers.R")
source("code/fit_ccr.R")

# load compiled survey data
alldat <- readRDS("data/survey-data-loaded.rds")

# filter survey data to MC
to_keep <- alldat$scientific_name == "Maccullochella peelii"
alldat <- alldat[to_keep, ]

# load otolith data
oti_data <- readRDS("data/otolith-data-loaded.rds")

# filter to MC
oti_data <- oti_data[oti_data$SPECIES == "Maccullochella peelii", ]

# need to load flow data
flow_data <- readRDS("data/flow-data-loaded.rds")

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

# fill NAs in length data for padded rows
survey_data$length_mm[survey_data$dataset == "PADDED_TO_GET_FLOW_YEARS"] <- 0

# remove any rows with missing data
flow_data <- flow_data[apply(survey_data, 1, function(x) !any(is.na(x))), ]
survey_data <- survey_data[apply(survey_data, 1, function(x) !any(is.na(x))), ]

# remove padded rows from survey data to calculate age matrix
survey_data_filtered <- survey_data[survey_data$dataset != "PADDED_TO_GET_FLOW_YEARS", ]

# bin otolith data by age and size: offset by -0.4 (0-0.6 = YOY, 0.6-1.6 = 1YO, etc.)
oti_data$age_class <- cut(oti_data$AGE, breaks = c(-0.4:ceiling(max(oti_data$AGE, na.rm = TRUE))),
                          labels = FALSE)
size_breaks <- 10 * length_from_age(c(-0.4:140.5), length_inf = 150, time_zero = 6, k_param = 0.0011, c_param = -103)

oti_data$len_class <- cut(oti_data$T_Length..mm., breaks = size_breaks, labels = FALSE)
length_age_matrix <- classify(oti_data$len_class, oti_data$age_class)
length_age_matrix <- rbind(length_age_matrix,
                           matrix(0,
                                  nrow = (length(size_breaks) - 1 - nrow(length_age_matrix)),
                                  ncol = ncol(length_age_matrix)))
length_age_matrix <- cbind(length_age_matrix,
                           matrix(0,
                                  nrow = nrow(length_age_matrix),
                                  ncol = (length(size_breaks) - 1 - ncol(length_age_matrix))))

# how many length classes do we want to include?
n_len <- 6

# how many age classes do we want to keep?
n_age <- 6

# let's just keep those from the age_length_matrix
length_age_matrix <- length_age_matrix[seq_len(n_len), seq_len(n_age)]

# need to bin the survey data by lengths
response_matrix <- do.call(
  rbind, tapply(survey_data_filtered$length_mm,
                list(survey_data_filtered$system, survey_data_filtered$year),
                hist_fn, breaks = size_breaks))

# pull out adult abundances
response_matrix_full <- do.call(
  rbind, tapply(survey_data$length_mm,
                list(survey_data$system, survey_data$year),
                hist_fn, breaks = size_breaks))
adult_catch <- apply(response_matrix_full[, 5:ncol(response_matrix_full)], 1, sum)

# filter length class survey data to the length classes we care about (1:n_len)
response_matrix <- response_matrix[, seq_len(n_len)]

# pad diagonal of length-age matrix
diag(length_age_matrix) <- diag(length_age_matrix) + 1

# define system/year predictors
system <- c(tapply(survey_data_filtered$system_id,
                   list(survey_data_filtered$system, survey_data_filtered$year),
                   unique))
system <- system[!is.na(system)]
year <- c(tapply(survey_data_filtered$year_id,
                 list(survey_data_filtered$system, survey_data_filtered$year),
                 unique))
year <- year[!is.na(year)]

# how much effort went into each survey? This should be summed over all sits in a given system
#   (but not double-count effort for each individual fish)
sys_site <- paste(survey_data_filtered$system, survey_data_filtered$site, sep = "_")
effort_site <- tapply(survey_data_filtered$effort,
                      list(sys_site, survey_data_filtered$year),
                      unique)
year_ids <- rep(colnames(effort_site), each = nrow(effort_site))
system_ids <- sapply(strsplit(rownames(effort_site), "_"), function(x) x[1])
system_ids <- rep(system_ids, times = ncol(effort_site))
effort_site <- sapply(c(effort_site), sum)
effort <- tapply(effort_site, list(system_ids, year_ids), sum)
effort <- effort[effort > 0]
sys_site_full <- paste(survey_data$system, survey_data$site, sep = "_")
effort_full <- tapply(survey_data$effort,
                      list(sys_site_full, survey_data$year_id),
                      unique)
year_ids_full <- rep(colnames(effort_full), each = nrow(effort_full))
system_ids_full <- sapply(strsplit(rownames(effort_full), "_"), function(x) x[1])
system_ids_full <- rep(system_ids_full, times = ncol(effort_full))
effort_full <- sapply(c(effort_full), sum)
effort_full <- tapply(effort_full, list(system_ids_full, as.integer(year_ids_full)), sum)

# want to identify years with synthetic data because effort and catch is incorrect
synth_data <- survey_data[survey_data$dataset == "PADDED_TO_GET_FLOW_YEARS", ]
all_years <- tapply(synth_data$year, synth_data$system, function(x) sort(unique(x)))
all_years$murray <- sort(c(all_years$murray, 2011:2018))
observed_years <- tapply(survey_data_filtered$year, survey_data_filtered$system, function(x) sort(unique(x)))
years_of_surveys <- mapply(`%in%`, all_years, observed_years)
years_of_surveys_matrix <- t(effort_full)
years_of_surveys_matrix[years_of_surveys_matrix > 0] <- do.call(c, years_of_surveys)
years_of_surveys_matrix <- t(years_of_surveys_matrix)

# compile flow predictors
all_vars <- c("spawning_variability", "prop_spring_lt",
              "prop_max_antecedent_lt", "prop_summer_lt",
              "prop_winter_lt",
              "spawning_temp")
var_sets <- c(
  list(all_vars),
  combn(all_vars, m = 5, simplify = FALSE),
  combn(all_vars, m = 4, simplify = FALSE),
  combn(all_vars, m = 3, simplify = FALSE)
)
include_adults <- rep(TRUE, length(var_sets))

for (mod_id in seq_along(var_sets)) {
  
  vars_to_include <- var_sets[[mod_id]]
  
  # prepare flow data
  sys_year <- data.frame(system = rep(sort(unique(survey_data$system_id)), times = length(unique(survey_data$year_id))),
                         year = rep(sort(unique(survey_data$year_id)), each = length(unique(survey_data$system_id))))
  flow_data$system_coded <- match(flow_data$system, systems_to_keep)
  flow_data$year_coded <- as.numeric(flow_data$year) - 1992
  idx <- match(paste(sys_year$system, sys_year$year, sep = "_"),
               paste(flow_data$system_coded, flow_data$year_coded, sep = "_"))
  flow_compiled <- flow_data[idx[!is.na(idx)], vars_to_include]
  flow_compiled <- apply(flow_compiled, 2, unlist)
  sys_year <- sys_year[!is.na(idx), ]
  
  #   - King = Ovens (2004-200? filled with average of following 4 years)
  #   - Murray 1993-2003 filled with average of Murray 2004-2008
  if ("spawning_temp" %in% vars_to_include) {
    flow_compiled[sys_year$system == 4 & sys_year$year <= 10, "spawning_temp"] <-
      mean(flow_compiled[sys_year$system == 4 & sys_year$year %in% c(11:16), "spawning_temp"])
    flow_compiled[sys_year$system == 5 & sys_year$year == 13, "spawning_temp"] <-
      mean(flow_compiled[sys_year$system == 5 & sys_year$year %in% c(12, 14), "spawning_temp"])
    flow_compiled[sys_year$system == 3, "spawning_temp"] <- 
      flow_compiled[match(paste("5", sys_year$year[sys_year$system == 3], sep = "_"),
                          paste(sys_year$system, sys_year$year, sep = "_")), "spawning_temp"]
    flow_compiled[sys_year$system == 3 & sys_year$year == 9, "spawning_temp"] <-
      flow_compiled[sys_year$system == 3 & sys_year$year == 11, "spawning_temp"]
  }
  
  # standardised flow data
  if (!is.matrix(flow_compiled))
    flow_compiled <- matrix(flow_compiled, ncol = 1)
  flow_std <- apply(flow_compiled, 2, scale)
  flow_scales <- apply(flow_compiled, 2, extract_standards)
  
  # add total abundance as a predictor
  years_of_surveys <- c(years_of_surveys_matrix)
  years_of_surveys <- years_of_surveys[effort_full > 0]
  effort_tmp <- effort_full[effort_full > 0]
  for (i in seq_len(max(sys_year$system))) {
    idx <- sys_year$system == i & years_of_surveys == 0
    idy <- sys_year$system == i & years_of_surveys == 1
    min_year <- min(sys_year$year[idy])
    idy <- sys_year$system == i & sys_year$year %in% c(min_year:(min_year + 3))
    effort_tmp[idx] <- mean(effort_tmp[idy])
    adult_catch[idx] <- median(adult_catch[idy])
  }
  if (include_adults[mod_id]) {
    flow_std <- cbind(flow_std, "adult_cpue" = c(scale(adult_catch / effort_tmp)))
    flow_scales <- cbind(flow_scales, "adult_cpue" = c(mean(adult_catch / effort_tmp), sd(adult_catch / effort_tmp)))
  }
  
  # save standardised flow values (need for final model, for plotting)
  if (mod_id == 1)
    saveRDS(flow_scales, file = "data/flow-standardisation.rds")
  
  # set future
  plan(multisession)
  
  # fit model
  mod <- fit_ccr(
    response = response_matrix,
    length_age_matrix = length_age_matrix,
    predictors = flow_std,
    effort = effort,
    system = system,
    year = year,
    sys_year_flow = sys_year,
    include = list(sys_flow = TRUE, survey = TRUE, predictors = TRUE),
    mcmc_settings = list(n_samples = 25000, warmup = 50000, chains = 12, thin = 1),
    optim_start = FALSE
  )
  
  # thin draws post sampling (errors in greta if thinned during)
  mod$draws <- lapply(mod$draws, function(x) x[seq(1, nrow(x), by = 25), ])
  
  # save outputs
  saveRDS(mod, file = paste0("outputs/fitted/mod_variant_", format(Sys.time(), "%Y%m%d_%H%M"), ".rds"))
  
}

