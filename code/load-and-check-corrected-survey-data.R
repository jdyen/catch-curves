# load in 2019 ARI fish data and work out survey dates available

# which files do we want?
file_list <- c("VEFMAP_ALL_FISH_190731.csv",
               "SNAGS_FISH_20190130.csv",
               "OVENS_NFRC_2019.csv",
               "OVENS_2018.csv",
               "OVENS_MID_2012_2018.csv",
               "OVENS_VBA_2008_2017.csv",
               "KING_2009_2016.csv",
               "GOULBURN_NON_VEFMAP_SITES_2014_2019.csv")

# load in sites one by one
data_list <- lapply(file_list, function(file) read.csv(paste0("data/", file), stringsAsFactors = FALSE))

# which columns do we need?
cols_to_keep <- c("system", "site", "reach",
                  "event_date",
                  "gear_type", "ef_seconds_total",
                  "scientific_name",
                  "no_collected", "no_observed",
                  "total_length_mm", "weight_g",
                  "temp_c", "do_mgl", "ph", "turbidity_ntu", "secchi_depth_m")

# standardise columns in all data sets
data_list <- lapply(data_list, col_match, cols_to_keep)

# some sites are integer, some are character, standardise them
for (i in seq_along(data_list))
  data_list[[i]]$site <- as.character(data_list[[i]]$site)

# dates need to be formatted correctly and consistentl
date_types <- c("dmy_HM", "dmy_HM", "dmy", "dmy", "dmy", "dmy", "dmy", "dmy")
for (i in seq_along(data_list))
  data_list[[i]]$date_formatted <- parse_date_time(data_list[[i]]$event_date, orders = date_types[i])

# collapse into one matrix
data_matrix <- do.call(rbind, data_list)
data_sets <- c("VEFMAP", "SNAGS", "NFRC", "OVENS", "OVENS_MID", "OVENS_VBA", "KING", "GOULBURN_EXTRA")
data_matrix$dataset <- rep(data_sets, times = sapply(data_list, nrow))

# fix species names
data_matrix$spp_formatted <- tolower(data_matrix$scientific_name)
data_matrix$spp_formatted <- gsub(" ", "", data_matrix$spp_formatted)
data_matrix$spp_formatted <- gsub("\\.", "", data_matrix$spp_formatted)
data_matrix$spp_formatted <- gsub("/", "", data_matrix$spp_formatted)
data_matrix$spp_formatted <- gsub("peeliipeelii", "peelii", data_matrix$spp_formatted)

# which genera do we want?
genera_to_keep <- c("maccullochella", "cyprinus",
                    "macquaria", "melanotaenia", "bidyanus")
idx <- unlist(lapply(genera_to_keep, function(name) grep(name, data_matrix$spp_formatted)))
data_matrix <- data_matrix[idx, ]

# which species do we need to ditch?
taxa_to_lose <- c("maccullochellasp", "auratus")
idx <- unlist(lapply(taxa_to_lose, function(name) grep(name, data_matrix$spp_formatted)))
data_matrix <- data_matrix[-idx, ]

# replace scientific names with correct (consistent) names
data_matrix$scientific_name <- add_scientific_names(data_matrix$spp_formatted)
data_matrix$common_name <- add_common_names(data_matrix$spp_formatted)

# tidy up the systems column
data_matrix$system <- tolower(data_matrix$system)
data_matrix$system <- gsub(" ", "", data_matrix$system)
data_matrix$system <- gsub("river", "", data_matrix$system)
data_matrix$system <- gsub("gouburn", "goulburn", data_matrix$system)

# create a matrix of survey dates and locations
survey_matrix <- do.call(rbind, tapply(year(data_matrix$date_formatted), data_matrix$site, function(x) ifelse(1999:2019 %in% x, 1, 0)))
colnames(survey_matrix) <- c(1999:2019)
survey_matrix <- as.data.frame(survey_matrix)
survey_matrix$system <- data_matrix$system[match(rownames(survey_matrix), data_matrix$site)]
survey_matrix$site <- rownames(survey_matrix)
survey_matrix <- survey_matrix[, c("system", "site", 1999:2019)]
effort_matrix <- tapply(data_matrix$ef_seconds_total,
                        list(data_matrix$date_formatted, data_matrix$site),
                        unique)

# how many ef values do we have for each survey?
num_ef <- sapply(c(effort_matrix), length)
multiple_ef_sites <- colnames(effort_matrix)[ceiling(which(num_ef > 1) / nrow(effort_matrix))]
multiple_ef_dates <- rownames(effort_matrix)[which(num_ef > 1) %% nrow(effort_matrix)]
survey_id <- paste(multiple_ef_sites, multiple_ef_dates, sep = "_")
multiple_ef_surveys <- data_matrix[paste(data_matrix$site, data_matrix$date_formatted, sep = "_") %in% survey_id, ]
multiple_ef_surveys <- multiple_ef_surveys[order(multiple_ef_surveys$date_formatted), ]

## some surveys have EF efforts based on two or more books or banks, so need to be summed
# which sites and dates have more than one EF value?
sites_to_sum <- colnames(effort_matrix)[ceiling(which(num_ef > 1) / nrow(effort_matrix))]
dates_to_sum <- rownames(effort_matrix)[which(num_ef > 1) %% nrow(effort_matrix)]

# let's pull these out of the full data set
survey_id <- paste(sites_to_sum, dates_to_sum, sep = "_")
tmp <- data_matrix[paste(data_matrix$site, data_matrix$date_formatted, sep = "_") %in% survey_id, ]

# how many seconds does each one have?
ef_seconds <- tapply(tmp$ef_seconds_total, list(tmp$site, tmp$date_formatted), unique)

# and to which sites/dates do these correspond?
idx <- rep(rownames(ef_seconds), times = ncol(ef_seconds))
idy <- rep(colnames(ef_seconds), each = nrow(ef_seconds))

# get rid of the invalid site/date combos
ef_seconds_flat <- c(ef_seconds)
to_keep <- !sapply(ef_seconds_flat, is.null)
ef_seconds_flat <- ef_seconds_flat[to_keep]
idx <- idx[to_keep]
idy <- idy[to_keep]

# sum over all EF attempts in a given survey
ef_seconds_sum <- sapply(ef_seconds_flat, sum)

# we need to back-fill the full data set with these values
survey_id <- paste(idx, idy, sep = "_")
row_ids <- match(paste(data_matrix$site, data_matrix$date_formatted, sep = "_"), survey_id)
data_matrix$ef_seconds_total[which(!is.na(row_ids))] <- ef_seconds_sum[row_ids[!is.na(row_ids)]]

# add a reach to all sites just because (default = r1)
data_matrix$reach[is.na(data_matrix$reach)] <- 1

# replace empty observed rows with zeros
data_matrix$no_observed[is.na(data_matrix$no_observed)] <- 0

# replace empty collected rows with zeros
data_matrix$no_collected[is.na(data_matrix$no_collected)] <- 0

# convert weights to numeric
data_matrix$weight_g <- as.numeric(data_matrix$weight_g)

# calculate missing lengths from weights
unique_species <- unique(data_matrix$spp_formatted)
for (i in seq_along(unique_species)) {
  data_sub <- data_matrix[data_matrix$spp_formatted == unique_species[i], ]
  data_lm <- data_sub[data_sub$total_length_mm > 0, ]
  data_lm <- data_lm[data_lm$weight_g > 0, ]
  weight_to_length <- lm(log(data_lm$total_length_mm) ~ log(data_lm$weight_g))
  coefs <- coef(weight_to_length)
  data_sub$total_length_mm[is.na(data_sub$total_length_mm)] <- exp(coefs[1] + coefs[2] * log(data_sub$weight_g[is.na(data_sub$total_length_mm)]))
  data_matrix$total_length_mm[data_matrix$spp_formatted == unique_species[i]] <- data_sub$total_length_mm
}

# tidy workspace
rm(cols_to_keep, data_list, date_types, dates_to_sum, ef_seconds, ef_seconds_flat,
   ef_seconds_sum, effort_matrix, file_list, genera_to_keep, i, idx, idy,
   multiple_ef_dates, multiple_ef_sites, multiple_ef_surveys, num_ef,
   row_ids, sites_to_sum, survey_id, survey_matrix, taxa_to_lose, tmp, to_keep,
   data_lm, data_sub, unique_species, weight_to_length, coefs)
