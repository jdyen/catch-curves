# R code to load flow data for each river reach
setwd("~/Dropbox/research/catch-curves/")

# need some packages
library(lubridate)
library(dplyr)

# load some helper functions
source("code/length_to_mass_calculations.R")
source("code/length_to_mass_calculations_sra.R")
source("code/data-helpers.R")
source("code/helpers.R")

# load the data
source("code/load-survey-data.R")

# save loaded data
saveRDS(alldat, file = "data/data-loaded-Jun19.rds")

# load predictor (flow) files into a list from all files with "to_use" suffix
pred_list <- dir("data/flow-data/")[grep("to_use", dir("data/flow-data/"))]
predictors <- vector("list", length = length(pred_list))

# make sure dates parse correctly
date_parse <- c(rep("dmy_HM", 3), "dmy", rep("dmy_HM", 4), "dmY", rep("dmy_HM", 9))

for (i in seq_along(pred_list)) {
  predictors[[i]] <- read.csv(paste0("./data/flow-data/", pred_list[i]),
                              stringsAsFactors = FALSE)
  if (i == 9) {
    predictors[[i]]$day <- day(parse_date_time(predictors[[i]]$date,
                                               orders = c("dmy_HM")))
    predictors[[i]]$day[c(21389:nrow(predictors[[i]]))] <-
      as.numeric(predictors[[i]]$date[c(21389:nrow(predictors[[i]]))])
    predictors[[i]]$date <- paste(predictors[[i]]$day,
                                  predictors[[i]]$month,
                                  predictors[[i]]$year,
                                  sep = "/")
  }
  predictors[[i]]$Event_Date <- parse_date_time2(predictors[[i]]$date,
                                                 orders = date_parse[i],
                                                 cutoff_2000 = 18)
  predictors[[i]] <- predictors[[i]][year(predictors[[i]]$Event_Date) > 1999, ]
  predictors[[i]] <- predictors[[i]][!is.na(predictors[[i]]$date), ]
}
names(predictors) <- sapply(strsplit(pred_list, "_"), function(x) paste(x[1], x[2], sep = "_"))

flow_std <- vector("list", length = length(predictors))
for (i in seq_along(predictors)) {

  x <- predictors[[i]]
  
  flow_data <- x[, grep("discharge_ml", colnames(x))]
  if (is.matrix(flow_data) | is.data.frame(flow_data)) {
    if (ncol(flow_data) > 1)
      flow_data <- flow_data[, -grep("qc", colnames(flow_data))]
  }
  
  lt_med_all <- median(flow_data, na.rm = TRUE)
  lt_medwin <- median(flow_data[month(x$Event_Date) %in% c(6:8)], na.rm = TRUE)

  flow_std[[i]] <- data.frame(date = x$Event_Date,
                              day = day(x$Event_Date),
                              month = month(x$Event_Date),
                              year = year(x$Event_Date),
                              flow_actual = flow_data,
                              flow_std_winter = flow_data / lt_medwin,
                              flow_std_all = flow_data / lt_med_all,
                              long_term_winter_median = rep(lt_medwin, length(flow_data)),
                              long_term_median = rep(lt_med_all, length(flow_data)))
  
}
names(flow_std) <- names(predictors)

for (i in seq_along(flow_std))
  write.csv(flow_std[[i]], paste0("data/flow-std/", names(flow_std)[i], ".csv"))



