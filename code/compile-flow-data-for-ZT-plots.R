# calc flow metrics all systems all years
# prepare all data
setwd("~/Dropbox/research/catch-curves/")

# need some packages
library(lubridate)
library(dplyr)

# load some helper functions
source("code/data-helpers.R")
source("code/helpers.R")

# R code to load flow data for each river reach
# load compiled survey data
alldat <- readRDS("data/data-loaded-Aug19.rds")
year_string <- as.character(1999:2019)
date_string <- ymd(paste0(year_string, "-05-23"))
sys_reach <- unique(paste(alldat$system, alldat$reach, sep = "_"))
data_matrix <- expand.grid(date_formatted = date_string,
                           system = sys_reach,
                           stringsAsFactors = FALSE)
sys_reach_tmp <- strsplit(data_matrix$system, "_")
data_matrix$system <- sapply(sys_reach_tmp, function(x) x[1])
data_matrix$reach <- sapply(sys_reach_tmp, function(x) x[2])

# load predictor (flow) files into a list from all files with "to_use" suffix
pred_list <- dir("data/flow-data/")[grep("to_use", dir("data/flow-data/"))]
predictors <- vector("list", length = length(pred_list))

# make sure dates parse correctly
date_parse <- c(rep("dmy_HM", 2), "dmy", rep("dmy_HM", 4), "dmY", rep("dmy_HM", 9))

for (i in seq_along(pred_list)) {
  predictors[[i]] <- read.csv(paste0("./data/flow-data/", pred_list[i]),
                              stringsAsFactors = FALSE)
  if (i == 8) {
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
  predictors[[i]] <- predictors[[i]][year(predictors[[i]]$Event_Date) > 1990, ]
  predictors[[i]] <- predictors[[i]][!is.na(predictors[[i]]$date), ]
}
names(predictors) <- sapply(strsplit(pred_list, "_"), function(x) paste(x[1], x[2], sep = "_"))

# backfill Ovens@Wang with temperature data from logger at Peechelba.
predictors$ovens_wangaratta$water_temp_c <- 
  predictors$ovens_peechelba$water_temp_c[match(predictors$ovens_wangaratta$Event_Date,
                                                predictors$ovens_peechelba$Event_Date)]

# calculate lag-0 flow metrics
out_test <- sapply(seq_len(nrow(data_matrix)),
                   calc_flow_fn,
                   data = data_matrix,
                   predictors = predictors,
                   na_thresh = 0.2)
flow_data <- t(out_test)

# repeat to get data lagged by 1 and 2 years
out_ym1 <- sapply(seq_len(nrow(data_matrix)),
                  calc_flow_fn,
                  data = data_matrix, 
                  predictors = predictors,
                  na_thresh = 0.2,
                  year_lag = 1)
out_ym2 <- sapply(seq_len(nrow(data_matrix)),
                  calc_flow_fn,
                  data = data_matrix,
                  predictors = predictors,
                  na_thresh = 0.2,
                  year_lag = 2) 
out_ym1 <- t(out_ym1)
out_ym2 <- t(out_ym2)
colnames(out_ym1) <- paste0(colnames(flow_data), "_ym1")
colnames(out_ym2) <- paste0(colnames(flow_data), "_ym2")
flow_data <- cbind(flow_data, out_ym1, out_ym2)

# easier to work with data.frame
flow_data <- as.data.frame(flow_data)

# add system and date
flow_data$system <- data_matrix$system
flow_data$reach <- data_matrix$reach
flow_data$date <- data_matrix$date_formatted

# reorder columns
flow_data <- flow_data[, c("system", "reach", "date", colnames(flow_data)[seq_len(ncol(flow_data) - 3)])]

# fix column names
colnames(flow_data) <- gsub("_win", "", colnames(flow_data))

# save outupts
write.csv(flow_data, file = "outputs/flow_data_catch_curves-v2.csv", row.names = FALSE)
