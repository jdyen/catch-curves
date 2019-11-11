# R code to load flow data for each river reach

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

# repeat to get data lagged by 1 year
out_ym1 <- sapply(seq_len(nrow(data_matrix)),
                  calc_flow_fn,
                  data = data_matrix, 
                  predictors = predictors,
                  na_thresh = 0.2,
                  year_lag = 1)
out_ym1 <- t(out_ym1)
colnames(out_ym1) <- paste0(colnames(flow_data), "_ym1")
flow_data <- cbind(flow_data, out_ym1)

# easier to work with data.frame
flow_data <- as.data.frame(flow_data)

# clean up workspace
rm(i, out_test, date_parse, predictors, pred_list)
