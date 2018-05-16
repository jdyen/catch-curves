# work out length-age relationship

oti_data <- vector("list", length = 4)
oti_data[[1]] <- read.csv("./data/MC.csv")
oti_data[[2]] <- read.csv("./data/TC.csv")
oti_data[[3]] <- read.csv("./data/YB.csv")
oti_data[[4]] <- read.csv("./data/SP.csv")
names(oti_data) <- c("murraycod", "troutcod", "goldenperch", "silverperch")

filter_fun <- function(x) {
  
  x_id <- unique(x$ID)
  x_len <- rep(NA, length(x_id))
  x_age <- rep(NA, length(x_id))
  
  for (i in seq_along(x_id)) {
    
    x_sub <- x[which(x$ID == x_id[i]), ]
    x_len[i] <- max(x_sub$Length)
    x_age[i] <- max(x_sub$Age)
    
  }
  
  x_filtered <- data.frame(id = x_id,
                           length = x_len,
                           age = x_age)
  
  x_filtered
  
}

coef_extract <- function(x_filtered) {
  
  params <- coef(lm(x_filtered$age ~ x_filtered$length))
  names(params) <- c("intercept", "slope")
  
  params
  
}

length_to_age <- function(x, params) {
  
  out <- round(params[1] + params[2] * x)
  out <- ifelse(out <= 0, NA, out)
  
  out

}

params_all <- t(sapply(oti_data, function(x) coef_extract(filter_fun(x))))
