# work out length-age relationship

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

