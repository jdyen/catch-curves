# helpers for catch curve analysis
remove_correlated <- function(flow) {
  
  # calculate correlations > 0.7 and remove one at a time, taking left-most of the columns
  #   correlated with the most other variables
  cor_vals <- cor(flow, use = "complete")
  sum_correlated <- apply(cor_vals, 2, function(x) sum(abs(x) > 0.7) - 1)
  while(any(sum_correlated > 0)) {
    cor_vals <- cor_vals[-max(which(sum_correlated == max(sum_correlated))), -max(which(sum_correlated == max(sum_correlated)))]
    sum_correlated <- apply(cor_vals, 2, function(x) sum(abs(x) > 0.7) - 1)
  }
  flow[, match(names(sum_correlated), colnames(flow))]
  
}

na_replace_fun <- function(x) {
  
  if (any(is.na(x))) 
    x[is.na(x)] <- mean(x, na.rm = TRUE)

  x
  
}

calc_flow_pc <- function(flow, scale = FALSE) {
  
  flow <- apply(flow, 2, na_replace_fun)
  if (scale)
    flow <- scale(flow)
  pc_out <- princomp(flow)
  
  pc_out
  
}
