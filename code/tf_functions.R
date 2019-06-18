# set up tf functions to calculate covariances quickly
tf <- tensorflow::tf
op <- greta::.internals$nodes$constructors$op

#' @export
hist.greta_array <- function(x, ...) {
  
  settings <- list(...)
  
  breaks <- settings$breaks
  
  # does x have appropriate dimensions?
  dims <- dim(x)
  if (length(dims) != 2 | dims[2] != 1) 
    stop("x must be a n x 1 column vector")
  
  # are breaks provided?
  if (is.null(breaks))
    stop("breaks must be provided", call. = FALSE)
  
  nbreak <- as.integer(length(breaks))
  nbins <- nbreak - 1L
  
  # we need the lower and upper bounds
  breaks_range <- c(min(breaks), max(breaks))
  
  # are the breaks appropriate integers?
  raw_breaks <- seq(min(breaks), max(breaks), by = 1)
  if (nbreak != length(raw_breaks))
    stop("breaks must be sequential integers", call. = FALSE)
  
  op("hist",
     x,
     operation_args = list(nbins = nbins, breaks_range = breaks_range),
     tf_operation = "tf_hist",
     dim = c(nbreak - 1L, 1))
  
}

tf_hist <- function(x, nbins, breaks_range) {
  
  tf$histogram_fixed_width(x,
                           value_range = tf$constant(breaks_range, dtype = tf$float64),
                           nbins = nbins,
                           dtype = tf$int64)
  
}
