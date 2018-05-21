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

# define some helper functions
max_fun <- function(x) {
  
  out <- NA
  
  if (any(!is.na(x)))
    out <- max(x, na.rm = TRUE)
  
  out
  
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

# do this for all ages, then just focus on 1-3 year olds for some analyses
catch_curve_fun <- function(x, sp) {
  
  x_sub <- x[which(x$species == sp), ]
  
  flow_tmp <- x_sub[, grep("mannf$", colnames(x_sub)):grep("cspwn_tm2", colnames(x_sub))]
  
  bins <- c(0, seq_len(max(x_sub$age, na.rm = TRUE) + 1)) + 0.5
  obs_unique <- paste(x_sub$system, paste0("reach", x_sub$reach), x_sub$year, sep = "_")
  
  obs_vec <- unique(obs_unique)
  
  out <- matrix(NA, nrow = length(obs_vec), ncol = (length(bins) - 1))
  flow_data <- matrix(NA, nrow = length(obs_vec), ncol = ncol(flow_tmp))
  for (i in seq_along(obs_vec)) {
    dat_sub <- x_sub[which(obs_unique == obs_vec[i]), ]
    out[i, ] <- hist(dat_sub$age, breaks = bins, plot = FALSE)$counts
    flow_data[i, ] <- apply(flow_tmp[which(obs_unique == obs_vec[i]), ], 2, mean, na.rm = TRUE)
  }
  colnames(flow_data) <- colnames(flow_tmp)
  colnames(out) <- seq_len(ncol(out))
  rownames(out) <- obs_vec
  site_split <- strsplit(obs_vec, "_")
  site_info <- data.frame(system = sapply(site_split, function(x) x[1]),
                          reach = sapply(site_split, function(x) as.numeric(substr(x[2],
                                                                                   start = nchar(x[2]),
                                                                                   stop = nchar(x[2])))),
                          year = sapply(site_split, function(x) as.numeric(x[3])))
  
  out <- list(age_dist = out,
              info = site_info,
              flow = flow_data,
              bins = bins)
  
}

# do this for all ages, then just focus on 1-3 year olds for some analyses
size_dist_fun <- function(x, sp, nbins = 5) {
  
  x_sub <- x[which(x$species == sp), ]
  
  bins <- c(0, 100, 250, 400, 600, 1000, 1500, 15000, max(x_sub$weight, na.rm = TRUE))
  obs_unique <- paste(x_sub$system, paste0("reach", x_sub$reach), x_sub$year, sep = "_")
  
  obs_vec <- unique(obs_unique)
  
  out <- matrix(NA, nrow = length(obs_vec), ncol = (length(bins) - 1))
  for (i in seq_along(obs_vec)) {
    dat_sub <- x_sub[which(obs_unique == obs_vec[i]), ]
    out[i, ] <- hist(dat_sub$weight, breaks = bins, plot = FALSE)$counts
  }
  colnames(out) <- seq_len(ncol(out))
  rownames(out) <- obs_vec
  site_split <- strsplit(obs_vec, "_")
  site_info <- data.frame(system = sapply(site_split, function(x) x[1]),
                          reach = sapply(site_split, function(x) as.numeric(substr(x[2],
                                                                                   start = nchar(x[2]),
                                                                                   stop = nchar(x[2])))),
                          year = sapply(site_split, function(x) as.numeric(x[3])))
  
  out <- list(size_dist = out,
              info = site_info,
              bins = bins)
  
}

create_greta_covar <- function(n_sp, n_age,
                               eta_sp = 1, eta_age = 4,
                               sigma_mean = 0.0, sigma_sd = 1.0) {

  if (length(sigma_mean) == 1)
    sigma_mean <- rep(sigma_mean, 2)
  if (length(sigma_sd) == 1)
    sigma_sd <- rep(sigma_sd, 2)
  
  species_corr_est <- lkj_correlation(eta = eta_sp, dim = n_sp)
  age_corr_est <- lkj_correlation(eta = eta_age, dim = n_age)
  species_sigma <- lognormal(mean = sigma_mean[1], sd = sigma_sd[1], dim = n_sp)
  age_sigma <- lognormal(mean = sigma_mean[2], sd = sigma_sd[2], dim = n_age)
  species_covar_est <- (species_sigma %*% t(species_sigma)) * species_corr_est
  age_covar_est <- age_sigma %*% t(age_sigma) * age_corr_est
  covar_est <- kronecker(species_covar_est, age_covar_est)
  
  covar_est

}

list_full_spp_names <- function() {
  
  # set up species names for plots
  sp_names <- data.frame(full = c("Australian bass",
                                  "Australian grayling",
                                  "Australian smelt",
                                  "Black bream",
                                  "Bony bream",
                                  "Brown trout",
                                  "Carp",
                                  "Carp gudgeon",
                                  "Common galaxias",
                                  "Dwarf flathead gudgeon",
                                  "Eastern gambusia",
                                  "Eel",
                                  "Estuary perch",
                                  "Flathead galaxias",
                                  "Flathead gudgeon",
                                  "Golden perch",
                                  "Long-finned eel",
                                  "Luderick",
                                  "Mountain galaxias",
                                  "Murray cod",
                                  "Murray river rainbowfish",
                                  "Obscure galaxias",
                                  "Oriental weather loach",
                                  "Pouched lamprey",
                                  "Pygmy perch",
                                  "Rainbow trout",
                                  "Redfin",
                                  "River blackfish",
                                  "River garfish",
                                  "Roach",
                                  "Sea mullet",
                                  "Short finned eel",
                                  "Short headed lamprey",
                                  "Silver perch",
                                  "Southern pygmy perch",
                                  "Tench",
                                  "Trout cod",
                                  "Tupong",
                                  "Two spined blackfish",
                                  "Unspecked hardyhead",
                                  "Variegated pygmy perch",
                                  "Western carp gudgeon",
                                  "Yarra pygmy perch",
                                  "Yellow eyed mullet"),
                         code = as.character(levels(alldat$species)))
  
  sp_names
  
}
