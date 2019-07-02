# set working directory
setwd("~/Dropbox/research/catch-curves/")

# load packages
library(rstanarm)

# load some helper functions
source("code/validate_glmer.R")

# load fitted models and data
all_mods <- dir("outputs/fitted")
all_mods <- all_mods[grep("-model.rds", all_mods)]
mod_list <- list()
for (i in seq_along(all_mods))
  mod_list[[i]] <- readRDS(paste0("outputs/fitted/", all_mods[i]))
data_set <- readRDS("outputs/fitted/data-set.rds")

# specify random effects
formula_tmp <- c(~ (rrang_vec + rrang_ym1_vec +
                      psprw_vec + psprw_ym1_vec +
                      psumw_vec + psumw_ym1_vec + 
                      minwin_vec + spwntmp_vec | system_vec) +
                   (-1 +
                      rrang_vec + rrang_ym1_vec +
                      psprw_vec + psprw_ym1_vec +
                      psumw_vec + psumw_ym1_vec +
                      minwin_vec + spwntmp_vec | age_factor),
                 ~ (rrang_vec + rrang_ym1_vec +
                      psprw_vec + psprw_ym1_vec +
                      psumw_vec + psumw_ym1_vec + 
                      minwin_vec + spwntmp_vec | system_vec),
                 ~ (1 | system_vec),
                 ~ (-1 + rrang_vec + rrang_ym1_vec +
                      psprw_vec + psprw_ym1_vec +
                      psumw_vec + psumw_ym1_vec +  
                      minwin_vec + spwntmp_vec | age_factor),
                 ~ 0)

# validate models
mod_cv <- list()
for (i in seq_along(mod_list)) {
  mod_cv[[i]] <- validate_glmer(mod_list[[i]], folds = 10,
                                settings = list(iter = 5000, re.form = formula_tmp[[i]]))
}
names(mod_cv) <- gsub("-", "_", sapply(strsplit(all_mods, "\\."), function(x) x[1]))
saveRDS(mod_cv, file = "outputs/fitted/mod_cv.rds")

# what about with site-based folds?
fold_list <- lapply(seq_along(unique(data_set$system_vec)), function(i) which(data_set$system_vec == i))
mod_transferability <- list()
for (i in seq_along(mod_list)) {
  mod_transferability[[i]] <- validate_glmer(mod_list[[i]], folds = fold_list,
                                             settings = list(iter = 5000, re.form = formula_tmp[[i]]))
}
names(mod_transferability) <- gsub("-", "_", sapply(strsplit(all_mods, "\\."), function(x) x[1]))
saveRDS(mod_transferability, file = "outputs/fitted/mod_transferability.rds")
