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

# validate models
mod_cv <- list()
for (i in seq_along(mod_list))
  mod_cv[[i]] <- validate_glmer(mod_list[[i]], folds = 10, settings = list(iter = 5000))

# what about with site-based folds?
fold_list <- lapply(seq_along(unique(data_set$system_vec)), function(i) which(data_set$system_vec == i))
mod_transferability <- list()
for (i in seq_along(mod_list))
  mod_transferability[[i]] <- validate_glmer(mod_list[[i]], folds = fold_list, settings = list(iter = 5000))
