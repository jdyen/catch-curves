# filter data down to a species and river subset

# add flow data
alldat <- cbind(alldat, flow_data)

# set systems of interest
system_sub <- c('BROKEN',
                'LOWERMURRAY',
                'CAMPASPE',
                'GOULBURN',
                'LODDON',
                'OVENS')
alldat <- alldat[-which(is.na(match(alldat$SYSTEM, system_sub))), ]

# keep spp of interest
alldat$Scientific.Name <- as.character(alldat$Scientific.Name)
alldat$Scientific.Name <- gsub(' ', '', alldat$Scientific.Name)
alldat$Scientific.Name <- gsub('Maccullochellapeeliipeelii', 'Maccullochellapeelii',
                               alldat$Scientific.Name)
alldat$Scientific.Name <- tolower(alldat$Scientific.Name)
sp_to_keep <- c('maccullochellapeelii',
                'maccullochellamacquariensis',
                'macquariaambigua',
                'bidyanusbidyanus',
                'melanotaeniafluviatilis',
                'retropinnasemoni',
                'gadopsismarmoratus',
                'cyprinuscarpio')
alldat <- alldat[alldat$Scientific.Name %in% sp_to_keep, ]
common_names <- c('murraycod',
                  'troutcod',
                  'goldenperch',
                  'silverperch',
                  'murrayriverrainbowfish',
                  'australiansmelt',
                  'riverblackfish',
                  'commoncarp')
names(common_names) <- sp_to_keep
alldat$Common.Name <- common_names[alldat$Scientific.Name]

# create SRA version of weights
alldat$WEIGHT_SRA <- alldat$WEIGHT
alldat$IMPUTED <- rep(FALSE, nrow(alldat))

# fill missing weights
sp_tmp <- unique(alldat$Common.Name)
for (i in seq_along(sp_tmp)) {
  
  # subset to single species
  dat_tmp <- alldat[which(alldat$Common.Name == sp_tmp[i]), ]
  
  # check to make sure some weights are missing      
  if (any(is.na(dat_tmp$WEIGHT))) {
    
    # use species-specific equation if it exists, generic otherwise
    if (length_weight_conversion[[sp_tmp[i]]]$n) {
      coefs <- length_weight_conversion[[sp_tmp[i]]]$coef
    } else {
      coefs <- length_weight_conversion$generic
    }
    
    if (length(length_weight_conversion_sra[[sp_tmp[i]]])) {
      coefs_sra <- length_weight_conversion_sra[[sp_tmp[i]]]
    } else {
      coefs_sra <- length_weight_conversion_sra$generic
    }
    
    # subset NA observations
    na_sub <- which(is.na(dat_tmp$WEIGHT))
    
    # estimate weight from length
    dat_tmp$WEIGHT[na_sub] <- exp(coefs[1] + coefs[2] * log(dat_tmp$totallength[na_sub]))
    dat_tmp$WEIGHT_SRA[na_sub] <- 10 ** (coefs_sra['intercept'] +
                                           coefs_sra['slope'] * log((dat_tmp$totallength[na_sub] / coefs_sra['fork_scale_fac']),
                                                                    base = 10))
    dat_tmp$IMPUTED[na_sub] <- TRUE
    
    # return estimated data to full data set
    alldat[alldat$Common.Name == sp_tmp[i], ] <- dat_tmp
    
  }
}

# arrange alldat to have more consistent set of variables
alldat <- data.frame(date = alldat$Event_Date,
                     year = alldat$YEAR,
                     site = alldat$SITE_CODE,
                     system = alldat$SYSTEM,
                     reach = alldat$Reach,
                     species = alldat$Common.Name,
                     sciname = alldat$Scientific.Name,
                     weight = alldat$WEIGHT,
                     abundance = alldat$Total.Sampled,
                     intensity = (alldat$total_no_passes * alldat$seconds),
                     mannf_mld = alldat$mannf_mld,
                     msprf_mld = alldat$msprf_mld,
                     msumf_mld = alldat$msumf_mld,
                     covaf_mld = alldat$covaf_mld,
                     mspwn_mld = alldat$mspwn_mld,
                     covsp_mld = alldat$covsp_mld,
                     maxan_mld = alldat$maxan_mld,
                     minan_mld = alldat$minan_mld,
                     madpth_m = alldat$madpth_m,
                     cvdpth_m = alldat$cvdpth_m,
                     maxdpth_m = alldat$maxdpth_m,
                     mindpth_m = alldat$mindpth_m)

# remove missing data
alldat <- alldat[!is.na(alldat$abundance), ]
alldat <- alldat[!is.na(alldat$intensity), ]
alldat <- alldat[!is.na(alldat$weight), ]

# set up species names for plots
sp_names <- data.frame(full = c('Australian smelt',
                                'Common Carp',
                                'Golden perch',
                                'Murray cod',
                                'Murray river rainbowfish',
                                'River blackfish',
                                'Silver perch',
                                'Trout cod'),
                       code = as.character(levels(alldat$species)))
