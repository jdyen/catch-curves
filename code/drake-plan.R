# load packages
library(drake)
library(greta)

# source functions
source("./code/helpers.R")
source("./code/load-data.R")

# create workflow
plan <- drake_plan(
  vefmap_data = get_data(file_in("./data/VEFMAP_FISH_20171024.csv")),
  snags_data = get_data(file_in("./data/SNAGS_FISH_20171205.csv")),
  ovens_data = get_data(file_in("./data/vba_ovens_2008_2017.csv")),
  all_data = compile_data(vefmap_data, snags_data, ovens_data),
  northern_data = subset_systems(all_data),
  clean_species_data = clean_spp(northern_data),
  length_conversions = calculate_length_conversions(clean_species_data),
  filled_data = impute_weights(clean_species_data, length_conversions),
  clean_reach_data = clean_reaches(filled_data),
  flow_data = calculate_flow_stats(dir("./data/flow-data/")[grep("to_use", dir("./data/flow-data"))]),
  data_with_flow = add_flow_data(clean_reach_data, flow_data),
  subsetted_data = subset_spp(data_with_flow, sp_sub <- c("goldenperch", "murraycod",
                                                          "troutcod", "silverperch")),
  otolith_data = load_otolith_data(file_in(c("./data/MC.csv", "./data/TC.csv",
                                             "./data/YB.csv", "./data/SP.csv"))),
  final_data = add_age_data(subsetted_data, otolith_data),
  catch_curves = create_catch_curves(final_data),
  strings_in_dots = "literals"
)

config <- drake_config(plan)
vis_drake_graph(config, targets_only = TRUE)

make(plan)
