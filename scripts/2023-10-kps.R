# load libraries and packages
library(tidyverse)

# import data
level_2_kgato <- readxl::read_excel("data/level_2_kgato.xlsx") %>% 
  janitor::clean_names()

# get duplicate names
# "Roll-up" values into one row per group (per "personID") 
rolled <- level_2_kgato %>% 
  # create groups by name
  group_by(name, formula) %>% 
  
  # order the rows within each group (e.g. by date)
  arrange(rt_min, .by_group = TRUE) %>% 
  
  # For each column, paste together all values within the grouped rows, separated by ";"
  summarise(
    across(everything(),                           # apply to all columns
           ~paste0(na.omit(.x), collapse = "; "))) # function is defined which combines non-NA values

# create split and mean function for now nested values
split_and_mean <- function(x) {
  sapply(strsplit(x, ";"), function(x) mean(as.numeric(x)))
}
split_and_sum <- function(x) {
  sapply(strsplit(x, ";"), function(x) sum(as.numeric(x)))
}

# calculate means for samples
rolled[23:41] <- lapply(rolled[23:41], split_and_mean)
rolled[7:11] <- lapply(rolled[7:11], split_and_mean)
rolled[13:14] <- lapply(rolled[13:14], split_and_mean)
rolled[12] <- lapply(rolled[12], split_and_sum)
rolled[42] <- lapply(rolled[42], split_and_sum)

# remove other duplicates
rolled$category <- gsub(';.*', '', rolled$category)
rolled$annot_source_predicted_compositions <- gsub(';.*', '', rolled$annot_source_predicted_compositions)
rolled$annot_source_mass_list_search <- gsub(';.*', '', rolled$annot_source_mass_list_search)
rolled$annot_source_mz_cloud_search <- gsub(';.*', '', rolled$annot_source_mz_cloud_search)
rolled$mass_list_match_antibiotics_itn_msca_answer_160616_w_dtxsi_ds <- gsub(';.*', '', rolled$mass_list_match_antibiotics_itn_msca_answer_160616_w_dtxsi_ds)
rolled$psychoactive_substances <- gsub(';.*', '', rolled$psychoactive_substances)
rolled$pharmaceuticals_oct22 <- gsub(';.*', '', rolled$pharmaceuticals_oct22)
rolled$mass_list_match_itn_kps <- gsub(';.*', '', rolled$mass_list_match_itn_kps)
rolled$mass_list_match_clews_list <- gsub(';.*', '', rolled$mass_list_match_clews_list)
rolled$mass_list_match_efs_hram_compound_database <- gsub(';.*', '', rolled$mass_list_match_efs_hram_compound_database)
rolled$ms2 <- gsub(';.*', '', rolled$ms2)
rolled$reference_ion <- gsub(';.*', '', rolled$reference_ion)

write.csv(rolled, "rolled_kgato.csv")

# convert to long format
long <- rolled %>% 
  pivot_longer(cols = c(starts_with("sc")) ,
               names_to = "sample",
               values_to = "intensity")

# drop NaNs
shortened <- long %>% drop_na(intensity)

write.csv(shortened, "long_kgato.csv")
