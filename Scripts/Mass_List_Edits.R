# bring data from mass list csv using janitor to clean names
library(tidyverse)

ITN_ML <- read_csv("Data/Antibiotics_ITN_MSCA_ANSWER_160616_wDTXSIDs.csv") %>%
  janitor::clean_names()
ITN_Metabolites_ML <- read_csv("Data/ITNANTIBIOTIC_CYP_Metabolites.csv") %>%
  janitor::clean_names()
Psychoactive_ML <- read_csv("Data/KPS_Psychoactive_substances_v2.csv") %>%
  janitor::clean_names()
Pharmaceuticals_ML <- read_csv("Data/KPS_Pharmaceuticals.csv") %>%
  janitor::clean_names()

# change column names in psychoactive to match for a large join. 
colnames(Psychoactive_ML)[2] = "name"

# make sure each mass list has a column with which mass list each compound is from before we join together.
ITN_ML2 <- add_column(ITN_ML, mass_list = "ITN", .after = 0)
ITN_Metabolites_ML2 <- add_column(ITN_Metabolites_ML, mass_list = "ITN_Metabolites", .after = 0)
Psychoactive_ML2 <- add_column(Psychoactive_ML, mass_list = "Psychoactive", .after = 0)
Pharmaceuticals_ML2 <- add_column(Pharmaceuticals_ML, mass_list = "Pharmaceuticals", .after = 0)