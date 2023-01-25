#SETUP/CLEAN----------
# install/load packages (working within a project so working directory is fine) 
library(tidyverse)

# import data from csv file, using the janitor package to clean the name format up.
Raw_Data <- readr::read_csv("~/Desktop/Charlie_CD/Charlie_CD/Raw_Data.csv") %>% 
  janitor::clean_names()

# to make sure we don't edit any raw data files, duplicate and call the original dataset "base".
Base <- Raw_Data

# drop unnecessary columns entirely, unless you have used them in the CD software.
Base$tags = NULL
Base$checked = NULL

# concatenate compound names and retention times to make a unique identifier, this will make things easier later on.
BaseUniqueID <- add_column(Base, unique_id = NA, .after = 0)
BaseUniqueID$unique_id <- str_c(BaseUniqueID$name, "_", BaseUniqueID$rt_min)

# add peak number in as another unique identifier.
BaseUniqueID_PeakNumber <- add_column(BaseUniqueID, peak_number = NA, .after = 0)
BaseUniqueID_PeakNumber$peak_number <- seq.int(nrow(BaseUniqueID_PeakNumber))

# starting with the group area measurements, lengthen the table to erase white space.
colnames(BaseUniqueID_PeakNumber) <- sub("*_raw_f\\d\\d*", "", colnames(BaseUniqueID_PeakNumber))

# pivot longer all in one
Longer <- BaseUniqueID_PeakNumber %>% 
  pivot_longer(cols = group_area_1_feedpump_a:peak_rating_qc3,
               names_to = "sample",
               values_to = "result")

# create a new column for sample name
SampleNames <- add_column(Longer, measurement = NA)

# fill in sample names
SampleNames <- mutate(SampleNames,
                      measurement = case_when(str_detect(sample, "group_area") ~ "group_area",
                                            str_detect(sample, "peak_rating") ~ "peak_rating"))

# clean up sample column
SampleNames$sample <- str_replace_all(SampleNames$sample, "group_area_", "")
SampleNames$sample <- str_replace_all(SampleNames$sample, "peak_rating_", "")

# pivot wider?
Wider <- SampleNames %>%
  pivot_wider(names_from = measurement, values_from = result)

#FILTER----------

# remove NAs and filter so that the peak_rating column only has values above 5.
NoNAs <- drop_na(Wider, group_area)
PeakRatingFiltered <- subset(NoNAs, peak_rating > 5)

# CHANGE THIS NUMBER DEPENDING ON INTENSITY FILTER
GroupAreaFiltered <- subset(PeakRatingFiltered, group_area > 100000)

# create column for replicate IDs.
FilteredReplicate <- add_column(GroupAreaFiltered, replicate = NA)

# fill replicate file based on end of string in sample column
FilteredReplicate <- mutate(FilteredReplicate,
                            replicate = case_when(
                              str_ends(sample, "a") ~ "1",
                              str_ends(sample, "b") ~ "2",
                              str_ends(sample, "c") ~ "3"))

# Add location column for each sample and remove numbers (DO NOT PUT NUMBERS IN LOCATION TITLES! e.g. if you're talking pipe_1/pipe_2, call them pipe_a/pipe_b)
FilteredReplicate$sample_location = FilteredReplicate$sample
FilteredReplicate$sample_location <- stringi::stri_replace_all_regex(FilteredReplicate$sample_location, "^\\d|\\d|_*", "")
FilteredReplicate$sample_location <- gsub('.{1}$', '', FilteredReplicate$sample_location)

# SPECIFIC FOR THIS DATASET
# correct the digester numbers
FilteredDigesterCorrect <- add_column(FilteredReplicate, digester_number = NA)
FilteredDigesterCorrect <- mutate(FilteredDigesterCorrect,
                            digester_number = case_when(
                              str_detect(sample, "digester1") ~ "A",
                              str_detect(sample, "digester2") ~ "B",
                              str_detect(sample, "digester3") ~ "C",
                              str_detect(sample, "digester4") ~ "D",
                              ))

# FIX THIS
#FilteredCorrectSamples <- add_column(FilteredDigesterCorrect, sample_names = NA)
#FilteredCorrectSamples$sample_names <- str_c(FilteredDigesterCorrect$sample_location, "_", FilteredDigesterCorrect$digester_number)

# Remove "solo" results.
SoloRemoved <- plyr::ddply(FilteredReplicate, c("unique_id", "sample"),
                      function(d) {if (nrow(d) > 1) d else NULL})
