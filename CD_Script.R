# install/load packages (working within a project so working directory is fine) 
library(tidyverse)

# import data from csv file, using the janitor package to clean the name format up.
Raw_Data <- readr::read_csv("~/Desktop/Charlie_CD/Charlie_CD/Raw_Data.csv") %>% 
  janitor::clean_names()

# drop unnecessary columns entirely, unless you have used them in the CD software.
Raw_Data$tags = NULL
Raw_Data$checked = NULL

# to make sure we don't edit any raw data files, duplicate and call the original dataset "base".
Base <- Raw_Data

# concatenate compound names and retention times to make a unique identifier, this will make things easier later on.
BaseUniqueID <- add_column(Base, unique_id = NA, .after = 0)
BaseUniqueID_Fill <- BaseUniqueID$unique_id <- str_c(BaseUniqueID$name, "_", BaseUniqueID$rt_min)

# add peak number in as another unique identifier.
BaseUniqueID_PeakNumber <- add_column(BaseUniqueID, peak_number = NA, .after = 0)
BaseUniqueID_PeakNumber$peak_number <- seq.int(nrow(BaseUniqueID_PeakNumber))

# starting with the group area measurements, lengthen the table to erase white space.
# this will also clean up names so we can see which sample each compound is found in.
Group_Area_Long <- pivot_longer(BaseUniqueID_PeakNumber, c(starts_with("group")),
                                names_to = "sample",
                                names_prefix = "group_area_",
                                values_to = "group_area")

# do the same with the peak rating data, as was done with the group area measurements.
Peak_Rating_Long <- pivot_longer(Group_Area_Long,
                                 c(contains("rating")),
                                 names_to = "sample_file",
                                 names_prefix = "peak_rating_",
                                 values_to = "peak_rating")

# clean up sample names column in the peak_rating_table, to remove the numbered file names from the CD software.
# this line uses regular expressions to tell us which part of the name to remove.
Peak_Rating_Long$sample_file <- stringi::stri_replace_all_regex(Peak_Rating_Long$sample_file, "^*_raw_f[\\dd]*", "")

# join both tables together so we have one big dataset.
# this line of code was too slow, needs making faster
# Joined <- merge(Group_Area_Long, Peak_Rating_Long, by = c("name", "sample"), all = TRUE)
# second attempt reformatting the data table.
#library(data.table)
#Group_Area_DT <- as.data.table(Group_Area_Long)
#Peak_Rating_DT <- as.data.table(Peak_Rating_Long)
#merge_on <- colnames(Group_Area_Long)[colnames(Group_Area_Long) != "group_area"]
#Joined <- merge(Group_Area_DT, Peak_Rating_DT, by = "sample")
#Joined <- merge(Group_Area_Long, Peak_Rating_Long, by = c("name", "sample"), all = TRUE)
#Joined <- merge(Group_Area_Long, Peak_Rating_Long, by = merge_on, all = TRUE)
