# install/load packages 
library(tidyverse)
library(dplyr)

# import data
Raw_Data <- readr::read_csv("~/Desktop/Charlie_CD/Charlie_CD/Raw_Data.csv") %>% 
  janitor::clean_names()

# drop unnecessary columns entirely
Raw_Data$tags = NULL
Raw_Data$checked = NULL

Base <- Raw_Data

# find way to separate measurement from file name
# lengthen the group areas to erase white space
Group_Area <- select(Base, -c(starts_with("peak")))
Group_Area_Long <- pivot_longer(Group_Area,
                                c(starts_with("group")),
                                names_to = "sample",
                                names_prefix = "group_area_",
                                values_to = "group_area")

Peak_Rating <- select(Base, -c(starts_with("group")))
Peak_Rating_Long <- pivot_longer(Peak_Rating,
                                 c(starts_with("peak")),
                                 names_to = "sample",
                                 names_prefix = "peak_rating_",
                                 values_to = "peak_rating")
# clean up sample names column
Peak_Rating_Long$sample <- str_replace_all(Peak_Rating_Long$sample, "^*_raw_f[\\dd]*", "")

# join both tables together
Joined <- merge(Group_Area_Long, Peak_Rating_Long, by = c("name":"sample"))