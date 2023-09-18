# PLEASE REVIEW AND CHANGE ANY FUNCTIONAL CODE WRITTEN IN CAPITAL LETTERS
# LOAD TIDYVERSE AND CHECK PACKAGE UPDATES
library(tidyverse)

# IMPORT YOUR CD DATA
batch_3 <- readxl::read_excel("Leeds/08_Sept_Batch3_v1.xlsx") %>% 
  janitor::clean_names()
batch_5 <- readxl::read_excel("Leeds/08Sept_batch5_v1.xlsx") %>% 
  janitor::clean_names()

# RENAME DATAFRAME FOR CODE TO WORK WITH MINIMAL CHANGES
basedata_3 <- batch_3
basedata_5 <- batch_5

# DROP THESE COLUMNS UNLESS THEY HAVE BEEN USED IN YOUR CD WORKFLOW
basedata_3$tags = NULL
basedata_3$checked = NULL
basedata_5$tags = NULL
basedata_5$checked = NULL

# FILTER SAMPLES WITH NO COMPOUND NAMES AND NO MS2 DATA (IF NEEDED)
basedatawithcompoundnames_3 <- with(basedata_3, basedata_3[!(name == "" | 
                                                         is.na(name)), ])
ms2dataonly_3 <- basedatawithcompoundnames_3[!grepl('No MS2', 
                                                basedatawithcompoundnames_3$ms2),]
basedatawithcompoundnames_5 <- with(basedata_5, basedata_5[!(name == "" | 
                                                         is.na(name)), ])
ms2dataonly_5 <- basedatawithcompoundnames_5[!grepl('No MS2', 
                                                basedatawithcompoundnames_5$ms2),]

# CREATE A UNIQUE IDENTIFIER FOR EACH FEATURE USING CONCATENATION
# CHANGE BASEDATAWITHCOMPOUNDNAMES TO MS2DATAONLY IF YOU CHOOSE TO RUN THAT CODE LINE
uniqueid_3 <- add_column(basedatawithcompoundnames_3, unique_id = NA, .after = 0)
uniqueid_3$unique_id <- str_c(uniqueid_3$name, "_", uniqueid_3$rt_min)
uniqueid_5 <- add_column(basedatawithcompoundnames_5, unique_id = NA, .after = 0)
uniqueid_5$unique_id <- str_c(uniqueid_5$name, "_", uniqueid_5$rt_min)

# PEAK NUMBERS CAN BE USED AS ANOTHER IDENTIFIER
peaknumber_3 <- add_column(uniqueid_3, peak_number = NA, .after = 0)
peaknumber_3$peak_number <- seq.int(nrow(peaknumber_3))
peaknumber_5 <- add_column(uniqueid_5, peak_number = NA, .after = 0)
peaknumber_5$peak_number <- seq.int(nrow(peaknumber_5))

# REMOVE CD FILE NUMBERS FROM THE END OF SAMPLE NAMES
colnames(peaknumber_3) <- sub("*_raw_f.*", "", colnames(peaknumber_3))
colnames(peaknumber_5) <- sub("*_raw_f.*", "", colnames(peaknumber_5))

# LENGTHEN THE TABLE TO REMOVE WHITESPACE
longer_3 <- peaknumber_3 %>% 
  pivot_longer(cols = group_area_batch_3_1_a:peak_rating_qc_6,
               names_to = "sample",
               values_to = "result")
longer_5 <- peaknumber_5 %>% 
  pivot_longer(cols = group_area_batch_5_25_a:peak_rating_qc_6,
               names_to = "sample",
               values_to = "result")

# CREATE A SAMPLE NAME COLUMN AND FILL, SO WE CAN GROUP PEAK RATING AND GROUP AREA
samplenames_3 <- add_column(longer_3, measurement = NA)
samplenames_3 <- mutate(samplenames_3,
                      measurement = case_when(str_detect(sample, 
                                                         "group_area") ~ 
                                                "group_area",
                                              str_detect(sample, 
                                                         "peak_rating") ~ 
                                                "peak_rating"))
samplenames_5 <- add_column(longer_5, measurement = NA)
samplenames_5 <- mutate(samplenames_5,
                      measurement = case_when(str_detect(sample, 
                                                         "group_area") ~ 
                                                "group_area",
                                              str_detect(sample, 
                                                         "peak_rating") ~ 
                                                "peak_rating"))

# CLEAN SAMPLE NAME COLUMN
samplenames_3$sample <- str_replace_all(samplenames_3$sample, "group_area_", "")
samplenames_3$sample <- str_replace_all(samplenames_3$sample, "peak_rating_", "")
samplenames_5$sample <- str_replace_all(samplenames_5$sample, "group_area_", "")
samplenames_5$sample <- str_replace_all(samplenames_5$sample, "peak_rating_", "")

# WIDEN TABLE
wider_3 <- samplenames_3 %>%
  pivot_wider(names_from = 'measurement', values_from = 'result')
wider_5 <- samplenames_5 %>%
  pivot_wider(names_from = 'measurement', values_from = 'result')

# REMOVE NAS
nona_3 <- drop_na(wider_3, group_area)
nona_5 <- drop_na(wider_5, group_area)

# FILTER DEPENDING ON PEAK RATING NUMBER
peakrating_3 <- subset(nona_3, peak_rating > 5)
peakrating_5 <- subset(nona_5, peak_rating > 5)

# FILTER DEPENDING ON INTENSITY
grouparea_3 <- subset(peakrating_3, group_area > 100000)
grouparea_5 <- subset(peakrating_5, group_area > 100000)

# FOR TECHNICAL REPLICATES ----
# CREATE A NEW COLUMN
replicates_3 <- add_column(grouparea_3, replicate = NA)
replicates_5 <- add_column(grouparea_5, replicate = NA)

# MAKE SURE TECHNICAL REPLICATES ARE AT THE END OF THE SAMPLE NAME AND CHANGE A/B/C ACCORDINGLY
replicates_3 <- mutate(replicates_3,
                     replicate = case_when(
                       str_ends(sample, "a") ~ "1",
                       str_ends(sample, "b") ~ "2",
                       str_ends(sample, "c") ~ "3"))
replicates_5 <- mutate(replicates_5,
                     replicate = case_when(
                       str_ends(sample, "a") ~ "1",
                       str_ends(sample, "b") ~ "2",
                       str_ends(sample, "c") ~ "3"))

# REMOVE REPLICATES FROM SAMPLE NAMES
replicates_3$sample <- str_sub(replicates_3$sample, end = -3)
replicates_5$sample <- str_sub(replicates_5$sample, end = -3)

# REMOVE PEAKS WITH RESULTS IN ONLY ONE REPLICATE
soloremoved_3 <- plyr::ddply(replicates_3, c("unique_id", "sample"),
                           function(d) {if (nrow(d) > 1) d else NULL})
soloremoved_5 <- plyr::ddply(replicates_5, c("unique_id", "sample"),
                           function(d) {if (nrow(d) > 1) d else NULL})

# MAKE COLUMNS FOR SAMPLES, EG PIVOT WIDER
wider_filtered_3 <- soloremoved_3 %>%
  pivot_wider(names_from = sample, values_from = group_area)
wider_filtered_5 <- soloremoved_5 %>%
  pivot_wider(names_from = sample, values_from = group_area)

# WRITE CSV
readr::write_csv(soloremoved_3, "Leeds/Filtered_Batch_3_Results.csv")
readr::write_csv(soloremoved_5, "Leeds/Filtered_Batch_5_Results.csv")
readr::write_csv(wider_filtered_3, "Leeds/Filtered_Wide_Batch_3_Results.csv")
readr::write_csv(wider_filtered_5, "Leeds/Filtered_Wide_Batch_5_Results.csv")

#----
# SPLIT RESULTS BASED ON MASS LIST VS MZCLOUD
split_3 <- split(soloremoved_3, soloremoved_3$annot_source_mass_list_search)
mzcloud_3 <- split_3$"No results"
masslists_3 <- split_3$"Full match"
split_5 <- split(soloremoved_5, soloremoved_5$annot_source_mass_list_search)
mzcloud_5 <- split_5$"No results"
masslists_5 <- split_5$"Full match"

# MERGE MASS LISTS INTO ONE COLUMN
masslistmerged_3 <- masslists_3 %>% 
  pivot_longer(cols = c(starts_with("mass_list_match")) ,
               names_to = "mass_list_name",
               names_prefix = "mass_list_match_",
               values_to = "mass_list_match")
masslistmerged_5 <- masslists_5 %>% 
  pivot_longer(cols = c(starts_with("mass_list_match")) ,
               names_to = "mass_list_name",
               names_prefix = "mass_list_match_",
               values_to = "mass_list_match")

# FILTER FOR NO MATCHES AND INVALID MASS RESULTS
filteredmasslist_3 <- masslistmerged_3[!grepl('No matches found', masslistmerged_3$mass_list_match),]
filteredmzcloud_3 <- mzcloud_3[!grepl('Invalid Mass', mzcloud_3$annot_source_mz_cloud_search),]
filteredmasslist_5 <- masslistmerged_5[!grepl('No matches found', masslistmerged_5$mass_list_match),]
filteredmzcloud_5 <- mzcloud_5[!grepl('Invalid Mass', mzcloud_5$annot_source_mz_cloud_search),]

# SPLIT THE INDIVIDUAL MASS LISTS
splitmasslist <- split(filteredmasslist, filteredmasslist$mass_list_name)
NAMEOFMASSLIST <- splitmasslist$"NAMEOFMASSLIST"

# PRODUCE A CSV OF RESULTS
write.csv(MASSLISTNAME, "Results/MASSLISTNAME _ DATASETNAME.csv", row.names = FALSE)

# PRODUCE A HEATMAP TO QUICKLY VISUALISE RESULTS
# CHANGE HEIGHT AND WIDTH AND MIDPOINT AS NEEDED
MASSLISTNAME %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", midpoint = 2e+08) +
  labs(x = "Sample", y = "Compound Name", colour = "Intensity") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_line(colour = "gray80"),
        panel.grid.minor = element_line(colour = "gray80"),
        axis.text.x = element_text(angle = 90),
        legend.text = element_text(family = "serif", 
                                   size = 10), 
        axis.text = element_text(family = "serif", 
                                 size = 10),
        axis.title = element_text(family = "serif",
                                  size = 10, face = "bold", colour = "gray20"),
        legend.title = element_text(size = 10,
                                    family = "serif"),
        plot.background = element_rect(colour = NA,
                                       linetype = "solid"), 
        legend.key = element_rect(fill = NA)) + labs(fill = "Intensity")
ggsave("Figures/SAMEASCSV.pdf", width = 15, height = 5)
