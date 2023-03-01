library(tidyverse)

# import data
Charlie_JAN26 <- readxl::read_excel("~/GitHub/cd_workflow_2023/Data/Charlie_KPS_26JAN_22.xlsx") %>% 
  janitor::clean_names()
Thailand_15Feb_RAW <- readxl::read_excel("Data/Thailand_15Feb_RAW.xlsx") %>% 
  janitor::clean_names()
KPS_15Feb_RAW <- readxl::read_excel("Data/SamplingCampaign_15Feb_UpdNames_RAW.xlsx") %>% janitor::clean_names()


Base1 <- Charlie_JAN26
Base2 <- Thailand_15Feb_RAW
Base3 <- KPS_15Feb_RAW

# filter to remove samples with no annotation, or no MS2 data.
Base1_NoNA <- with(Base1, Base1[!(name == "" | is.na(name)), ])
Base2_NoNA <- with(Base2, Base2[!(name == "" | is.na(name)), ])
Base3_NoNA <- with(Base3, Base3[!(name == "" | is.na(name)), ])
Base1_MS2 <- Base1_NoNA[!grepl('No MS2', Base1_NoNA$ms2),]
Base2_MS2 <- Base2_NoNA[!grepl('No MS2', Base2_NoNA$ms2),]
Base3_MS2 <- Base3_NoNA[!grepl('No MS2', Base3_NoNA$ms2),]

# drop unnecessary columns entirely, unless you have used them in the CD software.
Base1_MS2$tags = NULL
Base2_MS2$tags = NULL
Base3_MS2$tags = NULL
Base1_MS2$checked = NULL
Base2_MS2$checked = NULL
Base3_MS2$checked = NULL

# concatenate compound names and retention times to make a unique identifier, this will make things easier later on.
BaseUniqueID1 <- add_column(Base1_MS2, unique_id = NA, .after = 0)
BaseUniqueID2 <- add_column(Base2_MS2, unique_id = NA, .after = 0)
BaseUniqueID3 <- add_column(Base3_MS2, unique_id = NA, .after = 0)
BaseUniqueID1$unique_id <- str_c(BaseUniqueID1$name, "_", BaseUniqueID1$rt_min)
BaseUniqueID2$unique_id <- str_c(BaseUniqueID2$name, "_", BaseUniqueID2$rt_min)
BaseUniqueID3$unique_id <- str_c(BaseUniqueID3$name, "_", BaseUniqueID3$rt_min)

# add peak number in as another unique identifier.
BaseUniqueID_PeakNumber1 <- add_column(BaseUniqueID1, peak_number = NA, .after = 0)
BaseUniqueID_PeakNumber2 <- add_column(BaseUniqueID2, peak_number = NA, .after = 0)
BaseUniqueID_PeakNumber3 <- add_column(BaseUniqueID3, peak_number = NA, .after = 0)
BaseUniqueID_PeakNumber1$peak_number <- seq.int(nrow(BaseUniqueID_PeakNumber1))
BaseUniqueID_PeakNumber2$peak_number <- seq.int(nrow(BaseUniqueID_PeakNumber2))
BaseUniqueID_PeakNumber3$peak_number <- seq.int(nrow(BaseUniqueID_PeakNumber3))

# starting with the group area measurements, lengthen the table to erase white space.
colnames(BaseUniqueID_PeakNumber1) <- sub("*_raw_f\\d\\d*", "", colnames(BaseUniqueID_PeakNumber1))
colnames(BaseUniqueID_PeakNumber2) <- sub("*_raw_f\\d\\d*", "", colnames(BaseUniqueID_PeakNumber2))
colnames(BaseUniqueID_PeakNumber3) <- sub("*_raw_f\\d\\d*", "", colnames(BaseUniqueID_PeakNumber3))

# pivot longer all in one
Longer1 <- BaseUniqueID_PeakNumber1 %>% 
  pivot_longer(cols = group_area_1_feedpump_a:peak_rating_qc3,
               names_to = "sample",
               values_to = "result")
Longer2 <- BaseUniqueID_PeakNumber2 %>% 
  pivot_longer(cols = group_area_crh_fst1:peak_rating_thailand_sample_8,
               names_to = "sample",
               values_to = "result")
Longer3 <- BaseUniqueID_PeakNumber3 %>% 
  pivot_longer(cols = group_area_qc15x:peak_rating_solution_x2,
               names_to = "sample",
               values_to = "result")

# create a new column for sample name
SampleNames1 <- add_column(Longer1, measurement = NA)
SampleNames2 <- add_column(Longer2, measurement = NA)
SampleNames3 <- add_column(Longer3, measurement = NA)

# fill in sample names
SampleNames1 <- mutate(SampleNames1,
                      measurement = case_when(str_detect(sample, "group_area") ~ "group_area",
                                              str_detect(sample, "peak_rating") ~ "peak_rating"))
SampleNames2 <- mutate(SampleNames2,
                      measurement = case_when(str_detect(sample, "group_area") ~ "group_area",
                                              str_detect(sample, "peak_rating") ~ "peak_rating"))
SampleNames3 <- mutate(SampleNames3,
                      measurement = case_when(str_detect(sample, "group_area") ~ "group_area",
                                              str_detect(sample, "peak_rating") ~ "peak_rating"))

# clean up sample column
SampleNames1$sample <- str_replace_all(SampleNames1$sample, "group_area_", "")
SampleNames2$sample <- str_replace_all(SampleNames2$sample, "group_area_", "")
SampleNames3$sample <- str_replace_all(SampleNames3$sample, "group_area_", "")
SampleNames1$sample <- str_replace_all(SampleNames1$sample, "peak_rating_", "")
SampleNames2$sample <- str_replace_all(SampleNames2$sample, "peak_rating_", "")
SampleNames3$sample <- str_replace_all(SampleNames3$sample, "peak_rating_", "")

# pivot wider
Wider1 <- SampleNames1 %>%
  pivot_wider(names_from = measurement, values_from = result)
Wider2 <- SampleNames2 %>%
  pivot_wider(names_from = measurement, values_from = result)
Wider3 <- SampleNames3 %>%
  pivot_wider(names_from = measurement, values_from = result)

# remove NAs and filter so that the peak_rating column only has values above 5.
NoNAs1 <- drop_na(Wider1, group_area)
NoNAs2 <- drop_na(Wider2, group_area)
NoNAs3 <- drop_na(Wider3, group_area)
PeakRatingFiltered1 <- subset(NoNAs1, peak_rating > 5)
PeakRatingFiltered2 <- subset(NoNAs2, peak_rating > 5)
PeakRatingFiltered3 <- subset(NoNAs3, peak_rating > 5)

# CHANGE THIS NUMBER DEPENDING ON INTENSITY FILTER
GroupAreaFiltered1 <- subset(PeakRatingFiltered1, group_area > 100000)
GroupAreaFiltered2 <- subset(PeakRatingFiltered2, group_area > 100000)
GroupAreaFiltered3 <- subset(PeakRatingFiltered3, group_area > 100000)

# create column for replicate IDs if possible
FilteredReplicate1 <- add_column(GroupAreaFiltered1, replicate = NA)
# fill replicate file based on end of string in sample column
FilteredReplicate1 <- mutate(FilteredReplicate1,
                            replicate = case_when(
                              str_ends(sample, "a") ~ "1",
                              str_ends(sample, "b") ~ "2",
                              str_ends(sample, "c") ~ "3"))
# Add location column for each sample and remove numbers (DO NOT PUT NUMBERS IN LOCATION TITLES! e.g. if you're talking pipe_1/pipe_2, call them pipe_a/pipe_b)
FilteredReplicate1$sample_location = FilteredReplicate1$sample
FilteredReplicate1$sample_location <- stringi::stri_replace_all_regex(FilteredReplicate1$sample_location, "^\\d|\\d|_*", "")
FilteredReplicate1$sample_location <- gsub('.{1}$', '', FilteredReplicate1$sample_location)

# correct the digester numbers (or correct anything that has number separation)
FilteredDigesterCorrect1 <- add_column(FilteredReplicate1, digester_number = NA)
FilteredDigesterCorrect1 <- mutate(FilteredDigesterCorrect1,
                                  digester_number = case_when(
                                    str_detect(sample, "digester1") ~ "A",
                                    str_detect(sample, "digester2") ~ "B",
                                    str_detect(sample, "digester3") ~ "C",
                                    str_detect(sample, "digester4") ~ "D",
                                    !str_detect(sample, "digester") ~ ""))
FilteredMerge1 <- add_column(FilteredDigesterCorrect1, location = NA)
FilteredDigesterMerge1 <- FilteredMerge1 %>%
  unite("location", sample_location:digester_number)
# remove underscores
FilteredDigesterMerge1$location <- stringi::stri_replace_all_regex(FilteredDigesterMerge1$location, "_", "")
# change names back to original!
colnames(FilteredDigesterMerge1)[27] = "sample_location"
FilteredReplicate1 <- FilteredDigesterMerge1

# Remove "solo" results.
SoloRemoved1 <- plyr::ddply(FilteredReplicate1, c("unique_id", "sample_location"),
                           function(d) {if (nrow(d) > 1) d else NULL})

# rename if no filtering could be done.
SoloRemoved2 <- GroupAreaFiltered2
SoloRemoved3 <- GroupAreaFiltered3

# Split by the mass_list_search column, and make two tables for mzcloud results and mass_list results
Split1 <- split(SoloRemoved1, SoloRemoved1$annot_source_mass_list_search)
MZCloud1 <- Split1$"No results"
MassList1 <- Split1$"Full match"
Split2 <- split(SoloRemoved2, SoloRemoved2$annot_source_mass_list_search)
MZCloud2 <- Split2$"No results"
MassList2 <- Split2$"Full match"
Split3 <- split(SoloRemoved3, SoloRemoved3$annot_source_mass_list_search)
MZCloud3 <- Split3$"No results"
MassList3 <- Split3$"Full match"

# Bring together the mass lists so we can split by specific mass list.
MassListLonger1 <- MassList1 %>% 
  pivot_longer(cols = c(starts_with("mass_list_match")) ,
               names_to = "mass_list_name",
               names_prefix = "mass_list_match_",
               values_to = "mass_list_match")
MassListLonger2 <- MassList2 %>% 
  pivot_longer(cols = c(starts_with("mass_list_match")) ,
               names_to = "mass_list_name",
               names_prefix = "mass_list_match_",
               values_to = "mass_list_match")
MassListLonger3 <- MassList3 %>% 
  pivot_longer(cols = c(starts_with("mass_list_match")) ,
               names_to = "mass_list_name",
               names_prefix = "mass_list_match_",
               values_to = "mass_list_match")

# specific for Charlie dataset.
MixRemoved1 <- MassListLonger1[!grepl('mix', MassListLonger1$sample_location),]
MixRemoved2 <- MixRemoved1[!grepl('control', MixRemoved1$sample_location),]
MassListLonger1 <- MixRemoved2

# Filter for no matches, or invalid mass.
FilteredMassList1 <- MassListLonger1[!grepl('No matches found', MassListLonger1$mass_list_match),]
FilteredMassList2 <- MassListLonger2[!grepl('No matches found', MassListLonger2$mass_list_match),]
FilteredMassList3 <- MassListLonger3[!grepl('No matches found', MassListLonger3$mass_list_match),]
FilteredMZCloud1 <- MZCloud1[!grepl('Invalid Mass', MZCloud1$annot_source_mz_cloud_search),]
FilteredMZCloud2 <- MZCloud2[!grepl('Invalid Mass', MZCloud2$annot_source_mz_cloud_search),]
FilteredMZCloud3 <- MZCloud3[!grepl('Invalid Mass', MZCloud3$annot_source_mz_cloud_search),]

# split further into mass lists
SplitMassList1 <- split(FilteredMassList1, FilteredMassList1$mass_list_name)
SplitMassList2 <- split(FilteredMassList2, FilteredMassList2$mass_list_name)
SplitMassList3 <- split(FilteredMassList3, FilteredMassList3$mass_list_name)

ITN1 <- SplitMassList1$"itn_kps"
Cannabinoids1 <- SplitMassList1$"kps_cannabinoids"
ITNMetabolites1 <- SplitMassList1$"itn_cyp_metabolites"
Psychoactive1 <- SplitMassList1$"kps_psychoactive_substances_v2"
Pharmaceuticals1 <- SplitMassList1$"kps_pharmaceuticals_oct22"
NPL2 <- SplitMassList2$"kps_npl"
Psychoactive2 <- SplitMassList2$"kps_psychoactive_substances_v2"
Pharmaceuticals2 <- SplitMassList2$"kps_pharmaceuticals_oct22"
ITN3 <- SplitMassList3$"itn_kps"
ITNMetabolites3 <- SplitMassList3$"itn_cyp_metabolites"
Psychoactive3 <- SplitMassList3$"kps_psychoactive_substances_v2"
Pharmaceuticals3 <- SplitMassList3$"kps_pharmaceuticals_oct22"

# method to count unique compounds in each
# CHANGE NAME OF MASS LIST EACH TIME, NUMBER WILL PRINT IN CONSOLE.
length(unique(nnn$unique_id))
### go into a table or a flowchart?

# create a csv of filtered results.
write.csv(ITN1, "Results/ITN_Charlie_26JAN.csv", row.names = FALSE)
write.csv(Cannabinoids1, "Results/Cannabinoids_Charlie_26JAN.csv", row.names = FALSE)
write.csv(ITNMetabolites1, "Results/ITNMetabolites_Charlie_26JAN.csv", row.names = FALSE)
write.csv(Psychoactive1, "Results/Psychoactive_Charlie_26JAN.csv", row.names = FALSE)
write.csv(Pharmaceuticals1, "Results/Pharmaceuticals_Charlie_26JAN.csv", row.names = FALSE)
write.csv(FilteredMZCloud1, "Results/MZCloud_Filtered_Charlie_26JAN.csv", row.names = FALSE)
write.csv(NPL2, "Results/NPL_Thailand_15FEB.csv", row.names = FALSE)
write.csv(Psychoactive2, "Results/Psychoactive_Thailand_15FEB.csv", row.names = FALSE)
write.csv(Pharmaceuticals2, "Results/Pharmaceuticals_Thailand_15FEB.csv", row.names = FALSE)
write.csv(FilteredMZCloud2, "Results/MZCloud_Filtered_Thailand_15FEB.csv", row.names = FALSE)
write.csv(ITN3, "Results/ITN_Kgato_15FEB.csv", row.names = FALSE)
write.csv(ITNMetabolites3, "Results/ITNMetabolites_Kgato_15FEB.csv", row.names = FALSE)
write.csv(Psychoactive3, "Results/Psychoactive_Kgato_15FEB.csv", row.names = FALSE)
write.csv(Pharmaceuticals3, "Results/Pharmaceuticals_Kgato_15FEB.csv", row.names = FALSE)
write.csv(FilteredMZCloud3, "Results/MZCloud_Filtered_Kgato_15FEB.csv", row.names = FALSE)

# heatmap
#1--------
ITN1 %>% 
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
ggsave("Figures/ITN_Charlie_26JAN.pdf", width = 15, height = 5)

Cannabinoids1 %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", midpoint = 5e+07) +
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
ggsave("Figures/Cannabinoids_Charlie_26JAN.pdf", width = 15, height = 5)

ITNMetabolites1 %>% 
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
ggsave("Figures/ITNMetabolites_Charlie_26JAN.pdf", width = 15, height = 5)

Psychoactive1 %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", midpoint = 4e+09) +
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
ggsave("Figures/Psychoactive_Charlie_26JAN.pdf", width = 15, height = 6)

Pharmaceuticals1 %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", midpoint = 4e+09) +
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
ggsave("Figures/Psychoactive_Charlie_26JAN.pdf", width = 15, height = 40)

#2----------
NPL2 %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", midpoint = 8e+07) +
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
ggsave("Figures/NPL_Thailand_15FEB.pdf", width = 30, height = 70, limitsize = FALSE)

Psychoactive2 %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", midpoint = 1.6e+07) +
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
ggsave("Figures/Psychoactive_Thailand_15FEB.pdf", width = 20, height = 20)

Pharmaceuticals2 %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", midpoint = 1.5e+08) +
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
ggsave("Figures/Pharmaceuticals_Thailand_15FEB.pdf", width = 20, height = 20)

#3--------
ITN3 %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", midpoint = 1.2e+08) +
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
ggsave("Figures/ITN_Kgato_15FEB.pdf", width = 30, height = 20)

ITNMetabolites3 %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", midpoint = 1.2e+07) +
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
ggsave("Figures/ITNMetabolites_Kgato_15FEB.pdf", width = 20, height = 25)

Psychoactive3 %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", midpoint = 8e+07) +
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
ggsave("Figures/Psychoactive_Kgato_15FEB.pdf", width = 20, height = 25)

Pharmaceuticals3 %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", high = "orange", mid = "yellow", midpoint = 6e+08) +
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
ggsave("Figures/Pharmaceuticals_Kgato_15FEB.pdf", width = 20, height = 125, limitsize = FALSE)

