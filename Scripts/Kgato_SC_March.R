library(tidyverse)

# import data
KPS_March_RAW <- readxl::read_excel("Data/SC_March_RAW.xlsx") %>% 
  janitor::clean_names()

# rename for the code, might be unnecessary but it's easy
Base3 <- KPS_March_RAW

# filter to remove samples with no annotation, or no MS2 data
Base3_NoNA <- with(Base3, Base3[!(name == "" | is.na(name)), ])
Base3_MS2 <- Base3_NoNA[!grepl('No MS2', Base3_NoNA$ms2),]

# drop unnecessary columns entirely, 
# unless you have used them in the CD software
Base3_MS2$tags = NULL
Base3_MS2$checked = NULL

# concatenate compound names and retention times to make a unique identifier 
# this will make things easier later on
BaseUniqueID3 <- add_column(Base3_MS2, unique_id = NA, 
                            .after = 0)
BaseUniqueID3$unique_id <- str_c(BaseUniqueID3$name, "_", 
                                 BaseUniqueID3$rt_min)

# add peak number in as another unique identifier
BaseUniqueID_PeakNumber3 <- add_column(BaseUniqueID3, peak_number = NA, 
                                       .after = 0)
BaseUniqueID_PeakNumber3$peak_number <- seq.int(nrow(BaseUniqueID_PeakNumber3))

# starting with the group area measurements 
# lengthen the table to erase white space
colnames(BaseUniqueID_PeakNumber3) <- sub("*_raw_f\\d\\d*", "", 
                                          colnames(BaseUniqueID_PeakNumber3))

# pivot longer all in one
Longer3 <- BaseUniqueID_PeakNumber3 %>% 
  pivot_longer(cols = group_area_qc15x:peak_rating_solution_x2,
               names_to = "sample",
               values_to = "result")

# create a new column for sample name
SampleNames3 <- add_column(Longer3, measurement = NA)

# fill in sample names
SampleNames3 <- mutate(SampleNames3,
                       measurement = case_when(str_detect(sample, 
                                                          "group_area") ~ 
                                                 "group_area",
                                               str_detect(sample, 
                                                          "peak_rating") ~ 
                                                 "peak_rating"))

# clean up sample column
SampleNames3$sample <- str_replace_all(SampleNames3$sample, 
                                       "group_area_", "")
SampleNames3$sample <- str_replace_all(SampleNames3$sample, 
                                       "peak_rating_", "")

# pivot wider
Wider3 <- SampleNames3 %>%
  pivot_wider(names_from = measurement, 
              values_from = result)

# remove NAs and filter so that the peak_rating column only has values above 5
NoNAs3 <- drop_na(Wider3, group_area)
PeakRatingFiltered3 <- subset(NoNAs3, 
                              peak_rating > 5)

# CHANGE THIS NUMBER DEPENDING ON INTENSITY FILTER
GroupAreaFiltered3 <- subset(PeakRatingFiltered3, 
                             group_area > 100000)

# create column for replicate IDs if possible
FilteredReplicate3 <- add_column(GroupAreaFiltered3, 
                                 replicate = NA)

# fill replicate file based on end of string in sample column
FilteredReplicate3 <- mutate(FilteredReplicate3,
                             replicate = case_when(
                               str_ends(sample, "x") ~ "1",
                               str_ends(sample, "1") ~ "1",
                               str_ends(sample, "2") ~ "2"))

# Add location column for each sample and remove numbers
# (DO NOT PUT NUMBERS IN LOCATION TITLES! e.g. if you're talking pipe_1/pipe_2, call them pipe_a/pipe_b)
FilteredReplicate3$sample_location = FilteredReplicate3$sample
FilteredReplicate3$sample_location <- stringi::stri_replace_all_regex(
  FilteredReplicate3$sample_location, "_[^_]+$", "")

SoloRemoved3 <- plyr::ddply(FilteredReplicate3, c("unique_id", 
                                                  "sample_location"),
                            function(d) {if (nrow(d) > 1) d else NULL})

# rename if no filtering could be done.
SoloRemoved2 <- GroupAreaFiltered2

# calculate sum of intensity
Summary3_compound <- SoloRemoved3 %>% 
  group_by(name, sample_location) %>% 
  summarise(sum_of_intensity = sum(group_area))

# Split by the mass_list_search column
# and make two tables for mzcloud results and mass_list results
Split3 <- split(SoloRemoved3, SoloRemoved3$annot_source_mass_list_search)
MZCloud3 <- Split3$"No results"
MassList3 <- Split3$"Full match"

# Bring together the mass lists so we can split by specific mass list.
MassListLonger3 <- MassList3 %>% 
  pivot_longer(cols = c(starts_with("mass_list_match")) ,
               names_to = "mass_list_name",
               names_prefix = "mass_list_match_",
               values_to = "mass_list_match")

# Filter for no matches, or invalid mass.
FilteredMassList3 <- MassListLonger3[!grepl('No matches found', 
                                            MassListLonger3$mass_list_match),]
FilteredMZCloud3 <- MZCloud3[!grepl('Invalid Mass', 
                                    MZCloud3$annot_source_mz_cloud_search),]

# split further into mass lists
SplitMassList3 <- split(FilteredMassList3, FilteredMassList3$mass_list_name)

ITN3 <- SplitMassList3$"itn_kps"
ITNMetabolites3 <- SplitMassList3$"itn_cyp_metabolites"
Psychoactive3 <- SplitMassList3$"kps_psychoactive_substances_v2"
Pharmaceuticals3 <- SplitMassList3$"kps_pharmaceuticals_oct22"

# more stats for once we've split into mass lists
# calculate sum of intensity
Summary3_compound <- SoloRemoved3 %>% 
  group_by(name, sample_location) %>% 
  summarise(sum_of_intensity = sum(group_area))

# method to count unique compounds in each
# CHANGE NAME OF MASS LIST EACH TIME, NUMBER WILL PRINT IN CONSOLE.
length(unique(nnn$unique_id))
### go into a table or a flowchart?

# create a csv of filtered results.
write.csv(SummaryITN3, 
          "Results/SummaryITN_Kgato_15FEB.csv", 
          row.names = FALSE)
write.csv(SummaryITNMetabolites3, 
          "Results/SummaryITNMetabolites_Kgato_15FEB.csv", 
          row.names = FALSE)
write.csv(SummaryPsychoactive3, 
          "Results/SummaryPsychoactive_Kgato_15FEB.csv", 
          row.names = FALSE)
write.csv(SummaryPharmaceuticals3, 
          "Results/SummaryPharmaceuticals_Kgato_15FEB.csv", row.names = FALSE)
write.csv(SummaryMZCloud3, 
          "Results/SummaryMZCloud_Filtered_Kgato_15FEB.csv", row.names = FALSE)
write.csv(ITN3, 
          "Results/ITN_Kgato_15FEB.csv", row.names = FALSE)
write.csv(ITNMetabolites3, 
          "Results/ITNMetabolites_Kgato_15FEB.csv", row.names = FALSE)
write.csv(Psychoactive3, 
          "Results/Psychoactive_Kgato_15FEB.csv", row.names = FALSE)
write.csv(Pharmaceuticals3, 
          "Results/Pharmaceuticals_Kgato_15FEB.csv", row.names = FALSE)
write.csv(FilteredMZCloud3, 
          "Results/MZCloud_Filtered_Kgato_15FEB.csv", row.names = FALSE)


# heatmap
#1--------
ITN1 %>% 
  filter(!is.na(name)) %>% 
  ggplot(aes(y = name, 
             x = sample, 
             fill = group_area)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", 
                       high = "orange", 
                       mid = "yellow", 
                       midpoint = 2e+08) +
  labs(x = "Sample", 
       y = "Compound Name", 
       colour = "Intensity") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_line(colour = "gray80"),
        panel.grid.minor = element_line(colour = "gray80"),
        axis.text.x = element_text(angle = 90),
        legend.text = element_text(family = "serif", 
                                   size = 10), 
        axis.text = element_text(family = "serif", 
                                 size = 10),
        axis.title = element_text(family = "serif",
                                  size = 10, 
                                  face = "bold", 
                                  colour = "gray20"),
        legend.title = element_text(size = 10,
                                    family = "serif"),
        plot.background = element_rect(colour = NA,
                                       linetype = "solid"), 
        legend.key = element_rect(fill = NA)) + labs(fill = "Intensity")
ggsave("Figures/ITN_Charlie_26JAN.pdf", width = 15, height = 5)
