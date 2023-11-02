# PLEASE REVIEW AND CHANGE ANY FUNCTIONAL CODE WRITTEN IN CAPITAL LETTERS
# LOAD TIDYVERSE AND CHECK PACKAGE UPDATES
library(tidyverse)
library(fuzzyjoin)
library(stringi)
library(hrbrthemes)
library(viridis)
library(ggforce)

# IMPORT YOUR CD DATA
drying_study <- readr::read_delim("data/raw/2023-09-naburn-drying.csv", 
                                    delim = "\t", trim_ws = TRUE) %>% 
  janitor::clean_names()

# RENAME DATAFRAME FOR CODE TO WORK WITH MINIMAL CHANGES
basedata <- drying_study

# DROP THESE COLUMNS UNLESS THEY HAVE BEEN USED IN YOUR CD WORKFLOW
basedata$tags = NULL
basedata$checked = NULL

# FILTER SAMPLES WITH NO COMPOUND NAMES AND NO MS2 DATA (IF NEEDED)
basedatawithcompoundnames <- with(basedata, basedata[!(name == "" | 
                                                         is.na(name)), ])
ms2dataonly <- basedatawithcompoundnames[!grepl('No MS2', 
                                                basedatawithcompoundnames$ms2),]

# CREATE A UNIQUE IDENTIFIER FOR EACH FEATURE USING CONCATENATION
# CHANGE BASEDATAWITHCOMPOUNDNAMES TO MS2DATAONLY IF YOU CHOOSE TO RUN THAT CODE LINE
uniqueid <- add_column(basedatawithcompoundnames, unique_id = NA, .after = 0)
uniqueid$unique_id <- str_c(uniqueid$name, "_", uniqueid$rt_min)

# PEAK NUMBERS CAN BE USED AS ANOTHER IDENTIFIER
peaknumber <- add_column(uniqueid, peak_number = NA, .after = 0)
peaknumber$peak_number <- seq.int(nrow(peaknumber))

# REMOVE CD FILE NUMBERS FROM THE END OF SAMPLE NAMES
colnames(peaknumber) <- sub("*_raw_f\\d\\d*", "", colnames(peaknumber))

# LENGTHEN THE TABLE TO REMOVE WHITESPACE
longer <- peaknumber %>% 
  pivot_longer(cols = group_area_a1_a:peak_rating_qc_p,
               names_to = "sample",
               values_to = "result")

# CREATE A SAMPLE NAME COLUMN AND FILL, SO WE CAN GROUP PEAK RATING AND GROUP AREA
samplenames <- add_column(longer, measurement = NA)
samplenames <- mutate(samplenames,
                      measurement = case_when(str_detect(sample, 
                                                         "group_area") ~ 
                                                "group_area",
                                              str_detect(sample, 
                                                         "peak_rating") ~ 
                                                "peak_rating"))

# CLEAN SAMPLE NAME COLUMN
samplenames$sample <- str_replace_all(samplenames$sample, "group_area_", "")
samplenames$sample <- str_replace_all(samplenames$sample, "peak_rating_", "")

# WIDEN TABLE
wider <- samplenames %>%
  pivot_wider(names_from = measurement, values_from = result)

# REMOVE NAS
nona <- drop_na(wider, group_area)

# FILTER DEPENDING ON PEAK RATING NUMBER
peakrating <- subset(nona, peak_rating > 5)

# FILTER DEPENDING ON INTENSITY
grouparea <- subset(peakrating, group_area > 100000)

# FOR TECHNICAL REPLICATES ----
# CREATE A NEW COLUMN
replicates <- add_column(grouparea, replicate = NA)

# MAKE SURE TECHNICAL REPLICATES ARE AT THE END OF THE SAMPLE NAME AND CHANGE A/B/C ACCORDINGLY
replicates <- mutate(replicates,
                            replicate = case_when(
                              str_ends(sample, "a") ~ "a",
                              str_ends(sample, "b") ~ "b",
                              str_ends(sample, "c") ~ "c",
                              str_ends(sample, "d") ~ "d",
                              str_ends(sample, "e") ~ "e",
                              str_ends(sample, "f") ~ "f",
                              str_ends(sample, "g") ~ "g",
                              str_ends(sample, "h") ~ "h",
                              str_ends(sample, "i") ~ "i",
                              str_ends(sample, "j") ~ "j",
                              str_ends(sample, "k") ~ "k",
                              str_ends(sample, "l") ~ "l",
                              str_ends(sample, "m") ~ "m",
                              str_ends(sample, "n") ~ "n",
                              str_ends(sample, "o") ~ "o",
                              str_ends(sample, "p") ~ "p",
                              ))

# ADD A SAMPLE LOCATION COLUMN AND CLEAN TO REMOVE REPLICATE NAMES SO WE CAN REMOVE SOLOS
replicates$sample_location = replicates$sample
replicates$sample_location <- gsub('_.*', '', replicates$sample_location)

# REMOVE PEAKS WITH RESULTS IN ONLY ONE REPLICATE
soloremoved <- plyr::ddply(replicates, c("unique_id", "sample_location"),
                           function(d) {if (nrow(d) > 1) d else NULL})

# INCLUDE SAMPLE INFO IF NEEDED
sample_info <- readr::read_csv("data/samples/2023-09-naburn-drying-samples.csv")
sample_included <- soloremoved %>% left_join(sample_info, 
                                             by = "sample_location")
sample_included$sample_name <- paste(sample_included$day, 
                                     sample_included$height,
                                     sample_included$length, 
                                      sep = "-")
sample_included <- select(sample_included,
                           -day,
                           -height,
                           -length)

#----
# SPLIT RESULTS BASED ON MASS LIST VS MZCLOUD
split <- split(sample_included, sample_included$annot_source_mass_list_search)
mzcloud <- split$"No results"
masslists <- split$"Full match"

# MERGE MASS LISTS INTO ONE COLUMN
masslistmerged <- masslists %>% 
  pivot_longer(cols = c(starts_with("mass_list_match")) ,
               names_to = "mass_list_name",
               names_prefix = "mass_list_match_",
               values_to = "mass_list_match")

# FILTER FOR NO MATCHES AND INVALID MASS RESULTS
filteredmasslist <- masslistmerged[!grepl('No matches found', 
                                          masslistmerged$mass_list_match),]
filteredmzcloud <- mzcloud[!grepl('Invalid Mass', 
                                  mzcloud$annot_source_mz_cloud_search),]

# SPLIT THE INDIVIDUAL MASS LISTS
# FIRST FIND MASS LIST NAMES (APPEARING IN CONSOLE)
unique(filteredmasslist$mass_list_name)
# THEN SPLIT BY MASS LIST
splitmasslist <- split(filteredmasslist, filteredmasslist$mass_list_name)
antibiotics <- splitmasslist$"antibiotics_itn_msca_answer_160616_w_dtxsi_ds"
metabolites <- splitmasslist$"itnantibiotic_cyp_metabolites"
psychoactive <- splitmasslist$"kps_psychoactive_substances_v2"
pharmaceuticals <- splitmasslist$"kps_pharmaceuticals"

# wide view for samples and fully annotated view
antibiotics_means <- antibiotics %>%
  group_by(pick(peak_number, sample_name)) %>%
  summarise(mean = mean(group_area),
            std = sd(group_area),
            n = length(group_area),
            se = std/sqrt(n))
antibiotics_info <- select(antibiotics, 
                           -replicate,
                           -sample,
                           -peak_rating,
                           -group_area,
                           -sample_location,
                           -sample_name)
antibiotics_info <- unique(antibiotics_info)
antibiotics_annotated <- antibiotics_means %>% left_join(antibiotics_info, 
                                                         by = "peak_number")
antibiotics_annotated <- select(antibiotics_annotated,
                                peak_number,
                                name,
                                formula,
                                calc_mw,
                                m_z,
                                sample_name,
                                mean,
                                std,
                                n,
                                se)
# to see wide
antibiotics_wide <- antibiotics_annotated %>%
  group_by(name) %>%
  pivot_wider(names_from = sample_name, values_from = mean)

# split by day, height, location
day_location <- antibiotics_annotated %>% 
  separate(sample_name, c("day", "height", "location"), sep = "-")

# clean compound names in antibiotics dataset.
day_location$name <- day_location$name %>% 
  fedmatch::clean_strings()

# add antibiotic classes for common antibiotics
class_info <- read.csv("data/antimicrobial_classes.csv")
location_classes <- day_location %>% 
  fuzzy_left_join(class_info,
                  by = c("name" = "name"),
                  match_fun = str_detect)

# replace NA with 'unknown'
location_classes$class <- location_classes$class %>% 
  replace_na('unknown')

# PRODUCE A CSV OF RESULTS
write.csv(antibiotics, 
          "data/processed/2023-naburn-drying/itn_antibiotics.csv", 
          row.names = FALSE)
write.csv(location_classes, 
          "data/processed/2023-naburn-drying/antibiotics_class_organised.csv", 
          row.names = FALSE)
write.csv(metabolites, 
          "data/processed/2023-naburn-drying/itn_metabolites.csv", 
          row.names = FALSE)
write.csv(psychoactive, 
          "data/processed/2023-naburn-drying/psychoactive.csv", 
          row.names = FALSE)
write.csv(pharmaceuticals, 
          "data/processed/2023-naburn-drying/pharmaceuticals.csv", 
          row.names = FALSE)
write.csv(antibiotics_wide, 
          "data/processed/2023-naburn-drying/antibiotics_wide.csv")

# CHANGE HEIGHT AND WIDTH AND MIDPOINT AS NEEDED
antibiotics_annotated %>% 
  ggplot(aes(y = name, 
             x = sample_name, 
             fill = mean)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "turquoise3", 
                       high = "orange", 
                       mid = "yellow", 
                       midpoint = 2e+08) +
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
ggsave("figures/2023-naburn-drying/antibiotics.pdf", width = 15, height = 10)

# for classes as a whole.
location_classes %>% 
  group_by(location, day) %>% 
  ggplot(aes(x = height, y = mean, fill = class)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "height", y = "intensity") +
  facet_grid(location ~ day) +
  theme_ipsum(base_size = 10)
ggsave("figures/bar-classes.png", width = 4, height = 2)

# split the data into location and time studies
location_study <- location_classes[grepl('01|29',location_classes$day),]
time_study <- location_classes[grepl('bottom',location_classes$height),]
time_study <- time_study[grepl('half', time_study$location),]

# location study plots
# Create a list of target antibiotics.
antibiotic_classes <- unique(location_classes$class)

# Create a loop for each target antibiotic.
for (class in antibiotic_classes) {
  
  # Create a subset of the data for the current target antibiotic.
  data <- location_study[location_study$class == class, ]
  
  # Create a heatmap of the data.
  ggplot(data, aes(x = day, y = name.x, fill = mean)) +
    geom_tile() +
    scale_y_discrete(limits = rev) +
    scale_fill_viridis(discrete = F) +
    labs(x = "day", y = "compound", fill = "intensity") +
    facet_grid(location ~ height) +
    theme_ipsum(base_size = 10)
  
  # Save the linegraph to a file.
  ggsave(paste0("figures/heatmaps/location-heatmap-", class, ".png"), width = 7, height = 7)
}

# Create a loop for each target antibiotic.
for (class in antibiotic_classes) {
  
  # Create a subset of the data for the current target antibiotic.
  data2 <- location_study[location_study$class == class, ]
  
  # Create bar graphs of the data.
  data2 %>% 
    group_by(location, height) %>% 
    ggplot(aes(x = day, y = mean, fill = name.x)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_viridis(discrete = T) +
    labs(x = "day", y = "intensity", fill = "compound") +
    facet_grid(location ~ height) +
    theme_ipsum(base_size = 10)
  
  # Save the linegraph to a file.
  ggsave(paste0("figures/bargraph/location-bargraph-", class, ".png"), width = 7, height = 7)
}

# time series plots

# Create a loop for each target antibiotic.
for (class in antibiotic_classes) {
  
  # Create a subset of the data for the current target antibiotic.
  data3 <- time_study[time_study$class == class, ]
  
  # Create a heatmap of the data.
  ggplot(data3, aes(x = day, y = name.x, fill = mean)) +
    geom_tile() +
    scale_y_discrete(limits = rev) +
    scale_fill_viridis(discrete = F) +
    labs(x = "day", y = "compound", fill = "intensity") +
    theme_ipsum(base_size = 10)
  
  # Save the heatmap to a file.
  ggsave(paste0("figures/heatmaps/time-heatmap-", class, ".png"), width = 7, height = 7)
}

# Create a loop for each target antibiotic.
for (class in antibiotic_classes) {
  
  # Create a subset of the data for the current target antibiotic.
  data4 <- time_study[time_study$class == class, ]
  
  # Create line graphs of the data.
  ggplot(data4, aes(x = day, 
                    y = mean, 
                    colour = name.x)) +
    geom_point() +
    geom_errorbar(aes(x = day,
                      ymin = mean - se,
                      ymax = mean + se),
                  width = .6) +
    geom_line(aes(x = day,
                  y = mean,
                  colour =  name.x)) +
    labs(x = "day", y = "intensity", colour = "compound") +
    scale_color_viridis(discrete = TRUE) +
    theme_ipsum(base_size = 10)
  
  # Save the linegraph to a file.
  ggsave(paste0("figures/linegraph/time-error-linegraph-", class, ".png"), width = 7, height = 7)
}
# create new data for unknown categories
unknown_time <- time_study[grepl('unknown', time_study$class),]
unknown_location <- location_study[grepl('unknown', location_study$class),]

unknown_antibiotics <- location_classes[grepl('unknown', location_classes$class),]
write.csv(unknown_antibiotics, 
          "data/processed/2023-naburn-drying/unknown_class_antibiotics.csv", 
          row.names = FALSE)

# graphs as above
# Create a heatmap of the data.
  ggplot(unknown_location, aes(x = day, y = name.x, fill = mean)) +
    geom_tile() +
    scale_y_discrete(limits = rev, labels = scales::label_wrap(40)) +
    scale_fill_viridis(discrete = F) +
    labs(x = "day", y = "compound", fill = "intensity") +
    facet_wrap_paginate(~ height, nrow = 1, ncol = 3, page = 1) +
    theme_ipsum(base_size = 10)
  # Save the linegraph to a file.
  ggsave(paste0("figures/heatmaps/location-heatmap-unknown-location-1.png"), width = 10, height = 10)
  ggplot(unknown_location, aes(x = day, y = name.x, fill = mean)) +
    geom_tile() +
    scale_y_discrete(limits = rev, labels = scales::label_wrap(40)) +
    scale_fill_viridis(discrete = F) +
    labs(x = "day", y = "compound", fill = "intensity") +
    facet_wrap_paginate(~ height, nrow = 1, ncol = 3, page = 2) +
    theme_ipsum(base_size = 10)
  # Save the linegraph to a file.
  ggsave(paste0("figures/heatmaps/location-heatmap-unknown-location-2.png"), width = 10, height = 10)
  ggplot(unknown_location, aes(x = day, y = mean, fill = name.x)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_viridis(discrete = T) +
    labs(x = "day", y = "intensity", fill = "compound") +
    facet_wrap_paginate(~ height, nrow = 1, ncol = 3, page = 1) +
    theme_ipsum(base_size = 10)
  # Save the linegraph to a file.
  ggsave(paste0("figures/bargraph/location-heatmap-unknown-location-1.png"), width = 10, height = 10)
  ggplot(unknown_location, aes(x = day, y = mean, fill = name.x)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_viridis(discrete = T) +
    labs(x = "day", y = "intensity", fill = "compound") +
    facet_wrap_paginate(~ height, nrow = 1, ncol = 3, page = 2) +
    theme_ipsum(base_size = 10)
  # Save the linegraph to a file.
  ggsave(paste0("figures/bargraph/location-heatmap-unknown-location-.png"), width = 10, height = 10)
  # Create a heatmap of the data.
ggplot(unknown_time, aes(x = day, y = name.x, fill = mean)) +
  geom_tile() +
  scale_y_discrete(limits = rev, labels = scales::label_wrap(40)) +
  scale_fill_viridis(discrete = F) +
  labs(x = "day", y = "compound", fill = "intensity") +
  theme_ipsum(base_size = 10)
# Save the heatmap to a file.
ggsave(paste0("figures/heatmaps/time-heatmap-unknown.png"), width = 15, height = 15)
ggplot(unknown_time, aes(x = day, 
                  y = mean, 
                  colour = name.x)) +
  geom_point() +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line(aes(x = day,
                y = mean,
                colour =  name.x)) +
  labs(x = "day", y = "intensity", colour = "compound") +
  scale_color_viridis(discrete = T) +
  theme_ipsum(base_size = 10)
# Save the linegraph to a file.
ggsave(paste0("figures/linegraph/time-error-linegraph-unknown.png"), width = 15, height = 15)

# Create a subset of the data for the current target antibiotic.
trimethoprim <- time_study[grepl('trimethoprim',time_study$class),] 
# Create a heatmap of the data.
ggplot(trimethoprim, aes(x = day, y = name.x, fill = mean)) +
  geom_tile() +
  scale_y_discrete(limits = rev) +
  scale_fill_viridis(discrete = F) +
  labs(x = "day", y = "compound", fill = "intensity") +
  theme_ipsum(base_size = 10)
# Save the heatmap to a file.
ggsave(paste0("figures/heatmaps/time-heatmap-trimethoprim.png"), width = 7, height = 7)
