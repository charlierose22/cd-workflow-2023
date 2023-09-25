# PLEASE REVIEW AND CHANGE ANY FUNCTIONAL CODE WRITTEN IN CAPITAL LETTERS
# LOAD TIDYVERSE AND CHECK PACKAGE UPDATES
library(tidyverse)

# IMPORT YOUR CD DATA
drying_study <- readr::read_delim("Naburn Drying/Naburn_Sep2023.csv", 
                                    delim = "\t", trim_ws = TRUE) %>% 
  janitor::clean_names()

# RENAME DATAFRAME FOR CODE TO WORK WITH MINIMAL CHANGES
basedata <- drying_study

# DROP THESE COLUMNS UNLESS THEY HAVE BEEN USED IN YOUR CD WORKFLOW
basedata$tags = NULL
basedata$checked = NULL

# FILTER SAMPLES WITH NO COMPOUND NAMES AND NO MS2 DATA (IF NEEDED)
basedatawithcompoundnames <- with(basedata, basedata[!(name == "" | is.na(name)), ])
ms2dataonly <- basedatawithcompoundnames[!grepl('No MS2', basedatawithcompoundnames$ms2),]

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
  pivot_longer(cols = FIRST_GROUP_COLUMN_ON_LEFT:LAST_PEAK_RATING_COLUMN_ON_RIGHT,
               names_to = "sample",
               values_to = "result")

# CREATE A SAMPLE NAME COLUMN AND FILL, SO WE CAN GROUP PEAK RATING AND GROUP AREA
samplenames <- add_column(longer, measurement = NA)
samplenames <- mutate(samplenames,
                      measurement = case_when(str_detect(sample, "group_area") ~ "group_area",
                                              str_detect(sample, "peak_rating") ~ "peak_rating"))

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
                       str_ends(sample, "a") ~ "1",
                       str_ends(sample, "b") ~ "2",
                       str_ends(sample, "c") ~ "3"))

# ADD A SAMPLE LOCATION COLUMN AND CLEAN TO REMOVE REPLICATE NAMES SO WE CAN REMOVE SOLOS
FilteredReplicate1$sample_location = FilteredReplicate1$sample
FilteredReplicate1$sample_location <- stringi::stri_replace_all_regex(FilteredReplicate1$sample_location, "^\\d|\\d|_*", "")
FilteredReplicate1$sample_location <- gsub('.{1}$', '', FilteredReplicate1$sample_location)

# REMOVE PEAKS WITH RESULTS IN ONLY ONE REPLICATE
soloremoved <- plyr::ddply(replicates, c("unique_id", "sample_location"),
                           function(d) {if (nrow(d) > 1) d else NULL})
#----
# SPLIT RESULTS BASED ON MASS LIST VS MZCLOUD
split <- split(soloremoved, soloremoved$annot_source_mass_list_search)
mzcloud <- split$"No results"
masslists <- split$"Full match"

# MERGE MASS LISTS INTO ONE COLUMN
masslistmerged <- masslists %>% 
  pivot_longer(cols = c(starts_with("mass_list_match")) ,
               names_to = "mass_list_name",
               names_prefix = "mass_list_match_",
               values_to = "mass_list_match")

# FILTER FOR NO MATCHES AND INVALID MASS RESULTS
filteredmasslist <- masslistmerged[!grepl('No matches found', masslistmerged$mass_list_match),]
filteredmzcloud <- mzcloud[!grepl('Invalid Mass', mzcloud$annot_source_mz_cloud_search),]

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