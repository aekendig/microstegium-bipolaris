##### info ####

# file: fungal_isolation_data_processing_late_aug_2018_density_exp
# author: Amy Kendig
# date last edited: 11/24/20
# goal: combine two data sheets for fungal isolation data


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
samps <- read_csv("data/fungal_isolation_late_aug_2018_density_exp.csv")
isolates <- read_csv("data/fungal_isolates_late_aug_2018_density_exp.csv")


#### edit data ####

# look at columns
unique(isolates$notes)
unique(isolates$isolate)
unique(samps$gigantea)
unique(samps$plot)
unique(samps$sp)

# create columns
isolates2 <- isolates %>%
  mutate(Pyricularia_isolated = ifelse(str_detect(notes, "pyric") == T, 1, NA),
         small_Bipolaris_isolated = ifelse(str_detect(notes, "small") == T, 1, NA),
         Curvularia_isolated = ifelse(str_detect(notes, "curvularia") == T, 1, NA),
         gigantea_isolated = case_when(str_detect(notes, "pink") == T ~ NA_real_,
                              is.na(Pyricularia_isolated) & is.na(small_Bipolaris_isolated) & is.na(Curvularia_isolated) ~ 1,
                              TRUE ~ NA_real_),
         isolate = gsub("ev-", "", isolate),
         site = toupper(substr(isolate, 1, 2)),
         plot = gsub("[^0-9.]", "", substr(isolate, 3, 4)) %>%
           as.numeric(),
         treatment = ifelse(str_detect(isolate, "w") == T, "water", "fungicide"),
         sp = ifelse(host == "Microstegium", "Mv", "Ev"))

# check numbers
sum(isolates2$gigantea_isolated, na.rm = T)
sum(isolates2$Pyricularia_isolated, na.rm = T) # the "pink" may have been counted as pyricularia in original sheet, but it's unclear what it is
sum(isolates2$small_Bipolaris_isolated, na.rm = T)
sum(isolates2$Curvularia_isolated, na.rm = T)

# summarize by plot
isolates_plot <- isolates2 %>%
  group_by(site, plot, treatment, sp) %>%
  summarise(gigantea_isolated = sum(gigantea_isolated, na.rm = T),
            Pyricularia_isolated = sum(Pyricularia_isolated, na.rm = T),
            small_Bipolaris_isolated = sum(small_Bipolaris_isolated, na.rm = T),
            Curvularia_isolated = sum(Curvularia_isolated, na.rm = T))

# replace zeros with NAs
isolates_plot2 <- isolates_plot
isolates_plot2[isolates_plot2 == 0] <- NA

# combine datasheets
samps2 <- samps %>%
  full_join(isolates_plot2) %>%
  mutate(gigantea = ifelse(is.na(gigantea) & gigantea_isolated > 0, "Yes", gigantea))

# check for inconsistencies
samps2 %>%
  filter(!is.na(gigantea) & !is.na(gigantea_isolated)) %>%
  filter(gigantea == "No" & gigantea_isolated == 1)