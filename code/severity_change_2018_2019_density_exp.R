##### info ####

# file: severity_change_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/3/21
# goal: change in severity over growing season


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(car) # for logit

# import severity data
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
# leaf_scans_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") # leaf_scans_data_processing_2019_density_exp.R

#### edit data ####

sevD1Dat2 <- sevD1Dat %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity),
         logit_severity = logit(severity, adjust = 0.001),
         prop_leaves = leaves_infec/leaves_tot,
         prop_area = lesion_area.pix/leaf_area.pix,
         month_num = case_when(month == "jul" ~ 0,
                               month == "late_aug" ~ 1,
                               month == "sep" ~ 2))  %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age, sep = "_"),
         plant = paste(plant_group, site, plot, treatment, ID, sep = "_"))

sevD2Dat2 <- sevD2Dat %>%
  filter(!(month %in% c("may", "sep"))) %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity),
         logit_severity = logit(severity, adjust = 0.001),
         prop_leaves = leaves_infec/leaves_tot,
         prop_area = lesion_area.pix/leaf_area.pix,
         month_num = case_when(month == "jun" ~ 0,
                               month == "jul" ~ 1,
                               month == "early_aug" ~ 2,
                               month == "late_aug" ~ 3))  %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age, sep = "_"),
         plant = paste(plant_group, site, plot, treatment, ID, sep = "_"))


#### figure ####

ggplot(sevD1Dat2, aes(x = month_num, y = logit_severity, color = treatment, group = plant)) +
  stat_smooth(method = "lm", se = F, size = 0.5, alpha = 0.5) +
  facet_wrap(~ plant_group)
