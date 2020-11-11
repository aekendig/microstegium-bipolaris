##### info ####

# file: background_density_survival_2018_density_exp
# author: Amy Kendig
# date last edited: 11/11/20
# goal: update density values after survival has been accounted for


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
tagDat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
junDat <- read_csv("./data/bg_counts_jun_2018_density_exp.csv")
augDat <- read_csv("./data/bg_ev_counts_late_aug_2018_density_exp.csv")
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")


#### edit data ####

# make survival 1 if the plant produced seeds in summer
# remove NA's 
evSepSurvDat <- tagDat %>%
  filter(month == "September" & sp == "Ev" & focal == 0) %>%
  select(-month) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival)) %>%
  filter(!is.na(survival)) %>%
  group_by(site, plot, treatment, sp, age) %>%
  summarise(counted_sep = n(),
            surviving_sep = sum(survival))

# look at Ev notes
unique(augDat$field_notes)

# update late August counts if they are less than September
evSurvDat <- augDat %>%
  full_join(evSepSurvDat) %>%
  mutate(surviving_sep = replace_na(surviving_sep, 0),
         background_density = ifelse(surviving_sep > green, surviving_sep, green)) %>%
  select(site, plot, treatment, background_density)

# look at Mv notes
unique(junDat$field_notes)

# select Mv data
mvSurvDat <- junDat %>%
  filter(background == "Mv") %>%
  mutate(missing = replace_na(missing, 0),
         background_density = planted - missing - dead) %>%
  select(site, plot, treatment, background_density)

# combine data
plotsD1 <- evSurvDat %>%
  full_join(mvSurvDat) %>%
  full_join(plotsD %>%
              select(-background_density) %>%
              expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
              filter(!(site == "D4" & plot == 8 & treatment == "fungicide") & !(site == "D4" & plot == 10 & treatment == "fungicide"))) %>%
  mutate(background_density = replace_na(background_density, 0)) %>%
  arrange(site, treatment, plot)


#### output ####

write_csv(plotsD1, "intermediate-data/plot_densities_2018_density_exp.csv")
