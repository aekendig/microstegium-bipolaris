##### info ####

# file: mv_water_leaf_incidence.R
# author: Amy Kendig
# date last edited: 2/11/21
# goal: estimate of infection incidence of leaves from the field (for Kendig et al. 2021 Plos ONE)


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import severity data
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
# leaf_scans_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") # leaf_scans_data_processing_2019_density_exp.R


#### edit data ####

# combine data
sevDat <- sevD1Dat %>%
  filter(treatment == "water" & ID %in% c("1", "2", "3") & sp == "Mv") %>%
  mutate(year = 2018) %>%
  select(year, month, site, plot, ID, leaves_tot, leaves_infec) %>%
  full_join(sevD2Dat %>%
              filter(treatment == "water" & ID %in% c("1", "2", "3") & sp == "Mv") %>%
              mutate(year = 2019) %>%
              select(year, month, site, plot, ID, leaves_tot, leaves_infec))

# average leaf infection incidence
sevDat %>%
  group_by(year, month) %>%
  summarise(incidence = sum(leaves_infec)/sum(leaves_tot))

sevDat %>%
  group_by(year) %>%
  summarise(incidence = sum(leaves_infec)/sum(leaves_tot))

sevDat %>%
  summarise(incidence = sum(leaves_infec)/sum(leaves_tot))


#### output ####
write_csv(sevDat, "intermediate-data/mv_water_leaf_incidence.csv")
