##### info ####

# file: ev-survival-2018
# author: Amy Kendig
# date last edited: 4/7/19
# goal: Ev survival


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# run seed data files
source("./code/ev-seeds-data-processing-2018.R")

# clear everything except seed data
rm(list = setdiff(ls(), "eseeds"))




# Ev to remove
rem_ev <- eseeds %>%
  filter(remove == 1) %>%
  mutate(
    plant = paste(site, plot, treatment, age, ID, sep = ".")
  )

# combine Ev across dates
eseeds2 <- eseeds %>%
  group_by(site, plot, treatment, sp, age, ID, focal) %>%
  summarise(
    seeds = sum(seeds, na.rm = T)
  ) %>%
  mutate(
    plant = paste(site, plot, treatment, age, ID, sep = ".")
  ) %>%
  filter(!(plant %in% rem_ev$plant))