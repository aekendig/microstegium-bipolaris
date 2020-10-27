##### info ####

# file: elymus_adult_gs_survival_model_2018_litter_exp
# author: Amy Kendig
# date last edited: 10/25/20
# goal: estimate Elymus adult growing season survival based on litter amount and whether it was sterilized


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
survL1Dat <- read_csv("intermediate-data/ev_processed_survival_2018_litter_exp.csv")
plotsL <- read_csv("./data/plot_treatments_2018_litter_exp.csv")


#### edit data ####

# edit variables
# remove unnecessary variables
plotsL2 <- plotsL %>%
  mutate(sterilized = case_when(litter == "live" ~ 0,
                                TRUE ~ 1),
         sterilizedF = ifelse(sterilized == 0, "live", "sterilized") %>%
           fct_relevel("sterilized"),
         litter.g.m2 = litter_weight.g,
         litter.g.cm2 = litter.g.m2/10000) %>%
  select(-c(flag_color, justification, litter_weight.g))

# make survival 1 if the plant produced seeds in summer
# remove NA's 
evAGsSurvL1Dat <- survL1Dat %>%
  filter(month == "September") %>%
  select(-month) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, TRUE ~ survival)) %>%
  filter(!is.na(survival)) %>%
  left_join(plotsL2)
# only one plant did not survive out of 27