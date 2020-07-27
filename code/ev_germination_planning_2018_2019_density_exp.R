##### info ####

# file: ev_germination_planning_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 7/22/20
# goal: plan Elymus germination experiment

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
spike18 <- read_csv("./intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv")
spike19 <- read_csv("./intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
surv18 <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
treat <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

#### edit data ####

# focal plants 2019 (all were replaced and tracked)
foc19 <- tibble(site = rep(c("D1", "D2", "D3", "D4"), each = 4),
                ID = rep(c("1", "2", "3", "A"), 4),
                age = rep(c(rep("seedling", 3), "adult"), 4)) %>%
  mutate(sp = "Ev",
         focal = 1) %>%
  merge(treat, all = T) %>%
  as_tibble()

# 2018 plants
plants18 <- surv18 %>%
  filter(month == "September" & !is.na(survival) & !(sp == "Mv" & focal == 0))

# check that all focal are there
plants18 %>%
  filter(focal == 1) %>%
  group_by(sp, age, ID) %>%
  summarise(plots = n())

# focal plants 2018
ev18 <- plants18 %>%
  filter(sp == "Ev" & survival == 1) %>%
  select(site, plot, treatment, sp, ID, age, focal)

# check 2018 data
unique(spike18$ID_unclear)
unique(spike18$spikelet_notes)

# 2018 seed
seed18 <- spike18 %>%
  filter(ID_unclear == 0) %>%
  group_by(site, plot, treatment, sp, ID, focal, age) %>%
  summarise(spikelet_weight.g = sum(spikelet_weight.g),
            seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(ev18) %>%
  left_join(plots) %>%
  left_join(covar) %>%
  left_join(sevy1_1) %>%
  left_join(sevy1_bgev_1) %>%
  mutate(spikelet_weight.g = replace_na(spikelet_weight.g, 0),
         seeds = replace_na(seeds, 0),
         plant_type = paste(sp, age, sep = " "),
         fungicide = recode(treatment, water = 0, fungicide = 1),
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide"),
         seeds_round = round(seeds))

# separate focal
fseed18 <- seed18 %>%
  filter(focal == 1)

# check 2019 data
unique(spike19$spikelet_notes)
filter(spike19, spikelet_notes == "D2?") # does this need to be removed?
unique(spike19$processing_notes)

# 2019 seed
seed19 <- spike19 %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(spikelet_weight.g = sum(spikelet_weight.g),
            seeds = sum(seeds)) %>%
  mutate(age = case_when(ID == "A" ~ "adult",
                         TRUE ~ "seedling"),
         focal = 1) %>%
  ungroup() %>%
  full_join(foc19) %>%
  left_join(covar) %>%
  left_join(bgbio %>%
              rename(background_biomass = biomass.g)) %>%
  left_join(sevy2_1) %>%
  mutate(spikelet_weight.g = replace_na(spikelet_weight.g, 0),
         seeds = replace_na(seeds, 0),
         plant_type = paste(sp, age, sep = " "),
         fungicide = recode(treatment, water = 0, fungicide = 1),
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide"),
         background_pc_bio = background_biomass / background_density)

# combine years
seed <- fseed18 %>%
  mutate(yearf = "Year 1") %>%
  full_join(seed19 %>%
              mutate(yearf = "Year 2"))