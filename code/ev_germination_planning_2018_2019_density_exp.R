##### info ####

# file: ev_germination_planning_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 7/28/20
# goal: plan experimental design for Ev germination trials


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

# import seed model
# load("./output/ev_seeds_data_processing_2018_2019_seed_spikelet_weight_mod.rda")


#### edit data ####

# focal plants 2019 (all were replaced and tracked)
foc19 <- tibble(site = rep(c("D1", "D2", "D3", "D4"), each = 4),
                ID = rep(c("1", "2", "3", "A"), 4),
                age = rep(c(rep("seedling", 3), "adult"), 4)) %>%
  mutate(sp = "Ev",
         focal = 1) %>%
  merge(treat, all = T) %>%
  as_tibble()

4*20*4

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
  mutate(spikelet_weight.g = replace_na(spikelet_weight.g, 0),
         seeds = replace_na(seeds, 0),
         plant_type = paste(sp, age, sep = " "),
         fungicide = recode(treatment, water = 0, fungicide = 1),
         Treatment = recode(as.factor(fungicide), "0" = "control (water)", "1" = "fungicide"))


#### 2018 seeds ####

# no background plots
seed18_plot1 <- seed18 %>% 
  filter(plot == 1) %>%
  group_by(site, treatment, age) %>%
  summarise(seeds = sum(seeds)) %>%
  mutate(year = 2018)


# no background plots and high density plots
seed18_plot1hi <- fseed18 %>% 
  filter(plot %in% c(1, 4, 7, 10)) %>%
  group_by(site, treatment, plot, age) %>%
  summarise(seeds = sum(seeds)) %>%
  mutate(year = 2018)

# all plots
seed18_plotall <- fseed18 %>% 
  group_by(site, treatment, plot, age) %>%
  summarise(seeds = sum(seeds)) %>%
  mutate(year = 2018)


#### 2019 seeds ####

# no background plots
seed19_plot1 <- seed19 %>% 
  filter(plot == 1) %>%
  group_by(site, treatment, age) %>%
  summarise(seeds = sum(seeds)) %>%
  mutate(year = 2019)

# no background plots and high density plots
seed19_plot1hi <- seed19 %>% 
  filter(plot %in% c(1, 4, 7, 10)) %>%
  group_by(site, treatment, plot, age) %>%
  summarise(seeds = sum(seeds)) %>%
  mutate(year = 2019)

# all plots
seed19_plotall <- seed19 %>% 
  group_by(site, treatment, plot, age) %>%
  summarise(seeds = sum(seeds)) %>%
  mutate(year = 2019)


#### combined ####

# no background plots
seed_plot1 <- seed18_plot1 %>%
  full_join(seed19_plot1)

# separate containers needed
filter(seed_plot1, seeds > 0) # 23

# no background and high density plots
seed_plot1hi <- seed18_plot1hi %>%
  full_join(seed19_plot1hi) %>%
  mutate(seed_cat = case_when(seeds == 0 ~ "none",
                              seeds <= 10 ~ "low",
                              seeds > 10 & seeds <= 30 ~ "medium",
                              seeds > 30 ~ "high"))

# separate containers needed
filter(seed_plot1hi, seeds > 0) # 96

# number per group
seed_plot1hi %>%
  group_by(seed_cat) %>%
  count()

# all plots
seed_plotall <- seed18_plotall %>%
  full_join(seed19_plotall)

# separate containers needed
filter(seed_plotall, seeds > 0) # 242