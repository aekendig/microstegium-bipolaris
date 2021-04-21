#### info ####

# file: plot_data_processing_2018_density_exp
# author: Amy Kendig
# date last edited: 4/21/21
# goal: combine plot-scale data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
fSevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
evBgSevD1Dat <- read_csv("intermediate-data/ev_background_leaf_scans_2018_density_exp.csv")
mvBgSevD1Dat <- read_csv("intermediate-data/mv_background_leaf_scans_2018_density_exp.csv")
# all above: leaf_scans_data_processing_2018_density_exp.R

fDisJulD1Dat <- read_csv("data/focal_size_disease_jul_2018_density_exp.csv")
bgDisJulD1Dat <- read_csv("data/bg_disease_jul_2018_density_exp.csv")
fDisAugD1Dat <- read_csv("data/all_disease_seeds_late_aug_2018_density_exp.csv")
bgDisAugD1Dat <- read_csv("data/bg_disease_late_aug_2018_density_exp.csv")
evDisSepD1Dat <- read_csv("data/ev_disease_seeds_sep_2018_density_exp.csv")
mvDisSepD1Dat <- read_csv("data/mv_disease_sep_2018_density_exp.csv")

mvBioD1Dat <- read_csv("intermediate-data/mv_processed_biomass_oct_2018_density_exp.csv")
# mv_biomass_data_processing_2018_density_exp.R
mvSeedsD1Dat <- read_csv("intermediate-data/mv_processed_seeds_2018_density_exp.csv")
# mv_seeds_data_processing_2018_density_exp
evSeedsD1Dat <- read_csv("intermediate-data/ev_processed_seeds_2018_density_exp.csv")
# ev_seeds_data_processing_2018.R
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
# all_survival_data_processing_2018


#### combine leaf counts ####
disD1Dat <- fDisJulD1Dat %>%
  mutate(month = "jul", focal = 1) %>%
  select(month, site, plot, treatment, sp, ID, focal, leaves_tot, leaves_infec) %>%
  full_join(bgDisJulD1Dat %>%
              filter(!is.na(background)) %>%
              mutate(month = "jul", 
                     focal = 0,
                     ID = as.character(tiller),
                     sp = str_replace(background, " A", "")) %>% # make species Ev if Ev A
              select(month, site, plot, treatment, sp, ID, focal, leaves_tot, leaves_infec)) %>%
  full_join(fDisAugD1Dat %>%
              mutate(month = "late_aug", focal = 1) %>%
              select(month, site, plot, treatment, sp, ID, focal, leaves_tot, leaves_infec)) %>%
  full_join(bgDisAugD1Dat %>%
              mutate(month = "late_aug",
                     focal = 0,
                     ID = as.character(tiller),
                     sp = str_replace(background, " A", "")) %>% # make species Ev if Ev A
              select(month, site, plot, treatment, sp, ID, focal, leaves_tot, leaves_infec)) %>%
  full_join(evDisSepD1Dat %>%
              mutate(month = "sep", 
                     focal = case_when(ID %in% c("1", "2", "3", "A") ~ 1,
                                       TRUE ~ 0)) %>%
              select(month, site, plot, treatment, sp, ID, focal, leaves_tot, leaves_infec)) %>%
  full_join(mvDisSepD1Dat %>%
              mutate(month = "sep", 
                     focal = case_when(ID %in% c("1", "2", "3") ~ 1,
                                       TRUE ~ 0)) %>%
              select(month, site, plot, treatment, sp, ID, focal, leaves_tot, leaves_infec)) %>%
  filter(!is.na(leaves_tot))

# check levels
unique(disD1Dat$month)
unique(disD1Dat$site)
unique(disD1Dat$plot)
unique(disD1Dat$treatment)
unique(disD1Dat$sp)
unique(disD1Dat$ID)
unique(disD1Dat$focal)
unique(disD1Dat$leaves_tot)
unique(disD1Dat$leaves_infec)


#### Mv severity ####

# tiller counts by plot
mvDisD1Dat <- disD1Dat %>%
  filter(plot %in% c(2, 3, 4) & sp == "Mv") %>%
  group_by(month, site, treatment, plot) %>%
  summarise(leaves_tot = mean(leaves_tot),
            leaves_infec = mean(leaves_infec)) %>% # avg leaves infected per tiller for plot
  ungroup()

# check leaves
mvDisD1Dat %>%
  filter(leaves_infec > leaves_tot)

mvDisD1Dat %>%
  filter(leaves_infec == 0) # 2 plots in July

# severity
mvSevD1Dat <- mvBgSevD1Dat %>%
  filter(plot %in% c(2, 3, 4)) %>%
  full_join(fSevD1Dat %>%
              filter(plot %in% c(2, 3, 4) & sp == "Mv")) %>%
  group_by(site, treatment, plot) %>%
  summarise(leaf_area.pix = mean(leaf_area.pix / leaf_count), # manually checked leaf counts, need to do because some pictures have multiple leaves
            lesion_area.pix = mean(lesion_area.pix / leaf_count)) %>%
  ungroup() %>%
  full_join(mvDisD1Dat) %>%
  mutate(leaves_infec = case_when(leaves_infec == 0 & lesion_area.pix > 0 ~ 1, # applies to two plots with one leaf each
                                  TRUE ~ leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot)) %>%
  select(month, site, treatment, plot, severity)

# check
filter(mvSevD1Dat, severity > 1 | is.na(severity))

# figure
ggplot(mvSevD1Dat, aes(x = treatment, y = severity)) +
  geom_boxplot() +
  facet_wrap(~ month)

# make wide
mvSevD1DatW <- mvSevD1Dat %>%
  pivot_wider(names_from = month,
              values_from = severity,
              names_glue = "{month}_severity")


#### Mv biomass ####

# g per m2
# mv plots
# make wide
mvBioD1Dat2 <- mvBioD1Dat %>%
  mutate(biomass.g_m2 = bio.g * 0.25 * 0.49) %>%
  filter(plot %in% c(2, 3, 4)) %>%
  select(site, plot, treatment, biomass.g_m2)

# figure
ggplot(mvBioD1Dat2, aes(x = treatment, y = biomass.g_m2)) +
  geom_boxplot()


#### Mv seeds ####

# seeds per m2
# mv plots
# make wide
mvSeedsD1Dat2 <- mvSeedsD1Dat %>%
  mutate(seeds_soil = seeds_soil * (0.05 * 0.25),
         seeds_bio = seeds_bio * (0.49 * 0.25),
         seeds = seeds_soil + seeds_bio) %>%
  filter(plot %in% c(2, 3, 4)) %>%
  select(site, plot, treatment, seeds)

# figure
ggplot(mvSeedsD1Dat2, aes(x = treatment, y = seeds)) +
  geom_boxplot()


#### Ev severity ####

# leaf counts (average by plot because some are missing)
evDisD1Dat <- disD1Dat %>%
  filter(plot %in% 5:10 & sp == "Ev") %>%
  mutate(age = case_when(focal == 1 & ID == "A" ~ "adult",
                         focal == 0 & plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling")) %>%
  group_by(month, site, treatment, plot, age) %>%
  summarise(leaves_tot = mean(leaves_tot),
            leaves_infec = mean(leaves_infec)) %>% # avg leaves infected per tiller for plot
  ungroup()

# check leaves
evDisD1Dat %>%
  filter(leaves_infec > leaves_tot)

evDisD1Dat %>%
  filter(leaves_infec == 0) # 30 plots

# severity
evSevD1Dat <- evBgSevD1Dat %>%
  filter(plot %in% 5:10) %>%
  select(-c(leaves_tot, leaves_infec)) %>%
  full_join(fSevD1Dat %>%
              filter(plot %in% 5:10 & sp == "Ev") %>%
              select(-c(leaves_tot, leaves_infec))) %>%
  full_join(evDisD1Dat) %>%
  mutate(leaves_infec = case_when(leaves_infec == 0 & lesion_area.pix > 0 ~ 1, # applies to 3 plots, 4 leaves
                                  TRUE ~ leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = case_when(leaves_infec == 0 & is.na(severity) ~ 0, # 39 have no leaf scans, 12 had infected leaves (leave as NA)
                              TRUE ~ severity)) %>%
  filter(!is.na(severity)) %>%
  group_by(month, site, treatment, plot) %>%
  summarise(severity = mean(severity)) %>% # average over the ages
  ungroup()
# 138 out of 144 plots (only missing the ones that had a tree fall on them)

# check
filter(evSevD1Dat, severity > 1 | is.na(severity))

# figure
ggplot(evSevD1Dat, aes(x = treatment, y = severity)) +
  geom_boxplot() +
  facet_wrap(~ month)
# highest severity for water in July

# make wide
evSevD1DatW <- evSevD1Dat %>%
  pivot_wider(names_from = month,
              values_from = severity,
              names_glue = "{month}_severity")


#### Ev seeds ####

# need list of all plants
# make survival 1 if the plant produced seeds in summer
evSurvD1Dat <- survD1Dat %>%
  filter(month == "September" & sp == "Ev") %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival)) %>%
  filter(survival == 1) %>%
  select(-c(month, field_notes, seeds_produced, focal))

# need number of plants
evPlantsD1Dat <- bgDisJulD1Dat %>%
  select(plot, planted) %>%
  unique() %>%
  filter(plot %in% 5:10) %>%
  mutate(planted = case_when(plot %in% 5:7 ~ planted + 3,
                             plot %in% 8:10 ~ planted + 1))

# seeds
evSeedsD1Dat2 <- evSeedsD1Dat %>%
  filter(ID_unclear == 0) %>%
  group_by(site, plot, treatment, ID, age) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(evSurvD1Dat) %>%
  mutate(seeds = replace_na(seeds, 0)) %>% # checked that all plants in seeds were in survival
  filter((plot %in% 5:7 & age == "seedling") | (plot %in% 8:10 & age == "adult")) %>% # only use the type planted as background
  group_by(site, plot, treatment, age) %>%
  summarise(seeds_per_plant = mean(seeds)) %>%
  ungroup() %>%
  left_join(evPlantsD1Dat) %>%
  mutate(seeds = seeds_per_plant * planted) %>%
  select(-c(seeds_per_plant, planted))
# missing 5 plots

# figure
ggplot(evSeedsD1Dat2, aes(x = treatment, y = seeds)) +
  geom_boxplot() +
  facet_wrap(~ age, scales = "free")


#### combine data ####

# Mv
mvDat <- mvSevD1DatW %>%
  full_join(mvBioD1Dat2) %>%
  full_join(mvSeedsD1Dat2) %>%
  mutate(sp = "Mv",
         age = "seedling")

# Ev
evDat <- evSevD1DatW %>%
  full_join(evSeedsD1Dat2) %>%
  mutate(sp = "Ev",
         age = case_when(plot %in% 5:7 ~ "seedling",
                         TRUE ~ "adult"))

# both
dat <- mvDat %>%
  full_join(evDat)


#### output ####

write_csv(dat, "intermediate-data/plot_biomass_seeds_severity_2018_density_exp.csv")
