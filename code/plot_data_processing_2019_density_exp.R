#### info ####

# file: plot_data_processing_2019_density_exp
# author: Amy Kendig
# date last edited: 4/21/21
# goal: combine plot-scale data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
bgBioD2Dat <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv") 
# bg_biomass_data_processing_2019_density_exp.R
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
# mv_seeds_data_processing_2019_density_exp.R
# adj values substitute averages of other two plants in the plot when a plant is missing data
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv") 
# ev_seeds_data_processing_2019.R and ev_seeds_data_processing_2018.R
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") 
# leaf_scans_data_processing_2019_density_exp

evDisMayD2Dat <- read_csv("data/ev_disease_may_2019_density_exp.csv")
fDisJunD2Dat <- read_csv("data/focal_disease_jun_2019_density_exp.csv")
fDisJulD2Dat <- read_csv("data/focal_disease_jul_2019_density_exp.csv")
fDisEAugD2Dat <- read_csv("data/focal_disease_early_aug_2019_density_exp.csv")
fDisLAugD2Dat <- read_csv("data/focal_disease_late_aug_2019_density_exp.csv")


#### Ev seeds ####

# planted density
evPlantsD2Dat <- tibble(plot = 5:10,
                        planted = c(7, 11, 19, 3, 5, 9))

# format seeds
evSeedD2Dat2 <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(evBioD2Dat %>% # complete list of Ev
              select(site, plot, treatment, sp, ID)) %>%
  mutate(seeds = replace_na(seeds, 0),
         age = ifelse(ID == "A", "adult", "seedling")) %>%
  filter((plot %in% 5:7 & age == "seedling") | (plot %in% 8:10 & age == "adult")) %>% # only use the type planted as background
  group_by(site, plot, treatment, age) %>%
  summarise(seeds_per_plant = mean(seeds)) %>%
  ungroup() %>%
  left_join(evPlantsD2Dat) %>%
  mutate(seeds = seeds_per_plant * planted) %>%
  select(-c(seeds_per_plant, planted))

# check
filter(evSeedD2Dat2, is.na(seeds))

# figure
ggplot(evSeedD2Dat2, aes(treatment, seeds)) +
  geom_boxplot() +
  facet_wrap(~ age, scales = "free")


#### Mv seeds ####

# planted density
mvPlantsD2Dat <- tibble(plot = 2:4,
                        planted = c(11, 19, 67))

# format seeds
mvSeedD2Dat2 <- mvSeedD2Dat %>%
  mutate(age = "seedling")  %>%
  filter(plot %in% 2:4) %>%
  group_by(site, plot, treatment, sp, age) %>%
  summarise(seeds_per_plant = mean(seeds, na.rm = T)) %>% # plot average
  ungroup() %>%
  left_join(mvPlantsD2Dat) %>%
  mutate(seeds = seeds_per_plant * planted) %>%
  select(-c(seeds_per_plant, planted))

# check
filter(mvSeedD2Dat2, is.na(seeds))

# figure
ggplot(mvSeedD2Dat2, aes(treatment, seeds)) +
  geom_boxplot()


##### plot biomass ####

# use average of others in plot if plant is missing biomass
plotBioD2Dat <- bgBioD2Dat %>%
  mutate(sp = case_when(plot %in% 2:4 ~ "Mv",
                        plot %in% 5:10 ~ "Ev",
                        TRUE ~ "none"),
         age = case_when(plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling")) %>%
  filter(plot != 1) %>%
  left_join(mvBioD2Dat %>% # add focal biomass
              filter(plot %in% 2:4) %>%
              group_by(site, plot, treatment) %>%
              mutate(biomass_weight_adj.g = mean(biomass_weight.g, na.rm = T)) %>%
              ungroup() %>%
              mutate(biomass_weight.g = case_when(is.na(biomass_weight.g) ~ biomass_weight_adj.g,
                                                  TRUE ~ biomass_weight.g),
                     age = "seedling") %>%
              group_by(site, plot, treatment, sp, age) %>%
              summarise(biomass_foc = sum(biomass_weight.g)) %>%
              ungroup() %>%
              full_join(evBioD2Dat %>%
                          filter(ID %in% c("1", "2", "3") & plot %in% 5:7) %>%
                          group_by(site, plot, treatment) %>%
                          mutate(weight_adj = mean(weight, na.rm = T)) %>%
                          ungroup() %>%
                          mutate(weight = case_when(is.na(weight) ~ weight_adj,
                                                    TRUE ~ weight),
                                 sp = "Ev",
                                 age = "seedling") %>%
                          group_by(site, plot, treatment, sp, age) %>%
                          summarise(biomass_foc = sum(weight)) %>%
                          ungroup()) %>%
              full_join(evBioD2Dat %>%
                          filter(ID == "A" & plot %in% 8:10) %>%
                          select(site, plot, treatment, weight) %>%
                          rename(biomass_foc = weight) %>%
                          mutate(sp = "Ev",
                                 age = "adult"))) %>%
  mutate(biomass.g_m2 = biomass.g + biomass_foc) %>% # plot value
  select(-c(biomass.g, biomass_foc))

# check
filter(plotBioD2Dat, is.na(biomass.g_m2))

# figure
ggplot(plotBioD2Dat, aes(treatment, biomass.g_m2)) +
  geom_boxplot() +
  facet_wrap(sp ~ age, scales = "free")


#### plot severity ####

# leaf counts
plotDisD2Dat <- evDisMayD2Dat %>%
  mutate(month = "may") %>%
  full_join(fDisJunD2Dat %>%
              mutate(month = "jun")) %>%
  full_join(fDisJulD2Dat %>%
              mutate(month = "jul")) %>%
  full_join(fDisEAugD2Dat %>%
              mutate(month = "early_aug")) %>%
  full_join(fDisLAugD2Dat %>%
              mutate(month = "late_aug")) %>%
  filter(!is.na(leaves_tot) & !is.na(leaves_infec)) %>%
  filter((plot %in% 2:4 & sp == "Mv") | (plot %in% 5:10 & sp == "Ev")) %>%
  mutate(age = case_when(ID == "A" ~ "adult",
                         !(ID %in% c("1", "2", "3")) & plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling")) %>%
  group_by(month, site, treatment, plot, sp, age) %>%
  summarise(leaves_tot = mean(leaves_tot),
            leaves_infec = mean(leaves_infec)) %>% # avg leaves infected per tiller for plot
  ungroup()

# check leaves
plotDisD2Dat %>%
  filter(leaves_infec > leaves_tot)

plotDisD2Dat %>%
  filter(leaves_infec == 0) # 49 cases

# severity
plotSevD2Dat <- sevD2Dat %>%
  filter(month != "sep") %>% # too much data missing
  filter((plot %in% 2:4 & sp == "Mv") | (plot %in% 5:10 & sp == "Ev")) %>%
  mutate(age = case_when(ID == "A" ~ "adult",
                         !(ID %in% c("1", "2", "3")) & plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling")) %>%
  select(-c(leaves_tot, leaves_infec)) %>%
  full_join(plotDisD2Dat) %>%
  mutate(leaves_infec = case_when(leaves_infec == 0 & lesion_area.pix > 0 ~ 1, # applies to 66 leaves
                                  TRUE ~ leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = case_when(leaves_infec == 0 & is.na(severity) ~ 0, # 6 have no leaf scans and 0 infected, 33 had infected leaves (leave as NA)
                              TRUE ~ severity)) %>%
  filter(!is.na(severity)) %>%
  group_by(month, site, treatment, plot, sp) %>%
  summarise(severity = mean(severity)) %>% # average over the ages for Ev
  ungroup()

# check
filter(plotSevD2Dat, severity > 1 | is.na(severity))

# figure
ggplot(plotSevD2Dat, aes(x = treatment, y = severity)) +
  geom_boxplot() +
  facet_wrap(sp ~ month, scales = "free")
# highest in late Aug for both

# make wide
plotSevD2DatW <- plotSevD2Dat %>%
  pivot_wider(names_from = month,
              values_from = severity,
              names_glue = "{month}_severity")


##### combine data ####

dat <- evSeedD2Dat2 %>%
  mutate(sp = "Ev") %>%
  full_join(mvSeedD2Dat2) %>%
  full_join(plotBioD2Dat) %>%
  full_join(plotSevD2DatW)

# check values
dat %>%
  group_by(sp, age) %>%
  summarise(plot_types = length(unique(plot)),
            plots = n())


#### output ####

write_csv(dat, "intermediate-data/plot_biomass_seeds_severity_2019_density_exp.csv")
