##### info ####

# file: biomass_severity_correlations_2019_density_exp
# author: Amy Kendig
# date last edited: 3/15/21
# goal: test correlations of biomass-severity relationships


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
sevChgDat <- read_csv("intermediate-data/severity_change_model_ran_slopes_2018_2019_density_exp.csv")
# severity_change_2018_2019_density_exp.R
envD2Dat <- read_csv("intermediate-data/temp_humidity_monthly_2019_density_exp.csv") 
# temp_humidity_data_processing_2019_density_exp.R
plotBioD2Dat <- read_csv("intermediate-data/plot_biomass_density_2019_density_exp.csv")
# plot_biomass_Density_2019_density_exp.R


#### edit data ####

# plot-level dat
envBioDat <- envD2Dat %>%
  select(site, plot, treatment, month, dew_intensity2) %>%
  rename(dew_intensity = dew_intensity2) %>%
  filter(month == "late_aug") %>%
  left_join(plotBioD2Dat %>%
              mutate(total_biomass = Mv_seedling_biomass + Ev_seedling_biomass + Ev_adult_biomass))

disChgDat <- sevChgDat %>%
  filter(year == 2019 & treatment == "water") %>%
  group_by(plant_group, site, plot) %>%
  summarise(sev_chg = mean(slope)) %>%
  ungroup() %>%
  left_join(envBioDat)


#### correlations ####

cor.test(~ total_biomass + dew_intensity, data = envBioDat)
cor.test(~ total_biomass + sev_chg, data = disChgDat %>% filter(plant_group == "Mv_seedling"))
cor.test(~ total_biomass + sev_chg, data = disChgDat %>% filter(plant_group == "Ev_seedling"))
cor.test(~ total_biomass + sev_chg, data = disChgDat %>% filter(plant_group == "Ev_adult"))
cor.test(~ dew_intensity + sev_chg, data = disChgDat %>% filter(plant_group == "Mv_seedling"))
cor.test(~ dew_intensity + sev_chg, data = disChgDat %>% filter(plant_group == "Ev_seedling"))
cor.test(~ dew_intensity + sev_chg, data = disChgDat %>% filter(plant_group == "Ev_adult"))


#### visualizations ####

ggplot(disChgDat, aes(total_biomass, sev_chg, color = plant_group)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(envBioDat, aes(total_biomass, dew_intensity)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggplot(disChgDat, aes(dew_intensity, sev_chg, color = plant_group)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)
