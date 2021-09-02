##### info ####

# file: Simon_time_series
# author: Amy Kendig
# date last edited: 9/2/21
# goal: dataset of leaf area, leaf biomass, and temp/humidity over time

# desired data format
# month
# treatment
# plant group
# replicates


# avg temp/humidity/all available 


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
evDisMayD2Dat <- read_csv("data/ev_disease_may_2019_density_exp.csv")
fDisJunD2Dat <- read_csv("data/focal_disease_jun_2019_density_exp.csv")
fDisJulD2Dat <- read_csv("data/focal_disease_jul_2019_density_exp.csv")
fDisEAugD2Dat <- read_csv("data/focal_disease_early_aug_2019_density_exp.csv")
fDisLAugD2Dat <- read_csv("data/focal_disease_late_aug_2019_density_exp.csv")
fDisSepD2Dat <- read_csv("data/focal_disease_sep_2019_density_exp.csv")
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")
# leaf_scans_data_processing_2019_density_exp
edgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
# leaf_scans_data_processing_2019_density_exp
mvLeafD1Dat <- read_csv("data/mv_leaf_weight_2018_density_exp.csv")
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
# leaf_scans_data_processing_2018_density_exp.R
envD2Dat <- read_csv("intermediate-data/temp_humidity_monthly_2019_density_exp.csv") 
# temp_humidity_data_processing_2019_density_exp


#### edit data ####

# date dataframe
dateDat <- tibble(month = c("may", "jun", "jul", "early_aug", "late_aug", "sep", "oct"),
                  date = c("2019-05-07", "2019-06-04", "2019-07-02", "2019-07-29", "2019-08-28", "2019-09-24", "2019-10-21")) %>%
  mutate(date = as.Date(date))

# leaf counts
plotLeavesD2Dat <- evDisMayD2Dat %>%
  mutate(month = "may") %>%
  full_join(fDisJunD2Dat %>%
              mutate(month = "jun")) %>%
  full_join(fDisJulD2Dat %>%
              mutate(month = "jul")) %>%
  full_join(fDisEAugD2Dat %>%
              mutate(month = "early_aug")) %>%
  full_join(fDisLAugD2Dat %>%
              mutate(month = "late_aug")) %>%
  full_join(fDisSepD2Dat %>%
              mutate(month = "sep")) %>%
  filter(!is.na(leaves_tot) & !is.na(leaves_infec) & leaves_tot > 0) %>% # leaves_tot > 0 removes dead plants
  mutate(age = case_when(ID == "A" ~ "adult",
                         !(ID %in% c("1", "2", "3")) & plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling")) %>%
  group_by(month, treatment, sp, age, site, plot) %>%
  summarise(plant_replicates = n(),
            leaves_per_tiller = mean(leaves_tot),
            infected_leaves_per_tiller = mean(leaves_infec)) %>%
  ungroup()

# leaf area
plotAreaD2Dat <- sevD2Dat %>%
  mutate(age = case_when(ID == "A" ~ "adult",
                         !(ID %in% c("1", "2", "3")) & plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling")) %>%
  group_by(month, treatment, sp, age, site, plot) %>%
  summarise(leaf_replicates = n(),
            lesion_area_per_leaf.pix = mean(lesion_area.pix, na.rm = T),
            area_per_leaf.pix = mean(leaf_area.pix, na.rm = T)) %>% # leaf and lesion area averaged by plot
  ungroup()

# edge leaf area
edgeAreaD2Dat <- edgeSevD2Dat %>%
  mutate(age = "edge",
         leaf_replicates = leaf_count,
         lesion_area_per_leaf.pix = lesion_area.pix / leaf_count,
         area_per_leaf.pix = leaf_area.pix / leaf_count) %>%
  select(-c(focal, ID, leaf_area.pix, lesion_area.pix, leaf_count))

# leaf weight
leafWeightD1Dat <- mvLeafD1Dat %>%
  mutate(ID = as.character(ID)) %>%
  left_join(sevD1Dat %>%
              filter(sp == "Mv" & month == "sep") %>%
              select(site, plot, treatment, ID, leaf_area.pix))  %>%
  filter(!is.na(leaf_weight.g) & !(is.na(leaf_area.pix)) & (notes != "Leaf was mutiliated" | is.na(notes))) %>%
  rename(area_per_leaf.pix = leaf_area.pix)

# convert leaf area to weight
ggplot(leafWeightD1Dat, aes(area_per_leaf.pix, leaf_weight.g)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm")

leafMod <- lm(leaf_weight.g ~ area_per_leaf.pix, data = leafWeightD1Dat)
summary(leafMod)

# environmental data
envD2Dat2 <- envD2Dat %>%
  select(month, treatment, site, plot, hum_avg, dew_intensity2, temp_avg, temp_min, temp_max) %>%
  rename("dew_intensity" = "dew_intensity2")

# combine
dat <- plotLeavesD2Dat %>%
  full_join(plotAreaD2Dat) %>%
  full_join(edgeAreaD2Dat) %>%
  full_join(envD2Dat2) %>%
  left_join(dateDat) %>%
  mutate(plant_group = paste(sp, age) %>%
           fct_recode("Mv" = "Mv seedling"),
         weight_per_leaf.g = predict(leafMod, newdata = .),
         weight_per_leaf.g = if_else(sp == "Ev", NA_real_, weight_per_leaf.g)) %>%
  select(month, date, site, plot, treatment, sp, age, plant_group, leaves_per_tiller, infected_leaves_per_tiller, plant_replicates, area_per_leaf.pix, lesion_area_per_leaf.pix, weight_per_leaf.g, leaf_replicates, temp_avg, temp_min, temp_max, hum_avg, dew_intensity)


#### visualize ####

ggplot(dat, aes(date, leaves_per_tiller, col = plant_group, linetype = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  theme_bw()

ggplot(dat, aes(date, area_per_leaf.pix, col = plant_group, linetype = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  theme_bw()

ggplot(dat, aes(date, infected_leaves_per_tiller/leaves_per_tiller, col = plant_group, linetype = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  theme_bw()

ggplot(dat, aes(date, lesion_area_per_leaf.pix/area_per_leaf.pix, col = plant_group, linetype = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  theme_bw()

ggplot(dat, aes(date, weight_per_leaf.g, col = plant_group, linetype = treatment)) +
  stat_summary(geom = "line", fun = "mean") +
  theme_bw()

ggplot(dat, aes(date, hum_avg)) +
  stat_summary(geom = "line", fun = "mean") +
  theme_bw()


#### output ####
write_csv(dat, "intermediate-data/Simon_time_series_data.csv")
