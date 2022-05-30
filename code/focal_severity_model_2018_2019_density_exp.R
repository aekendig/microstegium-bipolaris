##### info ####

# file: focal_severity_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 2/16/22
# goal: severity model based on focals

### NOTE: this code left non-focal plants in dataset for year 1


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(car)
library(janitor)
library(emmeans)
library(lazerhawk)

# import data
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
# leaf_scans_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")
edgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv") 
# leaf_scans_data_processing_2019_density_exp

mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
bgBioD2Dat <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv") 
# bg_biomass_data_processing_2019_density_exp.R

fDisJulD1Dat <- read_csv("data/focal_size_disease_jul_2018_density_exp.csv")
fDisAugD1Dat <- read_csv("data/all_disease_seeds_late_aug_2018_density_exp.csv")
evDisSepD1Dat <- read_csv("data/ev_disease_seeds_sep_2018_density_exp.csv")
mvDisSepD1Dat <- read_csv("data/mv_disease_sep_2018_density_exp.csv")
evDisMayD2Dat <- read_csv("data/ev_disease_may_2019_density_exp.csv")
fDisJunD2Dat <- read_csv("data/focal_disease_jun_2019_density_exp.csv")
fDisJulD2Dat <- read_csv("data/focal_disease_jul_2019_density_exp.csv")
fDisEAugD2Dat <- read_csv("data/focal_disease_early_aug_2019_density_exp.csv")
fDisLAugD2Dat <- read_csv("data/focal_disease_late_aug_2019_density_exp.csv")

plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
plotsDexp <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")

envD1Dat <- read_csv("intermediate-data/covariates_2018_density_exp.csv") 
# covariate_data_processing_2018_density_exp
envD2Dat <- read_csv("intermediate-data/temp_humidity_monthly_2019_density_exp.csv") 
# temp_humidity_data_processing_2019_density_exp

# model functions
source("code/brms_model_fitting_functions.R")

# function to transform data to account for 0's and 1's
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}


#### edit data ####

# look at notes
filter(fDisJulD1Dat, !is.na(field_notes)) %>%
  data.frame()
filter(fDisJulD1Dat, !is.na(field_notes) & !is.na(leaves_tot)) %>%
  data.frame()
# all dead plants have NA values

filter(bgDisJulD1Dat, !is.na(field_notes))
# one plant has no leaves - remove

filter(fDisAugD1Dat, !is.na(field_notes))

filter(bgDisAugD1Dat, !is.na(field_notes)) %>%
  data.frame()

filter(evDisSepD1Dat, !is.na(field_notes)) %>%
  data.frame()

filter(mvDisSepD1Dat, !is.na(field_notes)) %>%
  data.frame()

disD1Dat <- fDisJulD1Dat %>%
  mutate(month = "jul", focal = 1) %>%
  select(month, site, plot, treatment, sp, ID, focal, leaves_tot, leaves_infec) %>%
  full_join(fDisAugD1Dat %>%
              mutate(month = "late_aug", focal = 1) %>%
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
  filter(!is.na(leaves_tot) & leaves_tot > 0 & focal == 1) %>%
  mutate(age = case_when(ID == "A" ~ "adult",
                         !(ID %in% c("1", "2", "3")) & plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling")) 

# check leaves
disD1Dat %>%
  filter(leaves_infec > leaves_tot)

# look at notes
unique(evDisMayD2Dat$field_notes)
filter(evDisMayD2Dat, !is.na(field_notes)) %>%
  data.frame()

unique(fDisJunD2Dat$field_notes)
filter(fDisJunD2Dat, !is.na(field_notes)) %>%
  data.frame()
# "dead" should have NA values

unique(fDisJulD2Dat$field_notes)
filter(fDisJulD2Dat, !is.na(field_notes)) %>%
  data.frame()
# "dead" should have NA values

unique(fDisEAugD2Dat$field_notes)
filter(fDisEAugD2Dat, !is.na(field_notes) & field_notes != "no 2 lesion leaf") %>%
  data.frame()
# "dead" should have NA values

unique(fDisLAugD2Dat$field_notes)
filter(fDisLAugD2Dat, !is.na(field_notes) & field_notes != "too few green leaves") %>%
  data.frame()
# "appears dead" should have NA values

# leaf counts
disD2Dat <- evDisMayD2Dat %>%
  mutate(month = "may") %>%
  full_join(fDisJunD2Dat %>%
              mutate(month = "jun")) %>%
  full_join(fDisJulD2Dat %>%
              mutate(month = "jul")) %>%
  full_join(fDisEAugD2Dat %>%
              mutate(month = "early_aug")) %>%
  full_join(fDisLAugD2Dat %>%
              mutate(month = "late_aug")) %>%
  filter(!is.na(leaves_tot) & !is.na(leaves_infec) & leaves_tot > 0) %>% # leaves_tot > 0 removes dead plants
  mutate(age = case_when(ID == "A" ~ "adult",
                         !(ID %in% c("1", "2", "3")) & plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling")) %>%
  select(month, site, plot, treatment, sp, ID, age, leaves_tot, leaves_infec)

# check leaves
disD2Dat %>%
  filter(leaves_infec > leaves_tot)

# severity
sevD1Dat2 <- sevD1Dat %>%
  select(month, site, plot, treatment, sp, ID, focal, age, leaf_area.pix, lesion_area.pix) %>%
  full_join(disD1Dat) %>%
  mutate(leaves_infec = case_when(leaves_infec == 0 & lesion_area.pix > 0 ~ 1, # add an infected leaf if scan found lesions
                                  TRUE ~ leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = case_when(leaves_infec == 0 & is.na(severity) ~ 0, # make severity zero if no leaves infected and no severity info
                              TRUE ~ severity),
         plant_group = case_when(sp == "Mv" ~ "Mv",
                                 sp == "Ev" & age == "seedling" ~ "EvS",
                                 TRUE ~ "EvA")) %>%
  filter(!is.na(severity))

sevD2Dat2 <- sevD2Dat %>%
  select(month, site, plot, treatment, sp, ID, focal, age, leaf_area.pix, lesion_area.pix) %>%
  full_join(disD2Dat) %>%
  filter(month != "sep" & focal == 1) %>% # too much data missing in sep
  mutate(leaves_infec = case_when(leaves_infec == 0 & lesion_area.pix > 0 ~ 1, # add an infected leaf if scan found lesions
                                  TRUE ~ leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = case_when(leaves_infec == 0 & is.na(severity) ~ 0, # make severity zero if no leaves infected and no severity info
                              TRUE ~ severity),
         plant_group = case_when(sp == "Mv" ~ "Mv",
                                 sp == "Ev" & age == "seedling" ~ "EvS",
                                 TRUE ~ "EvA")) %>%
  filter(!is.na(severity))

# severity in next month
sevNextD1Dat <- sevD1Dat2 %>%
  filter(month != "jul") %>%
  mutate(month = fct_recode(month,  # match prior month
                            "jul" = "late_aug",
                            "late_aug" = "sep"),
         plant_group = case_when(sp == "Mv" ~ "Mv",
                                 sp == "Ev" & age == "seedling" ~ "EvS",
                                 TRUE ~ "EvA")) %>%
  select(month, site, plot, treatment, sp, ID, plant_group, severity) %>%
  rename(next_severity = severity)

sevNextD2Dat <- sevD2Dat2 %>%
  filter(!(month %in% c("may", "jun"))) %>%
  mutate(month = fct_recode(month,  # match prior month
                            "jun" = "jul",
                            "jul" = "early_aug",
                            "early_aug" = "late_aug"),
         plant_group = case_when(sp == "Mv" ~ "Mv",
                                 sp == "Ev" & age == "seedling" ~ "EvS",
                                 TRUE ~ "EvA")) %>%
  select(month, site, plot, treatment, sp, ID, plant_group, severity) %>%
  rename(next_severity = severity)

# mv biomass
mvBioD2Dat2 <- mvBioD2Dat %>% # add focal biomass
  group_by(site, plot, treatment) %>%
  mutate(biomass_weight_adj.g = mean(biomass_weight.g, na.rm = T)) %>% # replace missing weights with neighbor weights
  ungroup() %>%
  mutate(biomass_weight.g = case_when(is.na(biomass_weight.g) ~ biomass_weight_adj.g,
                                      TRUE ~ biomass_weight.g),
         plant = as.character(plant)) %>%
  select(site, plot, treatment, sp, plant, biomass_weight.g) %>%
  rename(ID = plant)

# ev biomass
evBioD2Dat2 <- evBioD2Dat %>%
  filter(ID %in% c("1", "2", "3")) %>%
  group_by(site, plot, treatment) %>%
  mutate(weight_adj = mean(weight, na.rm = T)) %>% # replace missing weights with neighbor weights
  ungroup() %>%
  mutate(weight = case_when(is.na(weight) ~ weight_adj,
                            TRUE ~ weight)) %>%
  full_join(evBioD2Dat %>%
              filter(ID == "A")) %>%
  select(site, plot, treatment, sp, ID, weight) %>%
  rename(biomass_weight.g = weight)

# combine focal biomass
focBioD2Dat <- mvBioD2Dat2 %>%
  full_join(evBioD2Dat2) %>%
  mutate(age = if_else(ID == "A", "adult", "seedling"),
         plant_group = case_when(sp == "Mv" ~ "Mv",
                                 sp == "Ev" & age == "seedling" ~ "EvS",
                                 TRUE ~ "EvA"))

# add sp to bg bio
bgBioD2Dat2 <- bgBioD2Dat %>%
  mutate(Mv_bg_bio = ifelse(plot %in% 2:4, biomass.g, 0),
         EvS_bg_bio = ifelse(plot %in% 5:7, biomass.g, 0),
         EvA_bg_bio = ifelse(plot %in% 8:10, biomass.g, 0),
         Ev_bg_bio = EvS_bg_bio + EvA_bg_bio)

# all bio
allBioD2Dat <- focBioD2Dat  %>%
  group_by(site, plot, treatment, plant_group) %>%
  summarize(foc_bio = sum(biomass_weight.g)) %>%
  ungroup() %>%
  pivot_wider(names_from = plant_group,
              values_from = foc_bio,
              names_glue = "{plant_group}_foc_bio") %>%
  full_join(bgBioD2Dat2) %>%
  mutate(Mv_bio = Mv_bg_bio + Mv_foc_bio,
         EvS_bio = EvS_bg_bio + EvS_foc_bio,
         EvA_bio = EvA_bg_bio + EvA_foc_bio)

# density
plotDens <- plotsD %>%
  select(plot, treatment, background, background_density) %>%
  mutate(background = str_replace(background, " ", "_")) %>%
  pivot_wider(names_from = background,
              values_from = background_density) %>%
  mutate(Mv_dens = replace_na(Mv_seedling, 0) + 3,
         EvS_dens = replace_na(Ev_seedling, 0) + 3,
         EvA_dens = replace_na(Ev_adult, 0) + 1) %>%
  select(plot, treatment, Mv_dens, EvS_dens, EvA_dens)

# missing columns function
mis_cols <- function(dat, colname){
  
  if(!(colname %in% colnames(dat))) {
    dat[ , colname] <- NA_real_
  }
  
  return(dat)
}

# function to average by plot without focal
sev_neighbor_D1_fun <- function(foc_sp, foc_ID){
  
  # average severity without focal
  # can't do by plant group, because there's only one EvA measurement
  sev <- sevD1Dat2 %>%
    filter(!(sp == foc_sp & ID == foc_ID)) %>%
    group_by(month, site, plot, treatment, sp) %>%
    summarize(mean_sev = mean(severity)) %>%
    ungroup() %>%
    pivot_wider(names_from = sp,
                values_from = mean_sev,
                names_glue = "{sp}_sev")
  
  if(foc_sp == "Mv"){
    dens <- plotDens %>%
      mutate(Mv_dens = Mv_dens - 1)
  } else if(foc_ID == "A") {
    dens <- plotDens %>%
      mutate(EvA_dens = EvA_dens - 1)
  } else {
    dens <- plotDens %>%
      mutate(EvS_dens = EvS_dens - 1) 
  }
  
  dat_out <- sev %>%
    left_join(dens) %>%
    mutate(Mv_sev_dens = Mv_sev * Mv_dens,
           EvS_sev_dens = Ev_sev * EvS_dens,
           EvA_sev_dens = Ev_sev * EvA_dens,
           sp = foc_sp,
           ID = foc_ID)
  
  return(dat_out)
  
}

sev_neighbor_D2_fun <- function(foc_sp, foc_ID){
  
  # average severity without focal
  # can't do by plant group, because there's only one EvA measurement
  sev <- sevD2Dat2 %>%
    filter(!(sp == foc_sp & ID == foc_ID)) %>%
    group_by(month, site, plot, treatment, sp) %>%
    summarize(mean_sev = mean(severity)) %>%
    ungroup() %>%
    pivot_wider(names_from = sp,
                values_from = mean_sev,
                names_glue = "{sp}_sev")
  
  bio <- focBioD2Dat  %>%
    filter(!(sp == foc_sp & ID == foc_ID)) %>%
    group_by(site, plot, treatment, plant_group) %>%
    summarize(foc_bio = sum(biomass_weight.g)) %>%
    ungroup() %>%
    pivot_wider(names_from = plant_group,
                values_from = foc_bio,
                names_glue = "{plant_group}_foc_bio") %>%
    full_join(bgBioD2Dat2) %>%
    mis_cols("EvA_foc_bio") %>%
    mutate(EvA_foc_bio = replace_na(EvA_foc_bio, 0),
           EvS_foc_bio = replace_na(EvS_foc_bio, 0),
           Mv_foc_bio = replace_na(Mv_foc_bio, 0)) %>%
    mutate(Mv_bio = Mv_bg_bio + Mv_foc_bio,
           EvS_bio = EvS_bg_bio + EvS_foc_bio,
           EvA_bio = EvA_bg_bio + EvA_foc_bio)
  
  if(foc_sp == "Mv"){
    dens <- plotDens %>%
      mutate(Mv_dens = Mv_dens - 1)
  } else if(foc_ID == "A") {
    dens <- plotDens %>%
      mutate(EvA_dens = EvA_dens - 1)
  } else {
    dens <- plotDens %>%
      mutate(EvS_dens = EvS_dens - 1) 
  }
  
  dat_out <- sev %>%
    full_join(bio) %>%
    left_join(dens) %>%
    mutate(Mv_sev_bio = Mv_sev * Mv_bio,
           EvS_sev_bio = Ev_sev * EvS_bio,
           EvA_sev_bio = Ev_sev * EvA_bio,
           Mv_sev_dens = Mv_sev * Mv_dens,
           EvS_sev_dens = Ev_sev * EvS_dens,
           EvA_sev_dens = Ev_sev * EvA_dens,
           sp = foc_sp,
           ID = foc_ID)
  
  return(dat_out)
  
}

# apply function to all focals
sevDensD1Dat <- sev_neighbor_D1_fun("Mv", "1") %>%
  full_join(sev_neighbor_D1_fun("Mv", "2")) %>%
  full_join(sev_neighbor_D1_fun("Mv", "3")) %>%
  full_join(sev_neighbor_D1_fun("Ev", "1")) %>%
  full_join(sev_neighbor_D1_fun("Ev", "2")) %>%
  full_join(sev_neighbor_D1_fun("Ev", "3")) %>%
  full_join(sev_neighbor_D1_fun("Ev", "A"))

sevBioD2Dat <- sev_neighbor_D2_fun("Mv", "1") %>%
  full_join(sev_neighbor_D2_fun("Mv", "2")) %>%
  full_join(sev_neighbor_D2_fun("Mv", "3")) %>%
  full_join(sev_neighbor_D2_fun("Ev", "1")) %>%
  full_join(sev_neighbor_D2_fun("Ev", "2")) %>%
  full_join(sev_neighbor_D2_fun("Ev", "3")) %>%
  full_join(sev_neighbor_D2_fun("Ev", "A"))

# edge data
edgeSevD1Dat2 <- envD1Dat %>%
  mutate(edge_severity = 100 * mv_inf_jul.prop,
         edge_severity = ifelse(edge_severity > 100, 100, edge_severity),
         month = "late_aug") %>%
  select(month, site, plot, treatment, edge_severity)

edgeSevD2Dat2 <- edgeSevD2Dat %>%
  filter(month!= "sep") %>%
  mutate(edge_severity = 100 * (lesion_area.pix / leaf_area.pix),
         edge_severity = ifelse(edge_severity > 100, 100, edge_severity)) %>%
  select(month, site, plot, treatment, edge_severity)

# environmental data
envD1Dat2 <- envD1Dat %>%
  mutate(soil_moisture_c = soil_moisture_jun.prop*100 - mean(soil_moisture_jun.prop*100, na.rm = T),
         canopy_cover_c = canopy_cover.prop*100 - mean(canopy_cover.prop*100, na.rm = T)) %>%
  select(site, plot, treatment, soil_moisture_c, canopy_cover_c)

envD2Dat2 <- envD2Dat %>%
  filter(month %in% c("early_aug", "late_aug")) %>%
  mutate(dew_c = (dew_intensity2 - mean(dew_intensity2, na.rm = T)) / sd(dew_intensity2, na.rm = T),
         temp_avg_s = (temp_avg - mean(temp_avg, na.rm = T)) / sd(temp_avg, na.rm = T),
         hum_avg_s = (hum_avg - mean(hum_avg, na.rm = T)) / sd(hum_avg, na.rm = T),
         month = fct_recode(month,
                            jul = "early_aug",
                            early_aug = "late_aug")) %>%
    select(month, site, plot, treatment, dew_c, temp_avg_s, hum_avg_s)

# combine
sevD1Dat3 <- sevD1Dat2 %>%
  filter(month %in% c("jul", "late_aug")) %>% # only months we have next data for
  full_join(sevNextD1Dat) %>%
  full_join(plotDens %>%
              select(plot, treatment, Mv_dens, EvS_dens, EvA_dens)) %>%
  full_join(sevDensD1Dat %>%
              select(month, site, plot, treatment, sp, ID, Mv_sev, Ev_sev, 
                     Mv_sev_dens, EvS_sev_dens, EvA_sev_dens)) %>%
  full_join(edgeSevD1Dat2) %>%
  left_join(envD1Dat2) %>%
  filter(!is.na(severity) & !is.na(next_severity)) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         plant_group = fct_relevel(plant_group, "Mv", "EvS"),
         next_severity_l = logit(next_severity, adjust = 1e-6),
         next_severity_t = transform01(next_severity))

sevD2Dat3 <- sevD2Dat2 %>%
  filter(month %in% c("jun", "jul", "early_aug")) %>% # only months we have next data for
  full_join(sevNextD2Dat) %>%
  full_join(allBioD2Dat %>%
              select(site, plot, treatment, Mv_bio, EvS_bio, EvA_bio)) %>%
  full_join(plotDens %>%
              select(plot, treatment, Mv_dens, EvS_dens, EvA_dens)) %>%
  full_join(sevBioD2Dat %>%
              select(month, site, plot, treatment, sp, ID, Mv_sev, Ev_sev, 
                     Mv_sev_bio, EvS_sev_bio, EvA_sev_bio, Mv_sev_dens, EvS_sev_dens, EvA_sev_dens)) %>%
  full_join(edgeSevD2Dat2) %>%
  left_join(envD2Dat2) %>%
  filter(!is.na(severity) & !is.na(next_severity) & !is.na(edge_severity)) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         plant_group = fct_relevel(plant_group, "Mv", "EvS"),
         next_severity_l = logit(next_severity, adjust = 1e-6),
         next_severity_t = transform01(next_severity))

# data by month
sevD1Dat3_aug <- sevD1Dat3 %>% filter(month == "late_aug") %>%
  mutate(next_severity_t = transform01(next_severity))
sevD1Dat3_jul <- sevD1Dat3 %>% filter(month == "jul") %>%
  mutate(next_severity_t = transform01(next_severity))

sevD2Dat3_aug <- sevD2Dat3 %>% filter(month == "early_aug") %>%
  mutate(next_severity_t = transform01(next_severity))
sevD2Dat3_jul <- sevD2Dat3 %>% filter(month == "jul") %>%
  mutate(next_severity_t = transform01(next_severity))
sevD2Dat3_jun <- sevD2Dat3 %>% filter(month == "jun") %>%
  mutate(next_severity_t = transform01(next_severity))


#### initial visualizations ####

# Mv density
sevD1Dat3 %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

sevD2Dat3 %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

# Mv continuous
sevD1Dat3 %>%
  ggplot(aes(x = Mv_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = Mv_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

# EvS density
sevD1Dat3 %>%
  filter(plot %in% c(1, 5:7)) %>%
  ggplot(aes(x = EvS_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

sevD2Dat3 %>%
  filter(plot %in% c(1, 5:7)) %>%
  ggplot(aes(x = EvS_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

# EvS continuous
sevD1Dat3 %>%
  ggplot(aes(x = EvS_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = EvS_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

# EvA density
sevD1Dat3 %>%
  filter(plot %in% c(1, 8:10)) %>%
  ggplot(aes(x = EvA_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

sevD2Dat3 %>%
  filter(plot %in% c(1, 8:10)) %>%
  ggplot(aes(x = EvA_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

# EvA continuous
sevD1Dat3 %>%
  ggplot(aes(x = EvA_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = EvA_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

# edge
ggplot(sevD1Dat3, aes(x = edge_severity, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

ggplot(sevD2Dat3, aes(x = edge_severity, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

# severity over time
sevD1Dat3 %>%
  ggplot(aes(x = month, y = next_severity, color = plant_group)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = mean_cl_boot,
               position = position_dodge(0.2)) +
  stat_summary(geom = "point", size = 2, fun = mean,
               position = position_dodge(0.2)) +
  facet_wrap(~ treatment)

sevD2Dat3 %>%
  mutate(month = fct_relevel(month, "jun", "jul")) %>%
  ggplot(aes(x = month, y = next_severity, color = plant_group)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = mean_cl_boot,
               position = position_dodge(0.2)) +
  stat_summary(geom = "point", size = 2, fun = mean,
               position = position_dodge(0.2)) +
  facet_wrap(~ treatment)

# edge severity and dens_sev
sevD1Dat3 %>%
  ggplot(aes(x = edge_severity, y = Mv_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD1Dat3 %>%
  ggplot(aes(x = edge_severity, y = EvS_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD1Dat3 %>%
  ggplot(aes(x = edge_severity, y = EvA_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = edge_severity, y = Mv_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = edge_severity, y = EvS_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = edge_severity, y = EvA_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

# logit transformation
sevD2Dat3 %>%
  ggplot(aes(x = next_severity, y = next_severity_l)) +
  geom_point()
# really skews severity changes

sevD2Dat3 %>%
  ggplot(aes(x = next_severity, y = next_severity_t)) +
  geom_point()
# linear change

# environmental variables
sevD1Dat3 %>%
  ggplot(aes(x = soil_moisture_c, y = next_severity_t, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month)

sevD1Dat3 %>%
  ggplot(aes(x = canopy_cover_c, y = next_severity_t, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month)

sevD2Dat3 %>%
  ggplot(aes(x = dew_c, y = next_severity_t)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month)

# temperature/humidity
sevD2Dat3 %>%
  filter(treatment == "water") %>%
  select(site, plot, Mv_bio, temp_avg_s) %>%
  unique() %>%
  ggplot(aes(x = Mv_bio, y = temp_avg_s)) +
  geom_point() +
  geom_smooth(method = "lm")

sevD2Dat3 %>%
  filter(treatment == "water") %>%
  select(site, plot, Mv_bio, hum_avg_s) %>%
  unique() %>%
  ggplot(aes(x = Mv_bio, y = hum_avg_s)) +
  geom_point() +
  geom_smooth(method = "lm")


#### fit models ####

# remove missing data
sevD1Dat3_aug2 <- sevD1Dat3_aug %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens) & !is.na(edge_severity)) %>%
  mutate(next_severity_t = transform01(next_severity))
sevD1Dat3_jul2 <- sevD1Dat3_jul %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens)) %>%
  mutate(next_severity_t = transform01(next_severity))

sevD2Dat3_aug2 <- sevD2Dat3_aug %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens)) %>%
  mutate(next_severity_t = transform01(next_severity))
sevD2Dat3_jul2 <- sevD2Dat3_jul %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens)) %>%
  mutate(next_severity_t = transform01(next_severity))
sevD2Dat3_jun2 <- sevD2Dat3_jun %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens)) %>%
  mutate(next_severity_t = transform01(next_severity))

# fit models
sevD1Mod_sev_dens_jul <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens) + (1|plotf),
                             data = sevD1Dat3_jul2, family = "beta",
                             prior <- c(prior(normal(0, 1), class = "Intercept"),
                                        prior(normal(0, 1), class = "b")), # use default for others
                             iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD1Mod_sev_dens_aug <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + edge_severity) + (1|plotf),
                             data = sevD1Dat3_aug2, family = "beta",
                             prior <- c(prior(normal(0, 1), class = "Intercept"),
                                        prior(normal(0, 1), class = "b")), # use default for others
                             iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD2Mod_sev_dens_aug <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + edge_severity) + (1|plotf),
                         data = sevD2Dat3_aug2, family = "beta",
                         prior <- c(prior(normal(0, 1), class = "Intercept"),
                                    prior(normal(0, 1), class = "b")), # use default for others
                         iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD2Mod_sev_dens_jul <- update(sevD2Mod_sev_dens_aug, newdata = sevD2Dat3_jul2)

sevD2Mod_sev_dens_jun <- update(sevD2Mod_sev_dens_aug, newdata = sevD2Dat3_jun2)

# check models
mod_check_fun(sevD1Mod_sev_dens_jul)
mod_check_fun(sevD1Mod_sev_dens_aug)
mod_check_fun(sevD2Mod_sev_dens_jun)
mod_check_fun(sevD2Mod_sev_dens_jul)
mod_check_fun(sevD2Mod_sev_dens_aug)

# save models
save(sevD1Mod_sev_dens_aug, file = "output/focal_severity_model_aug_2018_dens_exp.rda")
save(sevD1Mod_sev_dens_jul, file = "output/focal_severity_model_jul_2018_dens_exp.rda")

save(sevD2Mod_sev_dens_aug, file = "output/focal_severity_model_aug_2019_dens_exp.rda")
save(sevD2Mod_sev_dens_jul, file = "output/focal_severity_model_jul_2019_dens_exp.rda")
save(sevD2Mod_sev_dens_jun, file = "output/focal_severity_model_jun_2019_dens_exp.rda")

# model summary tables
sevD1_sev_dens_sum <- brms_SummaryTable(sevD1Mod_sev_dens_jul) %>%
  rename(Estimate_16 = Estimate,
         Est.Error_16 = Est.Error,
         "l-95% CI_16" = "l-95% CI",
         "u-95% CI_16" = "u-95% CI") %>%
  full_join(brms_SummaryTable(sevD1Mod_sev_dens_aug) %>%
              rename(Estimate_20 = Estimate,
                     Est.Error_20 = Est.Error,
                     "l-95% CI_20" = "l-95% CI",
                     "u-95% CI_20" = "u-95% CI")) %>%
  mutate(Covariate = str_replace_all(Covariate, "plant_groupEvS", "1st yr comp"),
         Covariate = str_replace_all(Covariate, "plant_groupEvA", "adult comp"),
         Covariate = str_replace_all(Covariate, "Mv_sev_dens", "src invader"),
         Covariate = str_replace_all(Covariate, "EvS_sev_dens", "src 1st yr"),
         Covariate = str_replace_all(Covariate, "EvA_sev_dens", "src adult"),
         Covariate = str_replace_all(Covariate, "edge_severity", "src surr"),
         Covariate = str_replace_all(Covariate, "fungicide", "fung"))
write_csv(sevD1_sev_dens_sum, "output/focal_severity_model_2018_dens_exp.csv")

sevD2_sev_dens_sum <- brms_SummaryTable(sevD2Mod_sev_dens_jun) %>%
  rename(Estimate_8 = Estimate,
         Est.Error_8 = Est.Error,
         "l-95% CI_8" = "l-95% CI",
         "u-95% CI_8" = "u-95% CI") %>%
  full_join(brms_SummaryTable(sevD2Mod_sev_dens_jul) %>%
              rename(Estimate_12 = Estimate,
                     Est.Error_12 = Est.Error,
                     "l-95% CI_12" = "l-95% CI",
                     "u-95% CI_12" = "u-95% CI")) %>%
  full_join(brms_SummaryTable(sevD2Mod_sev_dens_aug) %>%
              rename(Estimate_16 = Estimate,
                     Est.Error_16 = Est.Error,
                     "l-95% CI_16" = "l-95% CI",
                     "u-95% CI_16" = "u-95% CI")) %>%
  mutate(Covariate = str_replace_all(Covariate, "plant_groupEvS", "1st yr comp"),
         Covariate = str_replace_all(Covariate, "plant_groupEvA", "adult comp"),
         Covariate = str_replace_all(Covariate, "Mv_sev_dens", "src invader"),
         Covariate = str_replace_all(Covariate, "EvS_sev_dens", "src 1st yr"),
         Covariate = str_replace_all(Covariate, "EvA_sev_dens", "src adult"),
         Covariate = str_replace_all(Covariate, "edge_severity", "src surr"),
         Covariate = str_replace_all(Covariate, "fungicide", "fung"))
write_csv(sevD2_sev_dens_sum, "output/focal_severity_model_2019_dens_exp.csv")


#### transmission coefficients ####

# marginal trends  
edge_sev_dens_D1_aug <- emtrends(sevD1Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "edge_severity", transform = "response")
Mv_sev_dens_D1_aug <- emtrends(sevD1Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "Mv_sev_dens", transform = "response")
EvS_sev_dens_D1_aug <- emtrends(sevD1Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "EvS_sev_dens", transform = "response")
EvA_sev_dens_D1_aug <- emtrends(sevD1Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "EvA_sev_dens", transform = "response")

Mv_sev_dens_D1_jul <- emtrends(sevD1Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "Mv_sev_dens", transform = "response")
EvS_sev_dens_D1_jul <- emtrends(sevD1Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "EvS_sev_dens", transform = "response")
EvA_sev_dens_D1_jul <- emtrends(sevD1Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "EvA_sev_dens", transform = "response")

edge_sev_dens_aug <- emtrends(sevD2Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "edge_severity", transform = "response")
Mv_sev_dens_aug <- emtrends(sevD2Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "Mv_sev_dens", transform = "response")
EvS_sev_dens_aug <- emtrends(sevD2Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "EvS_sev_dens", transform = "response")
EvA_sev_dens_aug <- emtrends(sevD2Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "EvA_sev_dens", transform = "response")

edge_sev_dens_jul <- emtrends(sevD2Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "edge_severity", transform = "response")
Mv_sev_dens_jul <- emtrends(sevD2Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "Mv_sev_dens", transform = "response")
EvS_sev_dens_jul <- emtrends(sevD2Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "EvS_sev_dens", transform = "response")
EvA_sev_dens_jul <- emtrends(sevD2Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "EvA_sev_dens", transform = "response")

edge_sev_dens_jun <- emtrends(sevD2Mod_sev_dens_jun, c("plant_group", "fungicide"), var = "edge_severity", transform = "response")
Mv_sev_dens_jun <- emtrends(sevD2Mod_sev_dens_jun, c("plant_group", "fungicide"), var = "Mv_sev_dens", transform = "response")
EvS_sev_dens_jun <- emtrends(sevD2Mod_sev_dens_jun, c("plant_group", "fungicide"), var = "EvS_sev_dens", transform = "response")
EvA_sev_dens_jun <- emtrends(sevD2Mod_sev_dens_jun, c("plant_group", "fungicide"), var = "EvA_sev_dens", transform = "response")

sev_dens_D1_coef <- as_tibble(summary(edge_sev_dens_D1_aug, point.est = mean)) %>% 
  mutate(source = "edge", weeks = 20) %>%
  rename(trend = edge_severity.trend) %>%
  full_join(as_tibble(summary(Mv_sev_dens_D1_aug, point.est = mean)) %>% 
                        mutate(source = "Invader (Mv)", weeks = 20) %>%
              rename(trend = Mv_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvS_sev_dens_D1_aug, point.est = mean)) %>% 
                        mutate(source = "1st yr comp. (Ev)", weeks = 20) %>%
              rename(trend = EvS_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvA_sev_dens_D1_aug, point.est = mean)) %>% 
                        mutate(source = "Adult comp. (Ev)", weeks = 20) %>%
              rename(trend = EvA_sev_dens.trend)) %>%
  full_join(as_tibble(summary(Mv_sev_dens_D1_jul, point.est = mean)) %>% 
              mutate(source = "Invader (Mv)", weeks = 16) %>%
              rename(trend = Mv_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvS_sev_dens_D1_jul, point.est = mean)) %>% 
              mutate(source = "1st yr comp. (Ev)", weeks = 16) %>%
              rename(trend = EvS_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvA_sev_dens_D1_jul, point.est = mean)) %>% 
              mutate(source = "Adult comp. (Ev)", weeks = 16) %>%
              rename(trend = EvA_sev_dens.trend)) %>%
  mutate(treatment = if_else(fungicide == 0, "control", "fungicide"),
         plant_group = fct_recode(plant_group, "Invader (Mv)" = "Mv",
                                  "1st yr comp. (Ev)" = "EvS",
                                  "Adult comp. (Ev)" = "EvA"),
         source = fct_relevel(source, "Invader (Mv)") %>%
           fct_recode("Surrounding invader" = "edge"),
         sig = case_when(lower.HPD > 0 & upper.HPD > 0 ~ "yes",
                         lower.HPD < 0 & upper.HPD < 0 ~ "yes",
                         TRUE ~ "no"))

sev_dens_D2_coef <- as_tibble(summary(edge_sev_dens_aug, point.est = mean)) %>% 
  mutate(source = "edge", weeks = 16) %>%
  rename(trend = edge_severity.trend) %>%
  full_join(as_tibble(summary(Mv_sev_dens_aug, point.est = mean)) %>% 
              mutate(source = "Invader (Mv)", weeks = 16) %>%
              rename(trend = Mv_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvS_sev_dens_aug, point.est = mean)) %>% 
              mutate(source = "1st yr comp. (Ev)", weeks = 16) %>%
              rename(trend = EvS_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvA_sev_dens_aug, point.est = mean)) %>% 
              mutate(source = "Adult comp. (Ev)", weeks = 16) %>%
              rename(trend = EvA_sev_dens.trend)) %>%
  full_join(as_tibble(summary(edge_sev_dens_jul, point.est = mean)) %>% 
              mutate(source = "edge", weeks = 12) %>%
              rename(trend = edge_severity.trend)) %>%
  full_join(as_tibble(summary(Mv_sev_dens_jul, point.est = mean)) %>% 
              mutate(source = "Invader (Mv)", weeks = 12) %>%
              rename(trend = Mv_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvS_sev_dens_jul, point.est = mean)) %>% 
              mutate(source = "1st yr comp. (Ev)", weeks = 12) %>%
              rename(trend = EvS_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvA_sev_dens_jul, point.est = mean)) %>% 
              mutate(source = "Adult comp. (Ev)", weeks = 12) %>%
              rename(trend = EvA_sev_dens.trend)) %>%
  full_join(as_tibble(summary(edge_sev_dens_jun, point.est = mean)) %>% 
              mutate(source = "edge", weeks = 8) %>%
              rename(trend = edge_severity.trend)) %>%
  full_join(as_tibble(summary(Mv_sev_dens_jun, point.est = mean)) %>% 
              mutate(source = "Invader (Mv)", weeks = 8) %>%
              rename(trend = Mv_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvS_sev_dens_jun, point.est = mean)) %>% 
              mutate(source = "1st yr comp. (Ev)", weeks = 8) %>%
              rename(trend = EvS_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvA_sev_dens_jun, point.est = mean)) %>% 
              mutate(source = "Adult comp. (Ev)", weeks = 8) %>%
              rename(trend = EvA_sev_dens.trend)) %>%
  mutate(treatment = if_else(fungicide == 0, "control", "fungicide"),
         plant_group = fct_recode(plant_group, "Invader (Mv)" = "Mv",
                                  "1st yr comp. (Ev)" = "EvS",
                                  "Adult comp. (Ev)" = "EvA"),
         source = fct_relevel(source, "Invader (Mv)") %>%
           fct_recode("Surrounding invader" = "edge"),
         sig = case_when(lower.HPD > 0 & upper.HPD > 0 ~ "yes",
                         lower.HPD < 0 & upper.HPD < 0 ~ "yes",
                         TRUE ~ "no"))

# save compiled trends
write_csv(sev_dens_D1_coef, "output/focal_severity_model_2018_dens_exp.csv")
write_csv(sev_dens_D2_coef, "output/focal_severity_model_2019_dens_exp.csv")


#### fungicide effects ####

# absolute difference
# see disease_treatment_effects_2018_2019_density_exp.R for % change

# hypotheses
mv_fung_hyp <- "inv_logit_scaled(Intercept + fungicide)*100 - inv_logit_scaled(Intercept)*100 = 0"
evS_fung_hyp <- "inv_logit_scaled(Intercept + fungicide + plant_groupEvS + plant_groupEvS:fungicide)*100 - inv_logit_scaled(Intercept + plant_groupEvS)*100  = 0"
evA_fung_hyp <- "inv_logit_scaled(Intercept + fungicide + plant_groupEvA + plant_groupEvA:fungicide)*100 - inv_logit_scaled(Intercept + plant_groupEvA)*100  = 0"

hypothesis(sevD2Mod_sev_dens_aug, c(mv_fung_hyp, evS_fung_hyp, evA_fung_hyp)) # Mv
hypothesis(sevD2Mod_sev_dens_jul, c(mv_fung_hyp, evS_fung_hyp, evA_fung_hyp)) # EvS and A
hypothesis(sevD2Mod_sev_dens_jun, c(mv_fung_hyp, evS_fung_hyp, evA_fung_hyp))


#### figures ####

fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0, unit = "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 7,
                                    margin = margin(0, 0, 0.1, 0, unit = "cm")),
        strip.text.y = element_text(size = 7,
                                    margin = margin(0, 0.1, 0, 0, unit = "cm")),
        strip.placement = "outside",
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))

col_pal <- c("#66a61e", "#d95f02", "#7570b3", "#e6ab02")

# figure
pdf("output/focal_severity_figure_2018_density_exp.pdf", width = 5.12, height = 3.94)
ggplot(sev_dens_D1_coef, aes(x = weeks, y = trend*100, color = source)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_errorbar(aes(ymin = lower.HPD*100, ymax = upper.HPD*100), 
                position = position_dodge(1), width = 0, size =0.35) +
  geom_point(aes(fill = source, shape = sig), position = position_dodge(1), size = 2) +
  facet_grid(treatment ~ plant_group, scales = "free") +
  scale_shape_manual(values = c(1, 19), guide = "none") +
  scale_color_manual(values = col_pal, name = "Disease source") +
  scale_fill_manual(values = col_pal, name = "Disease source") +
  scale_x_continuous(breaks = c(16, 20)) +
  labs(x = "Weeks post planting", y = "Change in disease severity (%)") +
  fig_theme
dev.off()

pdf("output/focal_severity_figure_2019_density_exp.pdf", width = 5.12, height = 3.54)
ggplot(sev_dens_D2_coef, aes(x = weeks, y = trend * 100, color = source)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_errorbar(aes(ymin = lower.HPD * 100, ymax = upper.HPD * 100), 
                position = position_dodge(2.25), width = 0, size = 0.35) +
  geom_point(aes(fill = source, shape = sig), position = position_dodge(2.25), size = 2) +
  facet_grid(treatment ~ plant_group, scales = "free") +
  scale_shape_manual(values = c(1, 19), guide = "none") +
  scale_color_manual(values = col_pal, name = "Disease source") +
  scale_fill_manual(values = col_pal, name = "Disease source") +
  scale_x_continuous(breaks = c(8, 12, 16)) +
  labs(x = "Weeks post planting", y = "Change in disease severity (%)") +
  fig_theme
dev.off()


#### environmental models ####

# remove missing data
sevD1Dat3_aug_env <- sevD1Dat3_aug %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens) & !is.na(edge_severity) & !is.na(soil_moisture_c) & !is.na(canopy_cover_c)) %>%
  mutate(next_severity_t = transform01(next_severity))
sevD1Dat3_jul_env <- sevD1Dat3_jul %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens) & !is.na(soil_moisture_c) & !is.na(canopy_cover_c)) %>%
  mutate(next_severity_t = transform01(next_severity))

sevD2Dat3_aug_env <- sevD2Dat3_aug %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens) & !is.na(edge_severity) & !is.na(dew_c)) %>%
  mutate(next_severity_t = transform01(next_severity))
sevD2Dat3_jul_env <- sevD2Dat3_jul %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens) & !is.na(edge_severity) & !is.na(dew_c)) %>%
  mutate(next_severity_t = transform01(next_severity))

# fit models
sevD1Mod_env_jul <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + soil_moisture_c + canopy_cover_c) + (1|plotf),
                             data = sevD1Dat3_jul_env, family = "beta",
                             prior <- c(prior(normal(0, 1), class = "Intercept"),
                                        prior(normal(0, 1), class = "b")), # use default for others
                             iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD1Mod_env_aug <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + edge_severity + soil_moisture_c + canopy_cover_c) + (1|plotf),
                             data = sevD1Dat3_aug_env, family = "beta",
                             prior <- c(prior(normal(0, 1), class = "Intercept"),
                                        prior(normal(0, 1), class = "b")), # use default for others
                             iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD2Mod_env_aug <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + edge_severity + dew_c) + (1|plotf),
                             data = sevD2Dat3_aug_env, family = "beta",
                             prior <- c(prior(normal(0, 1), class = "Intercept"),
                                        prior(normal(0, 1), class = "b")), # use default for others
                             iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD2Mod_env_jul <- update(sevD2Mod_env_aug, newdata = sevD2Dat3_jul_env)

# check models
mod_check_fun(sevD1Mod_env_aug)
mod_check_fun(sevD1Mod_env_jul)
mod_check_fun(sevD2Mod_env_aug)
mod_check_fun(sevD2Mod_env_jul)

# save models
save(sevD1Mod_env_aug, file = "output/focal_severity_model_env_aug_2018_dens_exp.rda")
save(sevD1Mod_env_jul, file = "output/focal_severity_model_env_jul_2018_dens_exp.rda")
save(sevD2Mod_env_aug, file = "output/focal_severity_model_env_aug_2019_dens_exp.rda")
save(sevD2Mod_env_jul, file = "output/focal_severity_model_env_jul_2019_dens_exp.rda")

# marginal trends  
soil_env_D1_aug <- emtrends(sevD1Mod_env_aug, c("plant_group", "fungicide"), var = "soil_moisture_c", transform = "response")
cano_env_D1_aug <- emtrends(sevD1Mod_env_aug, c("plant_group", "fungicide"), var = "canopy_cover_c", transform = "response")
soil_env_D1_jul <- emtrends(sevD1Mod_env_jul, c("plant_group", "fungicide"), var = "soil_moisture_c", transform = "response")
cano_env_D1_jul <- emtrends(sevD1Mod_env_jul, c("plant_group", "fungicide"), var = "canopy_cover_c", transform = "response")

dew_env_D2_aug <- emtrends(sevD2Mod_env_aug, "plant_group", var = "dew_c", transform = "response")
dew_env_D2_jul <- emtrends(sevD2Mod_env_jul, "plant_group", var = "dew_c", transform = "response")

# combine
env_D1_coef <- as_tibble(summary(soil_env_D1_aug, point.est = mean)) %>% 
  mutate(Variable = "soil moisture", Weeks = 20) %>%
  rename(trend = soil_moisture_c.trend) %>%
  full_join(as_tibble(summary(cano_env_D1_aug, point.est = mean)) %>% 
              mutate(Variable = "canopy cover", Weeks = 20) %>%
              rename(trend = canopy_cover_c.trend)) %>%
  full_join(as_tibble(summary(soil_env_D1_jul, point.est = mean)) %>% 
              mutate(Variable = "soil moisture", Weeks = 16) %>%
              rename(trend = soil_moisture_c.trend)) %>%
  full_join(as_tibble(summary(cano_env_D1_jul, point.est = mean)) %>% 
              mutate(Variable = "canopy cover", Weeks = 16) %>%
              rename(trend = canopy_cover_c.trend)) %>%
  mutate(Treatment = if_else(fungicide == 0, "control", "fungicide"),
         Focal = fct_recode(plant_group, "invader" = "Mv",
                                  "first-year comp." = "EvS",
                                  "adult comp." = "EvA"),
         Estimate = trend * 100,
         l95 = lower.HPD * 100,
         u95 = upper.HPD * 100) %>%
  select(Weeks, Focal, Variable, Treatment, Estimate, l95, u95) %>%
  arrange(Weeks, Variable, Focal, Treatment)

env_D2_coef <- as_tibble(summary(dew_env_D2_aug, point.est = mean)) %>% 
  mutate(Weeks = 16) %>%
  full_join(as_tibble(summary(dew_env_D2_jul, point.est = mean)) %>% 
              mutate(Weeks = 12)) %>%
  mutate(Focal = fct_recode(plant_group, "invader" = "Mv",
                            "first-year comp." = "EvS",
                            "adult comp." = "EvA"),
         Estimate = dew_c.trend * 100,
         l95 = lower.HPD * 100,
         u95 = upper.HPD * 100) %>%
  select(Weeks, Focal, Estimate, l95, u95) %>%
  arrange(Weeks, Focal)

# save
write_csv(env_D1_coef, "output/focal_severity_model_env_coefficients_2018_density_exp.csv")
write_csv(env_D2_coef, "output/focal_severity_model_env_coefficients_2019_density_exp.csv")

# check sig relationships
sevD1Dat3_jul_env %>%
  filter(plant_group == "EvS" & treatment == "fungicide") %>%
  ggplot(aes(x = soil_moisture_c, y = next_severity_t)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)

sevD2Dat3_jul_env %>%
  filter(plant_group == "EvS") %>%
  ggplot(aes(x = dew_c, y = next_severity_t)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)
# neither are driven by outliers

# check models with temp and humidity separate
sevD2Mod_temp_hum_aug <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + edge_severity + temp_avg_s + hum_avg_s) + (1|plotf),
                        data = sevD2Dat3_aug_env, family = "beta",
                        prior <- c(prior(normal(0, 1), class = "Intercept"),
                                   prior(normal(0, 1), class = "b")), # use default for others
                        iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)
sevD2Mod_temp_hum_jul <- update(sevD2Mod_temp_hum_aug, newdata = sevD2Dat3_jul_env)

summary(sevD2Mod_temp_hum_aug)
summary(sevD2Mod_temp_hum_jul)
# neither have statistically significant effects


#### water effect - compare to edge ####

# summarize disease severity by plot
sevD2Dat_sum <- sevD2Dat2 %>%
  group_by(month, site, plot, treatment, sp) %>%
  summarize(severity = mean(severity)) %>%
  ungroup() %>%
  filter(sp == "Mv") %>%
  inner_join(edgeSevD2Dat2 %>%
               mutate(edge_severity = edge_severity / 100)) %>%
  pivot_longer(cols = c(severity, edge_severity),
               names_to = "plant_type",
               values_to = "severity") %>%
  mutate(plant_type = if_else(plant_type == "edge_severity", "surrounding", "inside") %>%
           fct_relevel("inside"),
         severity_t = transform01(severity),
         month = fct_relevel(month, "jun", "jul", "early_aug"),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         treatment = fct_relevel(treatment, "water"))

# visualize
ggplot(sevD2Dat_sum, aes(x = month, y = severity_t, color = plant_type)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot", position = position_dodge(1)) +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(1)) +
  facet_wrap(~ treatment)

# model
ctrlD2Mod <- brm(severity_t ~ plant_type*treatment*month + (1|plotf),
                 data = sevD2Dat_sum, family = Beta,
                 prior <- c(prior(normal(0, 1), class = "Intercept"),
                            prior(normal(0, 1), class = "b")),
                 iter = 6000, warmup = 1000, chains = 3, cores = 3)

mod_check_fun(ctrlD2Mod)

# save model
save(ctrlD2Mod, file = "output/edge_experiment_mv_model_2019_density_exp.rda")

# do experimental plants have different severity than edge?
jun_ctrl_comp = "inv_logit_scaled(Intercept + plant_typesurrounding)*100 = inv_logit_scaled(Intercept)*100"
jun_fung_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + treatmentfungicide + plant_typesurrounding:treatmentfungicide)*100 = inv_logit_scaled(Intercept + treatmentfungicide)*100"
jul_ctrl_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + monthjul + plant_typesurrounding:monthjul)*100 = inv_logit_scaled(Intercept + monthjul)*100"
jul_fung_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + treatmentfungicide + monthjul + plant_typesurrounding:treatmentfungicide + plant_typesurrounding:monthjul + treatmentfungicide:monthjul + plant_typesurrounding:treatmentfungicide:monthjul)*100 = inv_logit_scaled(Intercept + treatmentfungicide + monthjul + treatmentfungicide:monthjul)*100"
eau_ctrl_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + monthearly_aug + plant_typesurrounding:monthearly_aug)*100 = inv_logit_scaled(Intercept + monthearly_aug)*100"
eau_fung_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + treatmentfungicide + monthearly_aug + plant_typesurrounding:treatmentfungicide + plant_typesurrounding:monthearly_aug + treatmentfungicide:monthearly_aug + plant_typesurrounding:treatmentfungicide:monthearly_aug)*100 = inv_logit_scaled(Intercept + treatmentfungicide + monthearly_aug + treatmentfungicide:monthearly_aug)*100"
lau_ctrl_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + monthlate_aug + plant_typesurrounding:monthlate_aug)*100 = inv_logit_scaled(Intercept + monthlate_aug)*100"
lau_fung_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + treatmentfungicide + monthlate_aug + plant_typesurrounding:treatmentfungicide + plant_typesurrounding:monthlate_aug + treatmentfungicide:monthlate_aug + plant_typesurrounding:treatmentfungicide:monthlate_aug)*100 = inv_logit_scaled(Intercept + treatmentfungicide + monthlate_aug + treatmentfungicide:monthlate_aug)*100"

# estimates
ctrlD2hyps <- hypothesis(ctrlD2Mod, c(jun_ctrl_comp, jun_fung_comp, jul_ctrl_comp, jul_fung_comp,
                                      eau_ctrl_comp, eau_fung_comp, lau_ctrl_comp, lau_fung_comp))

ctrlD2hypsOut <- ctrlD2hyps[[1]] %>%
  as_tibble() %>%
  mutate(Weeks = rep(c(4, 8, 12, 16), each = 2),
         Treatment = rep(c("control", "fungicide"), 4))

write_csv(ctrlD2hypsOut, "output/edge_experiment_mv_2019_density_exp.csv")
