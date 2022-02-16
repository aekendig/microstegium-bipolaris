##### info ####

# file: focal_severity_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 2/16/22
# goal: severity model based on focals


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(cowplot)
library(gridExtra)
library(car)
library(janitor)

# import data
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")
edgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv") 
# leaf_scans_data_processing_2019_density_exp
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
bgBioD2Dat <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv") 
# bg_biomass_data_processing_2019_density_exp.R
evDisMayD2Dat <- read_csv("data/ev_disease_may_2019_density_exp.csv")
fDisJunD2Dat <- read_csv("data/focal_disease_jun_2019_density_exp.csv")
fDisJulD2Dat <- read_csv("data/focal_disease_jul_2019_density_exp.csv")
fDisEAugD2Dat <- read_csv("data/focal_disease_early_aug_2019_density_exp.csv")
fDisLAugD2Dat <- read_csv("data/focal_disease_late_aug_2019_density_exp.csv")

# model functions
source("code/brms_model_fitting_functions.R")


#### edit data ####

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
sevD2Dat2 <- sevD2Dat %>%
  select(month, site, plot, treatment, sp, ID, focal, age, leaf_area.pix, lesion_area.pix) %>%
  full_join(disD2Dat) %>%
  filter(month != "sep" & focal == 1) %>% # too much data missing in sep
  mutate(leaves_infec = case_when(leaves_infec == 0 & lesion_area.pix > 0 ~ 1, # add an infected leaf if scan found lesions
                                  TRUE ~ leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = case_when(leaves_infec == 0 & is.na(severity) ~ 0, # make severity zero if no leaves infected and no severity info
                              TRUE ~ severity)) %>%
  filter(!is.na(severity))

# severity in next month
sevNextD2Dat <- sevD2Dat2 %>%
  filter(!(month %in% c("may", "june"))) %>%
  mutate(month = fct_recode(month,  # match prior month
                            "jun" = "jul",
                            "jul" = "early_aug",
                            "early_aug" = "late_aug")) %>%
  select(month, site, plot, treatment, sp, ID, severity) %>%
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
  full_join(evBioD2Dat2)

# add sp to bg bio
bgBioD2Dat2 <- bgBioD2Dat %>%
  mutate(Mv_bg_bio = ifelse(plot %in% 2:4, biomass.g, 0),
         Ev_bg_bio = ifelse(plot %in% 5:10, biomass.g, 0))

# function to average by plot without focal
sev_neighbor_fun <- function(foc_sp, foc_ID){
  
  # average severity without focal
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
    group_by(site, plot, treatment, sp) %>%
    summarize(foc_bio = sum(biomass_weight.g)) %>%
    ungroup() %>%
    pivot_wider(names_from = sp,
                values_from = foc_bio,
                names_glue = "{sp}_foc_bio") %>%
    full_join(bgBioD2Dat2) %>%
    mutate(Mv_bio = Mv_bg_bio + Mv_foc_bio,
           Ev_bio = Ev_bg_bio + Ev_foc_bio)
  
  dat_out <- sev %>%
    full_join(bio) %>%
    mutate(Mv_sev_bio = Mv_sev * Mv_bio,
           Ev_sev_bio = Ev_sev * Ev_bio,
           sp = foc_sp,
           ID = foc_ID)
  
  return(dat_out)
  
}

# apply function to all focals
sevBioD2Dat <- sev_neighbor_fun("Mv", "1") %>%
  full_join(sev_neighbor_fun("Mv", "2")) %>%
  full_join(sev_neighbor_fun("Mv", "3")) %>%
  full_join(sev_neighbor_fun("Ev", "1")) %>%
  full_join(sev_neighbor_fun("Ev", "2")) %>%
  full_join(sev_neighbor_fun("Ev", "3")) %>%
  full_join(sev_neighbor_fun("Ev", "A"))

# edge data
edgeSevD2Dat2 <- edgeSevD2Dat %>%
  filter(month!= "sep") %>%
  mutate(edge_severity = 100 * (lesion_area.pix / leaf_area.pix),
         edge_severity = ifelse(edge_severity > 100, 100, edge_severity)) %>%
  select(month, site, plot, treatment, edge_severity)

# combine
sevD2Dat3 <- sevD2Dat2 %>%
  filter(month %in% c("jun", "jul", "early_aug")) %>% # only months we have next data for
  full_join(sevNextD2Dat) %>%
  full_join(sevBioD2Dat %>%
              select(month, site, plot, treatment, sp, ID, Mv_sev_bio, Ev_sev_bio)) %>%
  full_join(edgeSevD2Dat2) %>%
  filter(!is.na(severity) & !is.na(next_severity) & !is.na(Mv_sev_bio) & !is.na(Ev_sev_bio) & !is.na(edge_severity)) %>%
  mutate(severity_change = log(next_severity / severity),
         severity_change = case_when(severity == 0 & next_severity == 0 ~ log(1),
                                     severity == 0 & next_severity > 0 ~ log(next_severity / 1e-6), # smallest is 1.3e-6
                                     severity > 0 & next_severity == 0 ~ log(1e-6 / severity),
                                     TRUE ~ severity_change),
         severity = severity * 100,
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         severity_c = severity - mean(severity),
         edge_severity_c = edge_severity - mean(edge_severity),
         Mv_sev_bio_s = (Mv_sev_bio - mean(Mv_sev_bio)) / sd(Mv_sev_bio),
         Ev_sev_bio_s = (Ev_sev_bio - mean(Ev_sev_bio)) / sd(Ev_sev_bio))


#### initial visualizations ####

ggplot(sevD2Dat3, aes(x = severity_c, y = severity_change, color = month)) +
  geom_point() +
  facet_grid(treatment ~ sp)

ggplot(sevD2Dat3, aes(x = Mv_sev_bio_s, y = severity_change, color = month)) +
  geom_point() +
  facet_grid(treatment ~ sp)

ggplot(sevD2Dat3, aes(x = Ev_sev_bio_s, y = severity_change, color = month)) +
  geom_point() +
  facet_grid(treatment ~ sp)

ggplot(sevD2Dat3, aes(x = edge_severity_c, y = severity_change, color = month)) +
  geom_point() +
  facet_grid(treatment ~ sp)


#### fit model ####

sevD2Mod <- brm(severity_change ~ month * sp * fungicide * (severity_c + Mv_sev_bio_s + Ev_sev_bio_s + edge_severity_c) + (1|plotf),
                data = sevD2Dat3, family = gaussian,
                prior <- c(prior(normal(0, 1), class = "Intercept"),
                           prior(normal(0, 1), class = "b")), # use default for sigma
                iter = 6000, warmup = 1000, chains = 3, cores = 3)

mod_check_fun(sevD2Mod)

save(sevD2Mod, "output/focal_severity_model_2019_density_exp.rda")


#### transmission coefficients ####

# Ev early August
beta_Ev_ctrl_aug_self <- "severity_c = 0"
beta_Ev_ctrl_aug_mv <- "Mv_sev_bio_s = 0"
beta_Ev_ctrl_aug_ev <- "Ev_sev_bio_s = 0"
beta_Ev_ctrl_aug_edge <- "edge_severity_c = 0"
beta_Ev_fung_aug_self <- "severity_c + fungicide:severity_c = 0"
beta_Ev_fung_aug_mv <- "Mv_sev_bio_s + fungicide:Mv_sev_bio_s = 0"
beta_Ev_fung_aug_ev <- "Ev_sev_bio_s + fungicide:Ev_sev_bio_s = 0"
beta_Ev_fung_aug_edge <- "edge_severity_c + fungicide:edge_severity_c = 0"

hypothesis(sevD2Mod, c(beta_Ev_ctrl_aug_self, beta_Ev_ctrl_aug_mv, beta_Ev_ctrl_aug_ev, beta_Ev_ctrl_aug_edge,
                       beta_Ev_fung_aug_self, beta_Ev_fung_aug_mv, beta_Ev_fung_aug_ev, beta_Ev_fung_aug_edge))
# self-limiting

# Ev July
beta_Ev_ctrl_jul_self <- "severity_c + monthjul:severity_c = 0"
beta_Ev_ctrl_jul_mv <- "Mv_sev_bio_s + monthjul:Mv_sev_bio_s = 0"
beta_Ev_ctrl_jul_ev <- "Ev_sev_bio_s + monthjul:Ev_sev_bio_s = 0"
beta_Ev_ctrl_jul_edge <- "edge_severity_c + monthjul:edge_severity_c = 0"
beta_Ev_fung_jul_self <- "severity_c + monthjul:severity_c + fungicide:severity_c = 0"
beta_Ev_fung_jul_mv <- "Mv_sev_bio_s + monthjul:Mv_sev_bio_s + fungicide:Mv_sev_bio_s = 0"
beta_Ev_fung_jul_ev <- "Ev_sev_bio_s + monthjul:Ev_sev_bio_s + fungicide:Ev_sev_bio_s = 0"
beta_Ev_fung_jul_edge <- "edge_severity_c + monthjul:edge_severity_c + fungicide:edge_severity_c = 0"

hypothesis(sevD2Mod, c(beta_Ev_ctrl_jul_self, beta_Ev_ctrl_jul_mv, beta_Ev_ctrl_jul_ev, beta_Ev_ctrl_jul_edge,
                       beta_Ev_fung_jul_self, beta_Ev_fung_jul_mv, beta_Ev_fung_jul_ev, beta_Ev_fung_jul_edge))
# self-limiting

# Ev June
beta_Ev_ctrl_jun_self <- "severity_c + monthjun:severity_c = 0"
beta_Ev_ctrl_jun_mv <- "Mv_sev_bio_s + monthjun:Mv_sev_bio_s = 0"
beta_Ev_ctrl_jun_ev <- "Ev_sev_bio_s + monthjun:Ev_sev_bio_s = 0"
beta_Ev_ctrl_jun_edge <- "edge_severity_c + monthjun:edge_severity_c = 0"
beta_Ev_fung_jun_self <- "severity_c + monthjun:severity_c + fungicide:severity_c = 0"
beta_Ev_fung_jun_mv <- "Mv_sev_bio_s + monthjun:Mv_sev_bio_s + fungicide:Mv_sev_bio_s = 0"
beta_Ev_fung_jun_ev <- "Ev_sev_bio_s + monthjun:Ev_sev_bio_s + fungicide:Ev_sev_bio_s = 0"
beta_Ev_fung_jun_edge <- "edge_severity_c + monthjun:edge_severity_c + fungicide:edge_severity_c = 0"

hypothesis(sevD2Mod, c(beta_Ev_ctrl_jun_self, beta_Ev_ctrl_jun_mv, beta_Ev_ctrl_jun_ev, beta_Ev_ctrl_jun_edge,
                       beta_Ev_fung_jun_self, beta_Ev_fung_jun_mv, beta_Ev_fung_jun_ev, beta_Ev_fung_jun_edge))
# self-limiting

# Mv early August
beta_Mv_ctrl_aug_self <- "severity_c + spMv:severity_c = 0"
beta_Mv_ctrl_aug_mv <- "Mv_sev_bio_s + spMv:Mv_sev_bio_s = 0"
beta_Mv_ctrl_aug_ev <- "Ev_sev_bio_s + spMv:Ev_sev_bio_s = 0"
beta_Mv_ctrl_aug_edge <- "edge_severity_c + spMv:edge_severity_c = 0"
beta_Mv_fung_aug_self <- "severity_c + fungicide:severity_c + spMv:severity_c + spMv:fungicide:severity_c = 0"
beta_Mv_fung_aug_mv <- "Mv_sev_bio_s + fungicide:Mv_sev_bio_s + spMv:Mv_sev_bio_s + spMv:fungicide:Mv_sev_bio_s = 0"
beta_Mv_fung_aug_ev <- "Ev_sev_bio_s + fungicide:Ev_sev_bio_s + spMv:Ev_sev_bio_s + spMv:fungicide:Ev_sev_bio_s = 0"
beta_Mv_fung_aug_edge <- "edge_severity_c + fungicide:edge_severity_c + spMv:edge_severity_c + spMv:fungicide:edge_severity_c = 0"

hypothesis(sevD2Mod, c(beta_Mv_ctrl_aug_self, beta_Mv_ctrl_aug_mv, beta_Mv_ctrl_aug_ev, beta_Mv_ctrl_aug_edge,
                       beta_Mv_fung_aug_self, beta_Mv_fung_aug_mv, beta_Mv_fung_aug_ev, beta_Mv_fung_aug_edge))
# self-limiting

# Mv July
beta_Mv_ctrl_jul_self <- "severity_c + spMv:severity_c + monthjul:severity_c + monthjul:spMv:severity_c = 0"
beta_Mv_ctrl_jul_mv <- "Mv_sev_bio_s + spMv:Mv_sev_bio_s + monthjul:Mv_sev_bio_s + monthjul:spMv:Mv_sev_bio_s = 0"
beta_Mv_ctrl_jul_ev <- "Ev_sev_bio_s + spMv:Ev_sev_bio_s + monthjul:Ev_sev_bio_s + monthjul:spMv:Ev_sev_bio_s = 0"
beta_Mv_ctrl_jul_edge <- "edge_severity_c + spMv:edge_severity_c + monthjul:edge_severity_c + monthjul:spMv:edge_severity_c = 0"
beta_Mv_fung_jul_self <- "severity_c + fungicide:severity_c + spMv:severity_c + spMv:fungicide:severity_c + monthjul:severity_c + monthjul:spMv:severity_c + monthjul:fungicide:severity_c + monthjul:spMv:fungicide:severity_c = 0"
beta_Mv_fung_jul_mv <- "Mv_sev_bio_s + fungicide:Mv_sev_bio_s + spMv:Mv_sev_bio_s + spMv:fungicide:Mv_sev_bio_s + monthjul:Mv_sev_bio_s + monthjul:spMv:Mv_sev_bio_s + monthjul:fungicide:Mv_sev_bio_s + monthjul:spMv:fungicide:Mv_sev_bio_s = 0"
beta_Mv_fung_jul_ev <- "Ev_sev_bio_s + fungicide:Ev_sev_bio_s + spMv:Ev_sev_bio_s + spMv:fungicide:Ev_sev_bio_s + monthjul:Ev_sev_bio_s + monthjul:spMv:Ev_sev_bio_s + monthjul:fungicide:Ev_sev_bio_s + monthjul:spMv:fungicide:Ev_sev_bio_s = 0"
beta_Mv_fung_jul_edge <- "edge_severity_c + fungicide:edge_severity_c + spMv:edge_severity_c + spMv:fungicide:edge_severity_c + monthjul:edge_severity_c + monthjul:spMv:edge_severity_c + monthjul:fungicide:edge_severity_c + monthjul:spMv:fungicide:edge_severity_c = 0"

hypothesis(sevD2Mod, c(beta_Mv_ctrl_jul_self, beta_Mv_ctrl_jul_mv, beta_Mv_ctrl_jul_ev, beta_Mv_ctrl_jul_edge,
                       beta_Mv_fung_jul_self, beta_Mv_fung_jul_mv, beta_Mv_fung_jul_ev, beta_Mv_fung_jul_edge))
# self-limiting

# Mv June
beta_Mv_ctrl_jun_self <- "severity_c + spMv:severity_c + monthjun:severity_c + monthjun:spMv:severity_c = 0"
beta_Mv_ctrl_jun_mv <- "Mv_sev_bio_s + spMv:Mv_sev_bio_s + monthjun:Mv_sev_bio_s + monthjun:spMv:Mv_sev_bio_s = 0"
beta_Mv_ctrl_jun_ev <- "Ev_sev_bio_s + spMv:Ev_sev_bio_s + monthjun:Ev_sev_bio_s + monthjun:spMv:Ev_sev_bio_s = 0"
beta_Mv_ctrl_jun_edge <- "edge_severity_c + spMv:edge_severity_c + monthjun:edge_severity_c + monthjun:spMv:edge_severity_c = 0"
beta_Mv_fung_jun_self <- "severity_c + fungicide:severity_c + spMv:severity_c + spMv:fungicide:severity_c + monthjun:severity_c + monthjun:spMv:severity_c + monthjun:fungicide:severity_c + monthjun:spMv:fungicide:severity_c = 0"
beta_Mv_fung_jun_mv <- "Mv_sev_bio_s + fungicide:Mv_sev_bio_s + spMv:Mv_sev_bio_s + spMv:fungicide:Mv_sev_bio_s + monthjun:Mv_sev_bio_s + monthjun:spMv:Mv_sev_bio_s + monthjun:fungicide:Mv_sev_bio_s + monthjun:spMv:fungicide:Mv_sev_bio_s = 0"
beta_Mv_fung_jun_ev <- "Ev_sev_bio_s + fungicide:Ev_sev_bio_s + spMv:Ev_sev_bio_s + spMv:fungicide:Ev_sev_bio_s + monthjun:Ev_sev_bio_s + monthjun:spMv:Ev_sev_bio_s + monthjun:fungicide:Ev_sev_bio_s + monthjun:spMv:fungicide:Ev_sev_bio_s = 0"
beta_Mv_fung_jun_edge <- "edge_severity_c + fungicide:edge_severity_c + spMv:edge_severity_c + spMv:fungicide:edge_severity_c + monthjun:edge_severity_c + monthjun:spMv:edge_severity_c + monthjun:fungicide:edge_severity_c + monthjun:spMv:fungicide:edge_severity_c = 0"

hypothesis(sevD2Mod, c(beta_Mv_ctrl_jun_self, beta_Mv_ctrl_jun_mv, beta_Mv_ctrl_jun_ev, beta_Mv_ctrl_jun_edge,
                       beta_Mv_fung_jun_self, beta_Mv_fung_jun_mv, beta_Mv_fung_jun_ev, beta_Mv_fung_jun_edge))
# self-limiting
# decrease with Ev bio sev