##### info ####

# file: severity_change_density_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/4/21
# goal: analyses of severity change as a function of density


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import severity change data
sevChgD1Dat <- read_csv("intermediate-data/severity_change_model_ran_slopes_2018_density_exp.csv")
sevChgD2Dat <- read_csv("intermediate-data/severity_change_model_ran_slopes_2019_density_exp.csv")
# both: severity_change_2018_2019_density_exp.R
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
# leaf_scans_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv") 
# leaf_scans_data_processing_2019_density_exp.R

# biomass
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
bgBioD2Dat <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv") 
# bg_biomass_data_processing_2019_density_exp.R

# function to transform data to account for 0's and 1's
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}

backtransform01 <- function(x) {
  (x * length(x) - 0.5) / (length(x) - 1)
}  

# function to transform logit to proportion
logit2prop <- function(x) {
  (exp(x) / (1 + exp(x)))
}


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(Mv_seedling_density = case_when(background == "Mv seedling" ~ background_density + 3, 
                                         TRUE ~ 3),
         Ev_seedling_density = case_when(background == "Ev seedling" ~ background_density + 3, 
                                         TRUE ~ 3),
         Ev_adult_density = case_when(background == "Ev adult" ~ background_density + 1, 
                                      TRUE ~ 1)) %>%
  select(plot, treatment, Mv_seedling_density, Ev_seedling_density, Ev_adult_density)

# plot biomass
# use average of other species in plot if plant is missing biomass
plotBioD2Dat <- bgBioD2Dat %>%
  mutate(sp = case_when(plot %in% 2:4 ~ "Mv",
                        plot %in% 5:10 ~ "Ev",
                        TRUE ~ "none"),
         age = case_when(plot %in% 8:10 ~ "adult",
                         TRUE ~ "seedling"),
         sp_age = paste(sp, age, sep = "_")) %>%
  pivot_wider(-c(sp, age),
              names_from = sp_age,
              names_glue = "{sp_age}_biomass",
              values_from = biomass.g) %>%
  mutate(Mv_seedling_biomass = replace_na(Mv_seedling_biomass, 0),
         Ev_seedling_biomass = replace_na(Ev_seedling_biomass, 0),
         Ev_adult_biomass = replace_na(Ev_adult_biomass, 0)) %>%
  select(-none_seedling_biomass) %>%
  left_join(mvBioD2Dat %>%
              group_by(site, plot, treatment) %>%
              mutate(biomass_weight_adj.g = mean(biomass_weight.g, na.rm = T)) %>%
              ungroup() %>%
              mutate(biomass_weight.g = case_when(is.na(biomass_weight.g) ~ biomass_weight_adj.g,
                                                  TRUE ~ biomass_weight.g)) %>%
              group_by(site, plot, treatment, ) %>%
              summarise(Mv_seedling_biomass_foc = sum(biomass_weight.g)) %>%
              ungroup()) %>%
  left_join(evBioD2Dat %>%
              filter(ID %in% c("1", "2", "3")) %>%
              group_by(site, plot, treatment) %>%
              mutate(weight_adj = mean(weight, na.rm = T)) %>%
              ungroup() %>%
              mutate(weight = case_when(is.na(weight) ~ weight_adj,
                                        TRUE ~ weight)) %>%
              group_by(site, plot, treatment) %>%
              summarise(Ev_seedling_biomass_foc = sum(weight)) %>%
              ungroup()) %>%
  left_join(evBioD2Dat %>%
              filter(ID == "A") %>%
              select(site, plot, treatment, weight) %>%
              rename(Ev_adult_biomass_foc = weight)) %>%
  mutate(Mv_seedling_biomass = Mv_seedling_biomass + Mv_seedling_biomass_foc,
         Ev_seedling_biomass = Ev_seedling_biomass + Ev_seedling_biomass_foc,
         Ev_adult_biomass = Ev_adult_biomass + Ev_adult_biomass_foc) %>%
  select(-c(Mv_seedling_biomass_foc, Ev_seedling_biomass_foc, Ev_adult_biomass_foc))

# severity change
# parse severity plant IDs
sevChgD1Dat2 <- sevChgD1Dat %>%
  rowwise() %>%
  mutate(sp = strsplit(plant, "_")[[1]][1],
         age = strsplit(plant, "_")[[1]][2],
         site = strsplit(plant, "_")[[1]][3],
         plot = strsplit(plant, "_")[[1]][4] %>%
           as.numeric(),
         treatment = strsplit(plant, "_")[[1]][5],
         ID = strsplit(plant, "_")[[1]][6]) %>%
  ungroup() %>%
  rename(sev_chg = slope) %>%
  select(-c(error, lower, upper))

sevChgD2Dat2 <- sevChgD2Dat %>%
  rowwise() %>%
  mutate(sp = strsplit(plant, "_")[[1]][1],
         age = strsplit(plant, "_")[[1]][2],
         site = strsplit(plant, "_")[[1]][3],
         plot = strsplit(plant, "_")[[1]][4] %>%
           as.numeric(),
         treatment = strsplit(plant, "_")[[1]][5],
         ID = strsplit(plant, "_")[[1]][6]) %>%
  ungroup() %>%
  rename(sev_chg = slope) %>%
  select(-c(error, lower, upper))

# raw severity
sevD1Dat2 <- sevD1Dat %>%
  filter(month == "sep") %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity),
         severity_01 = transform01(severity)) %>%
  select(site, plot, treatment, sp, age, ID, severity, severity_01)

sevD2Dat2 <- sevD2Dat %>%
  filter(month == "late_aug") %>%
  mutate(severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = ifelse(severity > 1, 1, severity),
         severity_01 = transform01(severity)) %>%
  select(site, plot, treatment, sp, age, ID, severity, severity_01)

# combine and separate data
mvD1Dat <- plotDens %>%
  filter(plot %in% 1:4) %>%
  left_join(sevChgD1Dat2 %>%
              full_join(sevD1Dat2)) %>%
  mutate(plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         fungicide = ifelse(treatment == "fungicide", 1, 0))

mvD2Dat <- plotDens %>%
  filter(plot %in% 1:4) %>%
  left_join(plotBioD2Dat) %>%
  left_join(sevChgD2Dat2 %>%
              full_join(sevD2Dat2)) %>%
  mutate(plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         fungicide = ifelse(treatment == "fungicide", 1, 0))


#### visualizations ####

# 2018 Mv density
ggplot(mvD1Dat, aes(Mv_seedling_density, sev_chg, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

ggplot(mvD1Dat, aes(Mv_seedling_density, severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

# 2019 Mv density
ggplot(mvD2Dat, aes(Mv_seedling_density, sev_chg, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

ggplot(mvD2Dat, aes(Mv_seedling_density, severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

# 2019 Mv density, fungicide only
ggplot(mvD2Dat %>% filter(treatment == "water"), aes(Mv_seedling_density, severity, color = plant_group)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  theme_bw()

# 2019 Mv biomass
ggplot(mvD2Dat, aes(Mv_seedling_biomass, sev_chg, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ plant_group, scales = "free")

ggplot(mvD2Dat, aes(Mv_seedling_biomass, severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ plant_group, scales = "free")

# 2019 site variation
ggplot(mvD2Dat, aes(site, sev_chg, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

ggplot(mvD2Dat, aes(site, severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

ggplot(plotBioD2Dat %>% filter(plot %in% 1:4), aes(site, Mv_seedling_biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)


#### Mv biomass model ####

# remove missing data
# add plot ID
# scale biomass
mvD2Dat2 <- mvD2Dat %>%
  filter(!is.na(severity)) %>%
  mutate(plot_ID = paste(site, plot, substr(treatment, 1, 1), sep = "_"),
         bio_scaled = scale(Mv_seedling_biomass))

# intercept prior
mvD2Dat2 %>%
  filter(treatment == "water" & plant_group == "Mv_seedling") %>%
  summarise(sev = mean(severity_01))
# 0.34

# 2019 model
mvSevD2Mod <- brm(severity_01 ~ bio_scaled * fungicide * plant_group + (1|plot_ID), 
                data = mvD2Dat2, family = "beta",
                prior <- c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 1), class = b)), # use default priors for phi and random effects
                iter = 6000, warmup = 1000, chains = 3)
prior_summary(mvSevD2Mod)
summary(mvSevD2Mod)
# no significant effects of biomass
pp_check(mvSevD2Mod, nsamples = 100)

#### start here : fit above model with severity change and normal distribution ####