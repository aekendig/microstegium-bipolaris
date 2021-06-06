##### info ####

# file: germination_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 6/5/21
# goal: analyses of germination


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(GGally)

# import data
mvGermD1Dat1 <- read_csv("data/mv_germination_disease_set_1_2018_density_exp.csv")
mvGermD1Dat2 <- read_csv("data/mv_germination_disease_set_2_2018_density_exp.csv")
evGermDat <- read_csv("data/ev_germination_2018_2019_density_exp.csv")
sevD1Dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
sevD2Dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R

# model functions
source("code/brms_model_fitting_functions.R")


#### edit data ####

# severity
sevD1Dat2 <- sevD1Dat %>%
  select(-lesions) %>%
  pivot_wider(names_from = month,
              values_from = severity,
              names_glue = "{month}_severity") %>%
  mutate(treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev())

sevD2Dat2 <- sevD2Dat %>%
  select(-lesions) %>%
  pivot_wider(names_from = month,
              values_from = severity,
              names_glue = "{month}_severity") %>%
  mutate(treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev())

# notes
unique(mvGermD1Dat1$notes_check_1) # seedlings with lesions in notes
unique(mvGermD1Dat1$notes_check_2) # may be a contaminate on plate
unique(mvGermD1Dat1$notes_check_3)
unique(mvGermD1Dat1$notes_germination_check_1) # some plates were put into fridge during one day
unique(mvGermD1Dat1$notes_germination_final) 
unique(mvGermD1Dat2$notes) # may be a contaminate on plate

# Mv data
# average across trials
mvGermD1Dat <- mvGermD1Dat1 %>%
  mutate(germination_final = ifelse(is.na(germination_final), germination_check_1, germination_final),
         seeds_dark_check_1 = rowSums(cbind(seeds_dark_check_1, seeds_pink_check_1, seeds_red_check_1, seeds_green_check_1), na.rm = T),
         seeds_dark_check_2 = rowSums(cbind(seeds_dark_check_2, seeds_pink_check_2, seeds_red_check_2, seeds_green_check_2), na.rm = T),
         seeds_seeds_dark_check_3 = rowSums(cbind(seeds_dark_check_3, seeds_red_check_3, seeds_green_check_3), na.rm = T),
         seeds_dark = pmax(seeds_dark_check_1, seeds_dark_check_2, seeds_dark_check_3, na.rm = T),
         seeds_light = pmax(seeds_light_check_1, seeds_light_check_2, seeds_light_check_3, na.rm = T)) %>%
  select(site_plot, trial, seeds, germination_final, seeds_dark, seeds_light) %>%
  full_join(mvGermD1Dat2 %>%
              mutate(seeds_dark = pmax(seeds_dark_check_1, seeds_dark_check_2, na.rm = T),
                     seeds_light = pmax(seeds_light_check_1, seeds_light_check_2, na.rm = T)) %>%
              select(site_plot, trial, seeds, germination_final, seeds_dark, seeds_light)) %>%
  mutate(seeds_infect = rowSums(cbind(seeds_dark, seeds_light), na.rm = T),
         prop_germ = germination_final / seeds,
         prop_dark = seeds_dark / seeds,
         prop_light = seeds_light / seeds,
         prop_infect = seeds_infect / seeds,
         site = gsub(" .*$", "", site_plot),
         plot = gsub(".* ","", site_plot) %>% 
           gsub("[^[:digit:]]", "", .) %>% 
           as.numeric(),
         treatment = gsub(".* ","", site_plot) %>% 
           gsub("[^[:alpha:]]", "", .) %>% 
           as.factor() %>%
           recode("F" = "fungicide", "W" = "control (water)") %>%
           fct_rev(),
         site = ifelse(site == "P1", "D1", site),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1))) %>%
  left_join(sevD1Dat2 %>%
              filter(sp == "Mv"))

# check
filter(mvGermD1Dat, prop_germ > 1 | prop_dark > 1 | prop_light > 1 | prop_infect > 1) %>%
  data.frame()
# don't use prop_infect
filter(mvGermD1Dat, is.na(jul_severity)) # 1 plot
filter(mvGermD1Dat, is.na(late_aug_severity)) # 3 plots
filter(mvGermD1Dat, is.na(sep_severity)) # 3 plots (one missing late aug)

# correct reduction in emergents between week 3 and 4
# correct the increase in cut-tops
# correct repair of cut tops between weeks 3 and 4
evGermDat2 <- evGermDat %>%
  filter(seeds_planted > 0) %>%
  mutate(week_4_emerg = case_when(week_4_emerg < week_3_emerg ~ week_3_emerg,
                                  TRUE ~ week_4_emerg),
         week_3_cut_tops = case_when(week_4_cut_tops > week_3_cut_tops & week_4_cut_tops <= week_2_emerg ~ week_4_cut_tops,
                                     TRUE ~ week_3_cut_tops),
         week_4_cut_tops = case_when(week_4_cut_tops > week_2_emerg ~ week_3_cut_tops,
                                     week_4_cut_tops < week_3_cut_tops ~ week_3_cut_tops,
                                     TRUE ~ week_4_cut_tops),
         week_3_new_emerg = week_3_emerg - week_3_cut_tops,
         week_4_new_emerg = week_4_emerg - week_4_cut_tops,
         germinants = week_2_emerg + week_4_new_emerg + week_4_soil_germ,
         germination = germinants/seeds_planted,
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         yearf = ifelse(year == 2018, 0, 1),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev()) %>%
  left_join(sevD1Dat2 %>%
              filter(sp == "Ev") %>%
              mutate(year = 2018) %>%
              full_join(sevD2Dat2 %>%
                          filter(sp == "Ev") %>%
                          mutate(year = 2019))) %>%
  mutate(final_severity = case_when(year == 2018 ~ sep_severity,
                                    year == 2019 ~ late_aug_severity))
  
# sample size
evGermDat2 %>%
  group_by(year, age) %>%
  count()
# 40 in 2018
# 56 in 2019

# number missing?
filter(evGermDat2, is.na(jul_severity)) %>% data.frame() # 3 all 2019
filter(evGermDat2, is.na(late_aug_severity)) %>% data.frame() # 27
filter(evGermDat2, is.na(final_severity)) %>% data.frame() # 31
filter(evGermDat2, year == 2018 & is.na(late_aug_severity)) %>% data.frame() # 12
filter(evGermDat2, year == 2018 & is.na(sep_severity)) %>% data.frame() # 16
filter(evGermDat2, year == 2019 & is.na(early_aug_severity)) %>% data.frame() # 7
filter(evGermDat2, year == 2019 & is.na(late_aug_severity)) %>% data.frame() # 15

# split by year
evGermD1Dat <- evGermDat2 %>% filter(year == 2018)
evGermD2Dat <- evGermDat2 %>% filter(year == 2019)


#### Mv models ####

# initial visualization
mvGermD1Dat %>%
  select(jul_severity, late_aug_severity, sep_severity,
         prop_dark, prop_light, prop_germ) %>%
  ggpairs()

# model
mvGermD1Mod <- brm(data = mvGermD1Dat, family = binomial,
                   germination_final | trials(seeds) ~ prop_dark + prop_light + (1|plotf),
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvGermD1Mod)
# prop dark decreases germination
# prop light increases germination

mvPropDarkMod <- brm(data = mvGermD1Dat, family = binomial,
                     seeds_dark | trials(seeds) ~ sep_severity + (1|plotf),
                     prior <- c(prior(normal(0, 10), class = "Intercept"),
                                prior(normal(0, 10), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvPropDarkMod)
# severity increases prop dark

mvPropLightMod <- brm(data = mvGermD1Dat, family = binomial,
                     seeds_light | trials(seeds) ~ sep_severity + (1|plotf),
                     prior <- c(prior(normal(0, 10), class = "Intercept"),
                                prior(normal(0, 10), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvPropLightMod)
# severity increases prop light

# save
save(mvGermD1Mod, file = "output/mv_germination_infection_model_2018_density_exp.rda")
save(mvPropDarkMod, file = "output/mv_seed_infection_dark_model_2018_density_exp.rda")
save(mvPropLightMod, file = "output/mv_seed_infection_light_model_2018_density_exp.rda")


#### Ev model ####

# initial visualization
evGermD1Dat %>%
  select(jul_severity, late_aug_severity, sep_severity, germination) %>%
  ggpairs()

evGermD2Dat %>%
  select(jul_severity, early_aug_severity, late_aug_severity, germination) %>%
  ggpairs()

evGermDat2 %>%
  select(late_aug_severity, final_severity, germination) %>%
  ggpairs()

# model
evGermMod <- brm(data = evGermDat2, family = binomial,
                   germinants | trials(seeds_planted) ~ final_severity*yearf*age + (1|site),
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3)
mod_check_fun(evGermMod)

# hypotheses
adult18 = "final_severity = 0"
seedling18 = "final_severity + ageseedling + final_severity:ageseedling = 0"
adult19 = "final_severity + yearf + final_severity:yearf = 0"
seedling19 = "final_severity + yearf + ageseedling  + final_severity:yearf + final_severity:ageseedling + yearf:ageseedling + final_severity:yearf:ageseedling = 0"

hypothesis(evGermMod, c(adult18, seedling18, adult19, seedling19))
# only the first is negative

# save
save(evGermMod, file = "output/ev_germination_severity_model_2018_2019_density_exp.rda")

