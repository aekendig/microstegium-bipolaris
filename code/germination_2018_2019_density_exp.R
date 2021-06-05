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
library(tidybayes) # for mean_hdi

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
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_rev()) %>%
  left_join(sevD1Dat2 %>%
              filter(sp == "Ev") %>%
              mutate(year = 2018) %>%
              full_join(sevD2Dat2 %>%
                          filter(sp == "Ev") %>%
                          mutate(year = 2019)))
  
# sample size
evGermDat2 %>%
  group_by(year, age) %>%
  count()

filter(evGermDat2, is.na(jul_severity)) %>% data.frame() # 3 plots
filter(evGermDat2, is.na(late_aug_severity)) %>% data.frame() # 27 missing
# use jul_severity as predictor
# or split by year and look at completeness of dataset, this reduces sample size though...

#### start here ####


#### 2018 Mv model ####

# initial visualization
ggplot(mvGermD1Dat, aes(treatment, germination, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

# model
mvGermD1Mod <- brm(data = mvGermD1Dat, family = binomial,
                   germinants | trials(seeds_planted) ~ fungicide + (1|site),
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3,
                   control = list(adapt_delta = 0.99999, max_treedepth = 15))
mod_check_fun(mvGermD1Mod)

# posterior means
post_pred_fun(mvGermD1Mod)

# save
save(mvGermD1Mod, file = "output/mv_germination_model_2018_density_exp.rda")


#### Ev model ####

# initial visualization
ggplot(evGermDat2, aes(treatment, germination, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ yearf)

# model
evGermMod <- brm(data = evGermDat2, family = binomial,
                   germinants | trials(seeds_planted) ~ fungicide * yearf + (1|site),
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3,
                 control = list(adapt_delta = 0.99999, max_treedepth = 15))
mod_check_fun(evGermMod)

# posterior means
post_pred_fun2(evGermMod)

# save
save(evGermMod, file = "output/ev_germination_model_2018_2019_density_exp.rda")

