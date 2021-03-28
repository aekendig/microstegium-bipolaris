##### info ####

# file: germination_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/28/21
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

# model functions
source("code/brms_model_fitting_functions.R")

# fungicide effect function
post_pred_fun <- function(mod, dat){

  posterior_samples(mod) %>%
    transmute(water = exp(b_Intercept) / (1 + exp(b_Intercept)),
              fung = exp(b_Intercept + b_fungicide) / (1 + exp(b_Intercept + b_fungicide)),
              perc = 100 * (fung - water) / water) %>%
    pivot_longer(cols = everything(),
                 names_to = "parameter",
                 values_to = "effect") %>%
    group_by(parameter) %>% mean_hdi(effect)

}

post_pred_fun2 <- function(mod, dat){
  
  posterior_samples(mod) %>%
    rename(b_yearf2019_fung = "b_fungicide:yearf2019") %>%
    transmute(water18 = exp(b_Intercept) / (1 + exp(b_Intercept)),
              fung18 = exp(b_Intercept + b_fungicide) / (1 + exp(b_Intercept + b_fungicide)),
              perc18 = 100 * (fung18 - water18) / water18,
              water19 = exp(b_Intercept + b_yearf2019) / (1 + exp(b_Intercept + b_yearf2019)),
              fung19 = exp(b_Intercept + b_fungicide + b_yearf2019 + b_yearf2019_fung) / (1 + exp(b_Intercept + b_fungicide + b_yearf2019 + b_yearf2019_fung)),
              perc19 = 100 * (fung19 - water19) / water19) %>%
    pivot_longer(cols = everything(),
                 names_to = "parameter",
                 values_to = "effect") %>%
    group_by(parameter) %>% mean_hdi(effect)
  
}


#### edit data ####

# Mv data
# average across trials
mvGermD1Dat <- mvGermD1Dat1 %>%
  mutate(germination_final = ifelse(is.na(germination_final), germination_check_1, germination_final)) %>%
  select(site_plot, trial, seeds, germination_final) %>%
  full_join(mvGermD1Dat2 %>%
              select(site_plot, trial, seeds, germination_final)) %>%
  mutate(site = gsub(" .*$", "", site_plot),
         plot = gsub(".* ","", site_plot) %>% 
           gsub("[^[:digit:]]", "", .) %>% 
           as.numeric(),
         treatment = gsub(".* ","", site_plot) %>% 
           gsub("[^[:alpha:]]", "", .) %>% 
           as.factor() %>%
           recode("F" = "fungicide", "W" = "water"),
         site = ifelse(site == "P1", "D1", site),
         fungicide = ifelse(treatment == "fungicide", 1, 0)) %>%
  group_by(site, plot, treatment, fungicide) %>%
  summarise(germinants = mean(germination_final) %>% round(),
            seeds_planted = mean(seeds) %>% round(),
            germination = germinants/seeds_planted) %>%
  ungroup()

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
         plotf = paste(site, plot, substr(treatment, 1, 1), sep = ""),
         yearf = as.factor(year))

# sample size
evGermDat2 %>%
  group_by(year, age) %>%
  count()


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

