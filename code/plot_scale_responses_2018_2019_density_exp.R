#### info ####

# file: plot_scale_responses_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 4/21/21
# goal: plot-level biomass, seeds, and severity


# t-test for each species in each year with plots paired (compares seedling to seedling and adult to adult for Ev)
# use all plots in which that species is planted as a background species
# logit-transform severity
# hedges d to display effect sizes


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(car)
library(tidyverse)
library(effsize)

# import data
d1dat <- read_csv("intermediate-data/plot_biomass_seeds_severity_2018_density_exp.csv")
d2dat <- read_csv("intermediate-data/plot_biomass_seeds_severity_2019_density_exp.csv")

# logit adjustment
log_adj = 0.001


#### functions ####

dat_format_fun <- function(res_var, sp_abb, year){

  # select dataset
  if(year == "2018"){
    dat <- d1dat2
  }else{
    dat <- d2dat2
  }

  # filter by species
  dat2 <- dat %>%
    filter(sp == sp_abb)

  # make data wide
  datw <- dat2 %>%
    select(site, plot, treatment, {{res_var}}) %>%
    pivot_wider(names_from = treatment,
                values_from = {{res_var}}) %>%
    drop_na()

  return(datw)

}

t_test_fun <- function(res_var, sp_abb, year){
  
  # format data
  datw <- dat_format_fun(res_var, sp_abb, year)

  # t-test
  mod <- t.test(datw$fungicide, datw$water, paired = T)

  # output
  dato <- tibble(t_val = as.numeric(mod$statistic),
                 t_df = as.numeric(mod$parameter),
                 t_p = mod$p.value)
  return(dato)
  
}

hedges_d_fun <- function(res_var, sp_abb, year){

  # format data
  datw <- dat_format_fun(res_var, sp_abb, year)
  
  # hedge's d
  d <- cohen.d(datw$fungicide, datw$water, hedges.correction = T, paired = T)
  
  # output
  dato <- tibble(hedges_d = d$estimate,
                 d_lower = as.numeric(d$conf.int[1]),
                 d_upper = as.numeric(d$conf.int[2]))
  return(dato)

}


#### edit data ####

# logit-transform proportions
d1dat2 <- d1dat %>%
  mutate(logit_jul_severity = logit(jul_severity, adjust = log_adj),
         logit_late_aug_severity = logit(late_aug_severity, adjust = log_adj),
         logit_sep_severity = logit(sep_severity, adjust = log_adj)) %>%
  select(-c(jul_severity, late_aug_severity, sep_severity))

d2dat2 <- d2dat %>%
  mutate(logit_may_severity = logit(may_severity, adjust = log_adj),
         logit_jun_severity = logit(jun_severity, adjust = log_adj),
         logit_jul_severity = logit(jul_severity, adjust = log_adj),
         logit_early_aug_severity = logit(early_aug_severity, adjust = log_adj),
         logit_late_aug_severity = logit(late_aug_severity, adjust = log_adj)) %>%
  select(-c(may_severity, jun_severity, jul_severity, early_aug_severity, late_aug_severity))


#### stats ####

# 2018
d1Sum <- d1dat2 %>%
  pivot_longer(cols = c(biomass.g_m2, seeds, logit_jul_severity:logit_sep_severity),
               names_to = "response",
               values_to = "values") %>%
  filter(!(sp == "Ev" & response == "biomass.g_m2")) %>%
  select(sp, response) %>%
  unique() %>%
  mutate(hedges_d = map2(response, as.character(sp), hedges_d_fun, year = "2018"),
         t_test = map2(response, as.character(sp), t_test_fun, year = "2018")) %>%
  unnest_wider(hedges_d) %>%
  unnest_wider(t_test)
  
# 2019
d2Sum <- d2dat2 %>%
  pivot_longer(cols = c(seeds, biomass.g_m2:logit_late_aug_severity),
               names_to = "response",
               values_to = "values") %>%
  filter(!(sp == "Mv" & response == "logit_may_severity")) %>%
  select(sp, response) %>%
  unique() %>%
  mutate(hedges_d = map2(response, as.character(sp), hedges_d_fun, year = "2019"),
         t_test = map2(response, as.character(sp), t_test_fun, year = "2019")) %>%
  unnest_wider(hedges_d) %>%
  unnest_wider(t_test)


#### format for figure ####

# combine data
# indicator of significance
# rename response
# remove unnecessary responses



temp <- dat_format_fun(biomass.g_m2, "Ev", "2019")
temp1 <- t.test(temp$fungicide, temp$water, paired = T)
temp2 <- cohen.d(temp$fungicide, temp$water, hedges.correction = T, paired = T)