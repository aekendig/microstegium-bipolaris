##### info ####

# file: plot_severity_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 4/24/21
# goal: effects of within and outside of plot severity


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(betareg)

# import data
d1dat <- read_csv("intermediate-data/plot_severity_2018_density_exp.csv")
# plot_data_processing_2018_density_exp.R
d2dat <- read_csv("intermediate-data/plot_severity_2019_density_exp.csv")
# plot_data_processing_2019_density_exp.R
edgeSevD2Dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
# leaf_scans_data_processing_2019_density_exp.R
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
envD2Dat <- read_csv("intermediate-data/temp_humidity_monthly_2019_density_exp.csv")
# temp_humidity_data_processing_2019_density_exp


#### functions ####

# transform data to account for 0's and 1's
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}


#### edit data ####

# format edge severity
edgeSevD2Dat2 <- edgeSevD2Dat %>%
  filter(!(month %in% c("sep"))) %>%
  mutate(edge_severity = lesion_area.pix / leaf_area.pix,
         edge_severity = ifelse(edge_severity > 1, 1, edge_severity)) %>%
  select(month, site, plot, treatment, edge_severity)

# one plot in one month is missing, use nearby plots
# make wide by month
edgeSevD2Dat3 <- tibble(site = "D1", plot = 6, treatment = "fungicide", month = "early_aug") %>%
  mutate(edge_severity = edgeSevD2Dat2 %>% 
           filter(site == "D1" & plot %in% c(4, 5) & treatment == "water" & month == "early_aug") %>%
           summarise(sev = mean(edge_severity)) %>%
           pull(sev)) %>%
  full_join(edgeSevD2Dat2) %>%
  pivot_wider(names_from = month,
              values_from = edge_severity,
              names_glue = "{month}_edge_severity") %>%
  drop_na() # use to check for missing values


#### intraspecific density regressions ####

# function
intra_dens_fun <- function(sev, SP, AGE, edge = "no", edge_sev = NULL, dat){
  
  # select group
  # tranform column
  # add density info
  dat2 <- dat %>%
    filter(sp == SP & age == AGE & treatment == "water") %>%
    mutate(severity = transform01({{sev}})) %>%
    left_join(plotsD)
  
  # add edge severity if available
  # center edge severity so that intercept is at mean value
  if(edge == "yes"){
    dat3 <- dat2 %>%
      left_join(edgeSevD2Dat3 %>%
                  select(site, plot, treatment, {{edge_sev}}) %>%
                  mutate(edge_severity = {{edge_sev}} - mean({{edge_sev}})))
  } else{
    dat3 <- dat2
  }
  
  # select plots
  # add focals to density
  if(SP == "Mv"){
    dat4 <- dat3 %>%
      filter(plot %in% 1:4) %>%
      mutate(planted_density = background_density + 3)
  }else if(AGE == "seedling"){
    dat4 <- dat3 %>%
      filter(plot %in% c(1, 5:7)) %>%
      mutate(planted_density = background_density + 3)
  }else{
    dat4 <- dat3 %>%
      filter(plot %in% c(1, 8:10)) %>%
      mutate(planted_density = background_density + 1)
  }
  
  # remove missing data
  dat5 <- dat4 %>%
    filter(!is.na(severity))
  
  # fit regression
  if(edge == "yes"){
    mod <- betareg(severity ~ planted_density + edge_severity, data = dat5)
  } else{
    mod <- betareg(severity ~ planted_density, data = dat5)
  }
  
  # extract values
  est <- summary(mod)[[1]]$mean[2, 1]
  p <- summary(mod)[[1]]$mean[2, 4]
  n <- nrow(dat5)
  
  return(c(est, p, n))
  
}

# year 1
d1IntraFits <- d1dat %>%
  pivot_longer(cols = c(jul_severity, late_aug_severity, sep_severity),
               names_to = "month_severity",
               values_to = "severity") %>%
  select(sp, age, month_severity) %>%
  unique() %>%
  mutate(sev_mod = map2())
#### start here: the data part might not work - had trouble with this in plot_scale_responses and used year instead ####
# arguments:
# sev, SP, AGE, edge = "no", edge_sev = NULL, dat