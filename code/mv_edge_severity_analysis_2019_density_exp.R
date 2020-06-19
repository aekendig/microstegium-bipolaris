##### info ####

# file: mv_edge_severity_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 6/17/20
# goal: analyze background disease severity


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(glmmTMB)

# import all raw data files
dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")


#### edit data ####

# add columns
dat2 <- dat %>%
  mutate(severity = lesion_area.pix / leaf_area.pix,
         severity_adjusted = case_when(severity == 0 ~ 0.001,
                                       severity >= 1 ~ 0.999,
                                             TRUE ~ severity),
         Month = recode(month, early_aug = "Early August", jul = "July", 
                        jun = "June", late_aug = "Late August", may = "May", sep = "September") %>%
           fct_relevel("May", "June", "July", "Early August", "Late August", "September"),
         Month_num = as.numeric(Month))


#### visualize ####

# all data
ggplot(dat2, aes(Month_num, severity, color = as.factor(plot))) +
  geom_point(aes(shape = treatment)) +
  geom_line(aes(linetype = treatment)) +
  facet_wrap(~ site, ncol = 1)
# some of the May values are higher than expected given June and July
# lesions detected on May scans are not visually obvious
# tried lightening May images and a slightly different setting for lesions, but these didn't change the lesion detection much
# two cases where every value is 1: D1 Sep and D4 Late Aug
# a lot of the scans in late August for D2-D4 are highly senesced
# will probably leave out September since there is so much missing data
# the sites look like they peak at different times

# remove D4 late August
dat2 %>%
  filter(!(Month == "Late August" & site == "D4")) %>%
  ggplot(aes(Month_num, severity, color = as.factor(plot))) +
  geom_point(aes(shape = treatment)) +
  geom_line(aes(linetype = treatment)) +
  facet_wrap(~ site, ncol = 1)

# only early August
dat2 %>%
  filter(Month == "Early August") %>%
  ggplot(aes(as.factor(plot), severity, color = treatment)) +
  stat_summary(fun = "mean", geom = "point", size = 3) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1)
# use these to account for background variation in disease pressure, which is disproportionate for the two treatments in some cases (1, 9, 10)


#### treatment analysis ####

# divide data
jul_dat <- dat2 %>% filter(Month == "July")
eau_dat <- dat2 %>% filter(Month == "Early August")

# model
jul_mod <- glmmTMB(severity_adjusted ~ treatment, data = jul_dat, family = "beta_family")
summary(jul_mod)
# no treatment effect

eau_mod <- glmmTMB(severity_adjusted ~ treatment, data = eau_dat, family = "beta_family")
summary(eau_mod)
# no treatment effect
