##### info ####

# file: ev_biomass_2019_fungicide_exp
# author: Amy Kendig
# date last edited: 5/12/21
# goal: fungicide effects on Ev biomass


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
dat <- read_csv("data/ev_biomass_dec_2019_fungicide_exp.csv")


#### edit data ####

# notes
unique(dat$notes)

# remove one with seeds
# make wide by biomass type
dat2 <- dat %>%
  filter(is.na(notes)) %>%
  pivot_wider(names_from = biomass_type,
              values_from = weight.g) %>%
  mutate(total_biomass = live + dead,
         prop_dead = dead / total_biomass)

#### stats ####

t.test(dead ~ treatment, data = dat2)
# no difference

t.test(live ~ treatment, data = dat2)
# no difference

t.test(total_biomass ~ treatment, data = dat2)
# no difference


#### litter ####

mean(dat2$prop_dead, na.rm = T)
sd(dat2$prop_dead, na.rm = T)
