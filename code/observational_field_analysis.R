##### info ####

# file: observational_field_analysis
# author: Amy Kendig
# date last edited: 5/5/20
# goal: evaluate the effects of Bipolaris infection on proportion Mv biomass


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
dat <- read_csv("./data/Mv disease survey data simple - pop to ship - 011618.csv")


#### edit data ####
dat2 <- dat %>%
  rename(TotBiomass = 'Total biomass') %>%
  mutate(prop_mv = TotMvBiomass / TotBiomass)


#### visualize ####
ggplot(dat2 %>%
         filter(!is.na(Bipolaris)), 
       aes(x = as.factor(Bipolaris), y = prop_mv)) +
  stat_summary(geom = "point", fun = "mean", size = 4) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0.2) +
  scale_x_discrete(labels = c("0" = "No Bipolaris",
                              "1" = "Bipolaris")) +
  ylab(expression(paste("Proportion biomass ", italic(Microstegium), sep = ""))) +
  theme(axis.title.x = element_blank())
  
