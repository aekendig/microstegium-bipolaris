##### info ####

# file: covariate_data_processing_2018_density_exp
# author: Amy Kendig
# date last edited: 3/6/20
# goal: create a dataset of covariates for the density experiment


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(GGally)

# import data
sj <- read_csv("./data/soil_moisture_jun_2018_density_exp.csv") # June soil moisture
so <- read_csv("./data/soil_moisture_oct_2018_density_exp.csv") # October soil moisture
cc <- read_csv("./data/canopy_cover_jul_2018_density_exp.csv") # canopy cover
bj <- read_csv("./data/plot_edge_mv_weight_jul_2018_density_exp.csv") # July biomass and infection
be <- read_csv("./data/plot_edge_mv_weight_early_aug_2018_density_exp.csv") # early August biomass
bl <- read_csv("./data/plot_edge_mv_weight_late_aug_2018_density_exp.csv") # late August biomass
bs <- read_csv("./data/plot_edge_mv_weight_sep_2018_density_exp.csv") # September biomass and infection


#### edit data ####

# June soil moisture
sj2 <- sj %>%
  mutate(soil_moisture_jun.prop = soil_moisture.vwc/100) %>%
  select(site, plot, treatment, soil_moisture_jun.prop)

# October soil moisture
so2 <- so %>%
  rowwise() %>%
  mutate(soil_moisture_oct.prop = mean(c(soil_moisture.vwc.1, soil_moisture.vwc.2, soil_moisture.vwc.3))/100) %>%
  select(site, plot, treatment, soil_moisture_oct.prop)

# Canopy cover
unique(cc$processing_notes)
cc2 <- cc %>%
  mutate(canopy_cover.prop = case_when(type == "o" ~ 1 - ((count * 1.04) / 100),
                                     type == "c" ~ (count * 1.04) / 100)) %>%
  select(site, plot, treatment, canopy_cover.prop)

# July plot edge biomass  
unique(bj$process_notes)
hist(bj$mv.g)
filter(bj, mv.g > 8) # D1 4F, D3 8F, Penny
hist(bj$mv_inf.g)
filter(bj, mv_inf.g > 1) # D1 4F, D3 8F, Penny

bj2 <- bj %>%
  select(site, plot, treatment, mv.g, mv_inf.g) %>%
  mutate(mv_jul.g = mv.g + mv_inf.g,
         mv_inf_jul.prop = mv_inf.g / mv_jul.g,
         mv_inf_jul.g = mv_inf.g) %>%
  select(site, plot, treatment, mv_jul.g, mv_inf_jul.prop, mv_inf_jul.g)

# Early August plot edge biomass
unique(be$process_notes)
hist(be$mv.g)
filter(be, mv.g > 15) # D1 4F, D2 4W

be2 <- be %>%
  select(site, plot, treatment, mv.g) %>%
  rename(mv_eau.g = mv.g)

# late August plot edge biomass
unique(bl$process_notes)
hist(bl$mv.g)

bl2 <- bl %>%
  select(site, plot, treatment, mv.g) %>%
  rename(mv_lau.g = mv.g)

# September plot edge biomass
unique(bs$process_notes)
hist(bs$mv.g)
filter(bs, mv.g > 12) # D4 7W, Penny
hist(bs$mv_inf.g)
filter(bs, mv_inf.g > 2.5) # D4 7W, Penny

bs2 <- bs %>%
  select(site, plot, treatment, mv.g, mv_inf.g) %>%
  mutate(mv_sep.g = mv.g + mv_inf.g,
         mv_inf_sep.prop = mv_inf.g / mv_sep.g,
         mv_inf_sep.g = mv_inf.g) %>%
  select(-c(mv.g, mv_inf.g))


# combine data
co <- full_join(sj2, so2) %>%
  full_join(cc2) %>%
  full_join(bj2) %>%
  full_join(be2) %>%
  full_join(bl2) %>%
  full_join(bs2)


#### visualize data ####

#co %>%
#  select(-c(site, plot, treatment)) %>%
#  ggpairs()
# the two soil moisture measurements are correlated - pick one: the October one has fewer extreme data and includes 3 data points
# July biomass correlated with early and late August
# early and late August biomass are correlated
# all biomass measurements have some extreme (high) values
# infected proportion in July and September are only correlated 0.3, the September data has fewer extreme values

# examine biomass values more
sum(is.na(co$mv_jul.g)) #3
sum(is.na(co$mv_eau.g)) #2
sum(is.na(co$mv_lau.g)) #3
sum(is.na(co$mv_sep.g)) #2
# Mv early August and September are uncorrelated and have all values

filter(co, mv_eau.g > 15) %>% 
  select(site, plot, treatment, mv_eau.g)
filter(co, mv_sep.g > 15) %>% 
  select(site, plot, treatment, mv_sep.g)
# these match the hand-written entries

co2 <- co %>%
  select(site:canopy_cover.prop, mv_eau.g, mv_sep.g, mv_inf_jul.prop, mv_inf_sep.prop)


#### output intermediate data ####

write_csv(co2, "intermediate-data/covariates_2018_density_exp.csv")

