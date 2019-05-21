##### info ####

# file: covariate-data-processing-2018
# author: Amy Kendig
# date last edited: 5/3/19
# goal: create a dataset of covariates


#### set up ####

# clear all existing data
#rm(list=ls())

# load packages
#library(tidyverse)
#library(GGally)

# import data
sj <- read_csv("./data/soil-moisture-jun-2018-density-exp.csv") # June soil moisture
so <- read_csv("./data/soil-moisture-oct-2018-density-exp.csv") # October soil moisture
cc <- read_csv("./data/canopy-cover-jul-2018-density-exp.csv") # canopy cover
bj <- read_csv("./data/plot-edge-mv-weight-jul-2018-density-exp.csv") # July biomass and infection
be <- read_csv("./data/plot-edge-mv-weight-early-aug-2018-density-exp.csv") # early August biomass
bl <- read_csv("./data/plot-edge-mv-weight-late-aug-2018-density-exp.csv") # late August biomass
bs <- read_csv("./data/plot-edge-mv-weight-sep-2018-density-exp.csv") # September biomass and infection


#### edit data ####

sj <- sj %>%
  select(-c(date, soil_moisture.period)) %>%
  rename(soil_moisture_jun = soil_moisture.vwc)

so <- so %>%
  rowwise() %>%
  mutate(soil_moisture_oct = mean(c(soil_moisture.vwc.1, soil_moisture.vwc.2, soil_moisture.vwc.3))) %>%
  select(site, plot, treatment, soil_moisture_oct)

unique(cc$processing_notes)
cc <- cc %>%
  select(-c(count, type, processor, process_date, processing_notes)) %>%
  rename(canopy_cover = percent)
  
unique(bj$process_notes)
bj <- bj %>%
  select(site, plot, treatment, mv.g, mv_inf.g) %>%
  mutate(mv_jul.g = mv.g + mv_inf.g,
         mv_inf_jul.prop = mv_inf.g / mv_jul.g,
         mv_inf_jul.g = mv_inf.g) %>%
  select(-c(mv.g, mv_inf.g))

unique(be$process_notes)
be <- be %>%
  select(site, plot, treatment, mv.g) %>%
  rename(mv_eau.g = mv.g)

unique(bl$process_notes)
bl <- bl %>%
  select(site, plot, treatment, mv.g) %>%
  rename(mv_lau.g = mv.g)

unique(bs$process_notes)
bs <- bs %>%
  select(site, plot, treatment, mv.g, mv_inf.g) %>%
  mutate(mv_sep.g = mv.g + mv_inf.g,
         mv_inf_sep.prop = mv_inf.g / mv_sep.g,
         mv_inf_sep.g = mv_inf.g) %>%
  select(-c(mv.g, mv_inf.g))

# combine data
co <- full_join(sj, so) %>%
  full_join(cc) %>%
  full_join(bj) %>%
  full_join(be) %>%
  full_join(bl) %>%
  full_join(bs)


#### visualize data ####

#co %>%
#  select(-c(site, plot, treatment)) %>%
#  ggpairs()
# the two soil moisture measurements are correlated - pick one: the October one has a fewer extreme data and includes 3 data points
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
# these match the entries

# create a mean biomass value
co <- co %>%
  rowwise() %>%
  mutate(mv_bio.g = mean(c(mv_eau.g, mv_sep.g, na.rm = T)))

# revised covariates
#co %>%
#  select(soil_moisture_oct, canopy_cover, mv_bio.g, mv_inf_sep.prop) %>%
#  ggpairs()

# save data
covar <- co %>%
  select(site, plot, treatment, soil_moisture_oct, canopy_cover, mv_bio.g, mv_inf_sep.prop)
