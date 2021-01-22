##### info ####

# file: mv_biomass_data_processing_2018_density_exp
# author: Amy Kendig
# date last edited: 1/21/21
# goal: edit Mv seed data and check for errors
# background: data collected in Keith's lab, notes indicate confusion over sample labels


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
mv_bio <- read_csv("./data/mv_biomass_oct_2018_density_exp.csv")


#### edit data ####

# look at notes
unique(mv_bio$processing_notes)

# edit biomass values and join with other data
micro <- mv_bio %>%
  mutate(site = case_when(processing_notes == "duplicate - one is D2, probably this one because the other was sorted with other D1 bags; Quynh guessed on which seeds go with which biomass (see email)" ~ "D2", 
                          TRUE ~ site),
         treatment = case_when(processing_notes == "duplicate - one is 7F, ID = A (not dated), probably this one" ~ "fungicide",
                               TRUE ~ treatment)) 

#### save data ####
write_csv(micro, "./intermediate-data/mv_processed_biomass_oct_2018_density_exp.csv")
