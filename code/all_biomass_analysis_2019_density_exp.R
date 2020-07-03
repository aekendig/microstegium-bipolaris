##### info ####

# file: all_biomass_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 7/3/20
# goal: analyze total biomass


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(DHARMa)
library(MASS)
library(MuMIn)
library(tidyverse)
library(glmmTMB)
library(cowplot)
library(lubridate)
library(GGally)

# import all raw data files
dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")
edge_dat <- read_csv("intermediate-data/mv_edge_leaf_scans_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
plots_simple <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
temp_hum <- read_csv("intermediate-data/temp_humidity_daily_2019_density_exp.csv")
bg_bio <- read_csv("intermediate-data/bg_processed_biomass_2019_density_exp.csv")
mv_bio <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
ev_bio <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")