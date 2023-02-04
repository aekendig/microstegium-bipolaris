##### info ####

# file: weather_summary_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 2/4/23
# goal: summarize weather for manuscript metadata


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(lubridate)

# import data
weatherD1Dat <- read_csv("data/BONWR_experiment_weather_2018.csv")
weatherD2Dat <- read_csv("data/BONWR_experiment_weather_2019.csv")


#### edit data ####

# summarize weather data by day
weatherD1Dat2 <- weatherD1Dat %>%
  filter(!is.na(TAIRGZ)) %>% # remove daily summaries
  mutate(date_time = mdy_hm(utc_valid),
         date = as_date(date_time)) %>%
  group_by(date) %>%
  summarize(temp_avg = mean(TAIRGZ),
            precip_cum = sum(PPHRGZ, na.rm = T),
            mis_vals = sum(is.na(PPHRGZ))) %>%
  ungroup()

hist(weatherD1Dat2$mis_vals) # <= 8

weatherD2Dat2 <- weatherD2Dat %>%
  filter(!is.na(TAIRGZ)) %>% # remove daily summaries
  mutate(date_time = mdy_hm(utc_valid),
         date = as_date(date_time)) %>%
  group_by(date) %>%
  summarize(temp_avg = mean(TAIRGZ),
            precip_cum = sum(PPHRGZ, na.rm = T),
            mis_vals = sum(is.na(PPHRGZ))) %>%
  ungroup()

hist(weatherD2Dat2$mis_vals) # <= 1

# growing season summaries
mean(weatherD1Dat2$temp_avg)
sd(weatherD1Dat2$temp_avg)
mean(weatherD2Dat2$temp_avg)
sd(weatherD2Dat2$temp_avg)

sum(weatherD1Dat2$precip_cum)
sum(weatherD2Dat2$precip_cum)
